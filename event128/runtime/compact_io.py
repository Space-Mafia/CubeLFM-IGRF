#!/usr/bin/env python3
"""CubeLFM compact HDF5 output. All array I/O runs on the allocated compute host."""
import argparse
from contextlib import contextmanager
import os
from pathlib import Path
import struct
import tempfile

import h5py
import numpy as np

HEADER = struct.Struct('<7i2d')
CHANNELS = 'rho vx vy vz sound_speed Bx_total By_total Bz_total'
FILTERS = dict(compression='gzip', compression_opts=4, shuffle=True, fletcher32=True)


@contextmanager
def new_h5(path, kind):
    """Publish only a completely closed HDF5 file, without replacing another writer."""
    fd, name = tempfile.mkstemp(dir=path.parent, prefix='.' + path.name + '.', suffix='.tmp')
    os.close(fd)
    with h5py.File(name, 'w') as out:
        out.attrs['format'] = 'CubeLFM-' + kind + '-1'
        yield out
    os.link(name, path)
    os.unlink(name)


def array(group, key, values):
    options = FILTERS if values.size >= 1024 else {}
    if values.ndim == 1 and options:
        options = dict(options, chunks=(min(8192, values.size),))
    return group.create_dataset(key, data=values, **options)


def read_n(stream, dtype, shape):
    count = int(np.prod(shape))
    values = np.fromfile(stream, dtype, count=count)
    if len(values) != count:
        raise EOFError(f'Truncated array {shape} in {stream.name}')
    return values.reshape(shape)


def write_mi(frame, path, step, epoch):
    with path.open('rb') as src:
        magic = src.read(8)
        if magic not in (b'LFMCUMI1', b'LFMCUMI2'):
            raise ValueError('Unknown MI format')
        nj, nk, nc, nh, actual_step, time, held, radius = struct.unpack('<5i3d', src.read(44))
        if actual_step != step or time != epoch:
            raise ValueError('MI and MHD epochs differ')
        g = frame.create_group('mi')
        g.attrs.update(magic=magic.decode(), step=step, time_s=time, held_s=held, radius_m=radius)
        if magic == b'LFMCUMI1':
            g.attrs['sigma'] = struct.unpack('<d', src.read(8))[0]
        array(g, 'cap', read_n(src,'<f8',(nh,nc,nk,nj)))
        if magic == b'LFMCUMI2':
            shape = struct.unpack('<4i',src.read(16))
            array(g,'aurora',read_n(src,'<f8',shape[::-1]))
            # Native edge count depends on wall dimensions, independently of cap shape.
            trace_bytes = 4*nh*nk*nj
            edge_bytes = path.stat().st_size - src.tell() - trace_bytes
            if edge_bytes <= 0 or edge_bytes % (8*nh):
                raise ValueError('Invalid edge array')
            array(g,'edge',read_n(src,'<f8',(nh,edge_bytes//(8*nh))))
            array(g,'trace',read_n(src,'<i4',(nh,nk,nj)))
        if src.read(1):
            raise ValueError('Trailing MI bytes')


def tiles(stream, ranks):
    for rank in range(ranks):
        header = HEADER.unpack(stream.read(HEADER.size))
        ni,nj,nk,face,i,j,k,t,rp = header
        yield str(rank), (ni,nj,nk), (face,i,j,k), t,rp


def publish_frame(prefix, stage, step, ranks, include_mi):
    grid_path = prefix.parent/(prefix.name+'-grid.h5')
    static = stage/'grid'
    if static.exists():
        with static.open('rb') as src, new_h5(grid_path,'grid') as grid:
            grid.attrs['ranks'] = ranks
            for rank,(ni,nj,nk),layout,t,rp in tiles(src,ranks):
                ds = array(grid,rank,read_n(src,'<f8',(nk+1,nj+1,ni+1,3)))
                ds.attrs.update(layout=layout,rp=rp)
            if src.read(1):
                raise ValueError('Trailing grid bytes')
        static.unlink()
    elif not grid_path.is_file():
        raise FileNotFoundError(grid_path)
    frame_path=prefix.parent/f'{prefix.name}-{step:07d}.h5'
    mi=stage/f'raw-mi-step{step:07d}.bin'
    with (stage/'fields').open('rb') as src, new_h5(frame_path,'frame') as out:
        out.attrs.update(step=step,ranks=ranks,grid=grid_path.name,channels=CHANNELS)
        for rank,(ni,nj,nk),layout,t,rp in tiles(src,ranks):
            group = out.create_group(rank)
            group.attrs.update(layout=layout,rp=rp)
            if rank == '0':
                out.attrs['time_s'] = t
            elif t != out.attrs['time_s']:
                raise ValueError('MHD rank clocks disagree')
            array(group,'primitive',read_n(src,'<f4',(8,nk,nj,ni)))
            array(group,'ct_flux',read_n(src,'<f8',(3,nk+1,nj+1,ni+1)))
        if src.read(1):
            raise ValueError('Trailing field bytes')
        if include_mi:
            write_mi(out,mi,step,t)
            from archive_global import write_global
            write_global(out,stage/'global-state')
    (stage/'fields').unlink()
    if include_mi:
        mi.unlink()
        (stage/'global-state').unlink()


def members(marker):
    lines = marker.splitlines()
    if len(lines)<5 or lines[2].split()[0]!=b'LFMCPC03':
        raise ValueError('Expected committed CPC03 checkpoint')
    groups=int(lines[2].split()[2])
    if len(lines)!=5+groups:
        raise ValueError('Incomplete checkpoint manifest')
    for index,line in enumerate(lines[5:]):
        group,first,ranks,length,crc=map(int,line.split())
        if group!=index:
            raise ValueError('Invalid group order')
        yield f'-g{group:06d}.cp',164,length,b'LFMCPS03'
    count=int(lines[4].split()[0])
    if count:
        yield '.mi',96,104+8*count,b'LFMCPI03'


def publish_checkpoint(prefix,stage,step):
    stem=f'raw-step{step:07d}'
    marker=stage/(stem+'.complete')
    contents=marker.read_bytes()
    paths=[marker]
    target=prefix.parent/f'{prefix.name}-step{step:07d}.h5'
    with new_h5(target,'checkpoint') as out:
        array(out,'manifest',np.frombuffer(contents,'u1'))
        for suffix,offset,length,magic in members(contents):
            path=stage/(stem+suffix)
            if path.stat().st_size!=length:
                raise ValueError(f'Truncated CPC03 member {path}')
            with path.open('rb') as src:
                header=src.read(offset)
                if header[:8]!=magic or header[8:72]!=contents.splitlines()[0]:
                    raise ValueError(f'Invalid identity for {path}')
                group=out.create_group(suffix)
                group.attrs['header']=np.void(header)
                array(group,'state',read_n(src,'<f8',((length-offset-8)//8,)))
                group.attrs['crc']=np.void(src.read(8))
            paths.append(path)
    for path in paths:
        path.unlink()


def restore_checkpoint(source,stem):
    """Recreate the exact native CPC03 records passed to CubeLFM's own CRC/physics reader."""
    with h5py.File(source,'r') as src:
        if src.attrs['format']!='CubeLFM-checkpoint-1':
            raise ValueError('Not a compact CubeLFM checkpoint')
        marker=src['manifest'][:].tobytes()
        for suffix,offset,length,magic in members(marker):
            group=src[suffix]
            path=stem.parent/(stem.name+suffix)
            with path.open('xb') as out:
                out.write(bytes(group.attrs['header']))
                for sl in group['state'].iter_chunks() if group['state'].chunks else [slice(None)]:
                    out.write(group['state'][sl].tobytes())
                out.write(bytes(group.attrs['crc']))
            if path.stat().st_size!=length:
                raise ValueError(f'Bad restored length for {path}')
        with (stem.parent/(stem.name+'.complete')).open('xb') as out:
            out.write(marker)


def read_frame_tile(frame,rank):
    """Return one tile in native Fortran axis order, coordinates in stored model units."""
    with h5py.File(frame,'r') as f, h5py.File(Path(frame).parent/f.attrs['grid'],'r') as grid:
        if f.attrs['format']!='CubeLFM-frame-1' or grid.attrs['format']!='CubeLFM-grid-1':
            raise ValueError('Unknown compact frame or grid')
        tile=f[str(rank)]
        xyz=grid[str(rank)]
        if not np.array_equal(tile.attrs['layout'],xyz.attrs['layout']) or tile.attrs['rp']!=xyz.attrs['rp']:
            raise ValueError('Frame and grid layouts disagree')
        return dict(time=f.attrs['time_s'],rp=tile.attrs['rp'],layout=tile.attrs['layout'],
                    xyz=xyz[:].T,q=tile['primitive'][:].T,flux=tile['ct_flux'][:].T)


def export_frame(frame, prefix):
    """Export a single HDF frame for existing native plotters (primitives widened from F32)."""
    with h5py.File(frame,'r') as f, h5py.File(frame.parent/f.attrs['grid'],'r') as grid:
        if f.attrs['format']!='CubeLFM-frame-1' or grid.attrs['format']!='CubeLFM-grid-1':
            raise ValueError('Expected a CubeLFM frame and matching grid')
        step=int(f.attrs['step'])
        for rank in range(int(f.attrs['ranks'])):
            tile=f[str(rank)]
            xyz=grid[str(rank)]
            if not np.array_equal(tile.attrs['layout'],xyz.attrs['layout']) or tile.attrs['rp']!=xyz.attrs['rp']:
                raise ValueError('Grid and frame rank layouts differ')
            shape=tile['primitive'].shape[:0:-1]
            path=prefix.parent/f'{prefix.name}-{step:05d}-r{rank:03d}.bin'
            with path.open('xb') as out:
                out.write(b'LFMCUB01')
                out.write(HEADER.pack(*shape,*tile.attrs['layout'],f.attrs['time_s'],tile.attrs['rp']))
                out.write(xyz[:].tobytes())
                out.write(tile['primitive'][:].astype('<f8').tobytes())
                out.write(tile['ct_flux'][:].tobytes())
        if 'mi' in f:
            export_mi(f,prefix)


def export_mi(f,prefix):
    step=int(f.attrs['step'])
    g=f['mi']
    path=prefix.parent/f'{prefix.name}-mi-step{step:07d}.bin'
    with path.open('xb') as out:
        out.write(g.attrs['magic'].encode())
        out.write(struct.pack('<5i3d',*g['cap'].shape[::-1],step,g.attrs['time_s'],g.attrs['held_s'],g.attrs['radius_m']))
        if 'sigma' in g.attrs:
            out.write(struct.pack('<d',g.attrs['sigma']))
        out.write(g['cap'][:].tobytes())
        if 'aurora' in g:
            out.write(struct.pack('<4i',*g['aurora'].shape[::-1]))
            for key in ('aurora','edge','trace'):
                out.write(g[key][:].tobytes())


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    commands=ap.add_subparsers(dest='cmd',required=True)
    for name in ('init','frame','checkpoint','restore','clean-restore','end'):
        sub=commands.add_parser(name)
        sub.add_argument('prefix',type=Path)
        if name in ('checkpoint','frame'):
            sub.add_argument('step',type=int)
        if name=='frame':
            sub.add_argument('ranks',type=int)
            sub.add_argument('mi',type=int,choices=[0,1])
        if name=='restore':
            sub.add_argument('source',type=Path)
    sub=commands.add_parser('export-checkpoint')
    sub.add_argument('source',type=Path)
    sub.add_argument('prefix',type=Path)
    sub=commands.add_parser('export-frame')
    sub.add_argument('frame',type=Path)
    sub.add_argument('prefix',type=Path)
    args=ap.parse_args()
    if args.cmd=='export-checkpoint':
        restore_checkpoint(args.source,args.prefix)
        return
    if args.cmd=='export-frame':
        export_frame(args.frame,args.prefix)
        return
    stage=args.prefix.parent/(args.prefix.name+'.io-tmp')
    if args.cmd=='init':
        stage.mkdir()
    elif args.cmd=='frame':
        publish_frame(args.prefix,stage,args.step,args.ranks,args.mi)
    elif args.cmd=='checkpoint':
        publish_checkpoint(args.prefix,stage,args.step)
    elif args.cmd=='restore':
        restore_checkpoint(args.source,stage/'restart')
    elif args.cmd=='clean-restore':
        marker=stage/'restart.complete'
        for suffix,*unused in members(marker.read_bytes()):
            (stage/('restart'+suffix)).unlink()
        marker.unlink()
    else:
        stage.rmdir()

if __name__=='__main__':
    main()
