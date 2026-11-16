"""Private IGRF magnetic-cap v5 state adapter around frozen original frame/CRC/CT readers.
Numerical operations run only on allocated compute, never local/login.
"""
from pathlib import Path
import numpy as np,hashlib,json,sys
MI_NJ,MI_NK=48,64
MI_SHAPE=(49,35,10,2)
def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()


def array(f,dtype,count,label):
    a=np.fromfile(f,dtype,int(count));require(a.size==count,'truncated '+label)
    return a


def execution_layout(manifest):
    """Explicit new layouts; only the historical12-rank default is inferred."""
    cube=manifest.get('solver','cube_LFM')=='cube_LFM';grid=list(map(int,manifest['grid'][:3]))
    declared=manifest.get('execution_layout')
    if declared is None:
        ranks=manifest.get('physical',{}).get('mpi_ranks',12)
        require(ranks==12,'nonhistorical MPI count requires frozen execution_layout')
        dims=[2,1,1] if cube else [2,3,2]
    else:
        ranks=declared['mpi_count'];dims=declared['dims']
        require(declared.get('global_shape')==grid,'execution_layout/global grid mismatch')
    require(isinstance(ranks,int) and not isinstance(ranks,bool) and 0<ranks<=192,'explicit supported MPI count must be in1..192')
    require(len(dims)==3 and all(isinstance(d,int) and not isinstance(d,bool) and d>0 for d in dims),'three positive integer partition dimensions')
    require(all(g>0 for g in grid),'positive global cell counts')
    require(manifest.get('physical',{}).get('mpi_ranks',ranks)==ranks,'physical/layout MPI count mismatch')
    require((6 if cube else 1)*int(np.prod(dims))==ranks,'partition/MPI count mismatch')
    if cube:
        require(len(manifest['grid'])==4 and manifest['grid'][3]==6 and grid[1]==grid[2] and dims[1]==dims[2], 'six-patch square angular topology')
        require(all(g%d==0 for g,d in zip(grid,dims)),'cube partition must divide global cells')
    return dict(mpi_count=ranks,dims=list(dims),global_shape=grid,explicit=declared is not None)


def rank_tile(rank,layout,cube=True):
    dims=layout['dims'];grid=layout['global_shape'];ranks=layout['mpi_count']
    require(0<=rank<ranks,'rank outside frozen execution layout')
    perface=int(np.prod(dims));face=rank//perface if cube else None;q=rank%perface
    k=q%dims[2];q//=dims[2];j=q%dims[1];i=q//dims[1]
    coords=[i,j,k]
    local=[g//d+int(c<g%d) for g,d,c in zip(grid,dims,coords)]
    first=[c*(g//d)+min(c,g%d) for g,d,c in zip(grid,dims,coords)]
    return face,first,local


def crc64_bytes(raw):
    """Independent slicing-by-eight CRC64-ECMA; wire bytes, no model state."""
    import struct
    poly=0x42F0E1EBA9EA3693;mask=(1<<64)-1
    tables=[];first=[]
    for b in range(256):
        c=b<<56
        for _ in range(8):c=((c<<1)^poly if c>>63 else c<<1)&mask
        first.append(c)
    tables.append(first)
    for _ in range(7):tables.append([((c<<8)&mask)^first[c>>56] for c in tables[-1]])
    t0,t1,t2,t3,t4,t5,t6,t7=tables;c=0;n=len(raw)//8*8
    for (w,) in struct.iter_unpack('>Q',memoryview(raw)[:n]):
        x=c^w;c=t7[x>>56]^t6[(x>>48)&255]^t5[(x>>40)&255]^t4[(x>>32)&255]^t3[(x>>24)&255]^t2[(x>>16)&255]^t1[(x>>8)&255]^t0[x&255]
    for b in memoryview(raw)[n:]:c=((c<<8)&mask)^first[(c>>56)^b]
    return c if c<(1<<63) else c-(1<<64)


def validate_mi_spatial(a):
    require(a.shape==MI_SHAPE and np.isfinite(a).all(),'MI spatial schema/finite state')
    require(np.array_equal(a[:,:,:2,0],a[:,:,:2,1]),'MI hemisphere coordinate axes')
    require(np.array_equal(a[:,:,4,:],a[:,:,8,:]+a[:,:,9,:]),'MI total/FAC/corotation potential closure')
    require(np.max(abs(a[:,:,6,:]-5.))<=1e-12,'MI SigmaP default5S')
    require(np.max(abs(a[:,:,7,:]))<=1e-12,'MI Hall default0S')
    require(np.all(a[:,:,9,:]==0),'rotation EMF is separate; scalar corotation potential remains zero')
    theta=np.sqrt(a[:,:,0,:]**2+a[:,:,1,:]**2)
    lat=(90.-180.*theta/np.pi)*np.array([1.,-1.])
    mlt=np.mod(12.+12.*np.arctan2(a[:,:,1,:],a[:,:,0,:])/np.pi,24.)
    require(np.allclose(a[:,:,2,:],lat,rtol=0,atol=5e-12),'MI signed magnetic latitude schema')
    require(np.all((a[:,:,3,:]>=0)&(a[:,:,3,:]<24)),'MI MLT range')
    require(np.max(abs(np.mod(a[:,:,3,:]-mlt+12,24)-12))<=5e-12,'MI MLT/native axes schema')
    return dict(CPCP_kV=(np.ptp(a[:,:,8,:],axis=(0,1))*1e-3).tolist(),
                mapped_solver_input_FAC_max_A_m2=float(np.max(abs(a[:,:,5,:]))))


def read_mi_snapshot(path):
    with path.open('rb') as f:
        require(f.read(8)==b'LFMCUMI1','MI snapshot magic')
        require(array(f,'<i4',4,'MI shape').tolist()==list(MI_SHAPE),'MI snapshot shape')
        step=int(array(f,'<i4',1,'MI step')[0]);time,last,rion,sigma=array(f,'<f8',4,'MI snapshot clock/config')
        require(np.isfinite([time,last,rion,sigma]).all() and 0<=last<=time and rion==6500000. and sigma==5.,'MI snapshot epoch/default controls')
        spatial=array(f,'<f8',int(np.prod(MI_SHAPE)),'MI snapshot').reshape(MI_SHAPE,order='F')
        require(not f.read(),'unexpected MI snapshot trailing bytes')
    return dict(step=step,time=float(time),last_time=float(last),spatial=spatial,
                spatial_bytes=spatial.tobytes(),**validate_mi_spatial(spatial))


def csv(path,required):
    rows=np.atleast_1d(np.genfromtxt(path,delimiter=',',names=True))
    require(rows.size>0 and set(required)<=set(rows.dtype.names or ()),'CSV columns/rows: '+str(path))
    require(np.isfinite(np.array(rows.tolist())).all(),'finite CSV: '+str(path))
    return rows


def solve_indices(clock,cadence=10.):
    """Cold initial solve, then first natural MHD step crossing each MI event."""
    require(clock['step'][0]==0 and clock['time_s'][0]==0,'cold MI schedule begins at0')
    indices=[0];next_time=cadence
    for i,time in enumerate(clock['time_s'][1:],1):
        tolerance=64*np.finfo(float).eps*max(1.,abs(float(time)),abs(next_time))
        if time>=next_time-tolerance:
            indices.append(i)
            while time>=next_time-64*np.finfo(float).eps*max(1.,abs(float(time)),abs(next_time)):next_time+=cadence
    return np.array(indices,int)


def read_frame(path):
    with path.open('rb') as f:
        assert f.read(8)==b'LFMCUB01',path
        n=np.fromfile(f,'<i4',3);face=int(np.fromfile(f,'<i4',1)[0]);first=np.fromfile(f,'<i4',3)
        time,rp=np.fromfile(f,'<f8',2)
        assert np.isfinite([time,rp]).all() and time>=0 and rp>0
        assert np.all(n>0) and n[0]<=128 and np.all(n[1:]<=32)
        xyz=np.fromfile(f,'<f8',3*np.prod(n+1)).reshape((3,*tuple(n+1)),order='F')
        q=np.fromfile(f,'<f8',8*np.prod(n)).reshape((*tuple(n),8),order='F')
        flux=np.fromfile(f,'<f8',3*np.prod(n+1)).reshape((*tuple(n+1),3),order='F')
        assert f.read()==b'',path
    assert np.isfinite(q).all() and np.isfinite(xyz).all() and np.isfinite(flux).all()
    assert q[...,0].min()>0 and q[...,4].min()>0
    div=(flux[1:,1:,1:,0]-flux[:-1,1:,1:,0]+flux[1:,1:,1:,1]-flux[1:,:-1,1:,1]
         +flux[1:,1:,1:,2]-flux[1:,1:,:-1,2])
    # Valid initial/unswept tiles outside the tagged taper can have exactly
    # zero flux. Match the native runtime's max(tiny, local_flux_scale).
    scale=float(np.max(abs(flux)))
    norm=float(np.max(abs(div))/max(np.finfo(float).tiny,scale))
    return dict(n=n,face=face,first=first,time=time,rp=rp,xyz=xyz,q=q,flux=flux,div=norm,zero_flux=scale==0.)


def read_compact_epoch(marker,manifest):
    """Parse every shard independently; scientific output is never a restart input.

    The audit uses its already-required initial MI diagnostic only as a static
    coordinate/sign-zero oracle. The Fortran roundtrip separately proves that
    all40929 MI values reconstruct without any scientific output files.
    """
    import io,struct
    lines=marker.read_text().splitlines();fingerprint=manifest['fingerprint'];layout=execution_layout(manifest)
    require(lines[0]==fingerprint,'CPC03 manifest fingerprint')
    step,epoch,ranks=lines[1].split();step=int(step);epoch=float(epoch);ranks=int(ranks)
    magic,g,groups,mode=lines[2].split();g=int(g);groups=int(groups);mode=int(mode)
    require(magic=='LFMCPC03' and 1<=g<=128 and ranks==layout['mpi_count'] and groups==(ranks+g-1)//g,'CPC03 group policy')
    require(list(map(int,lines[3].split()))==layout['global_shape']+layout['dims'],'CPC03 grid/layout')
    count,mi_crc=map(int,lines[4].split());require(len(lines)==5+groups,'CPC03 manifest entries')
    stem=marker.with_suffix('');local=rank_tile(0,layout)[2];cells=int(np.prod(local));n=np.asarray(local)
    words=10*cells+2*((n[0]+1)*n[1]*n[2]+n[0]*(n[1]+1)*n[2]+n[0]*n[1]*(n[2]+1))
    mi=None
    if mode==3:
        cap=MI_SHAPE[0]*MI_SHAPE[1]
        require(count==82+6*cap,'IGRF magnetic-cap v5 complete held state count')
        path=Path(str(stem)+'.mi');raw=path.read_bytes();require(len(raw)==104+8*count,'CPC03 MI exact byte length')
        require(raw[:8]==b'LFMCPI03' and raw[8:72].decode()==fingerprint,'CPC03 MI schema/provenance')
        ms,mt,mn=struct.unpack_from('<qdq',raw,72)
        require((ms,mt,mn)==(step,epoch,count),'CPC03 MI epoch/count')
        state=np.frombuffer(raw,dtype='<f8',count=count,offset=96).copy();stored=struct.unpack_from('<q',raw,len(raw)-8)[0]
        require(stored==mi_crc and crc64_bytes(raw[96:-8])==stored,'CPC03 MI CRC64')
        require(np.isfinite(state).all(),'IGRF magnetic-cap v5 finite complete state')
        require(state[:4].tolist()==[5.,100.,0.001,1.],'IGRF magnetic-cap v5 exact warning policy')
        require(state[14:18].tolist()==[6380000.,6500000.,manifest['rotation']['omega_rad_s'],manifest['rotation']['phase0_rad']],'IGRF magnetic-cap v5 field/rotation identity')
        require(''.join(chr(int(c))for c in state[18:82])=='e07b15a704841612ae7f0787e365bc75da46cb969720fd9e53143b9df9bd9bcd','IGRF coefficient SHA')
        require(0<=state[4]<=step and state[4]==int(state[4]) and 0<=state[5]<=epoch,'IGRF held solve epoch')
        require(np.all(state[8:10]>=0) and np.all(state[10:12]>=0),'native SER/ITER values')
        require(np.array_equal(state[12:14],(state[8:10]<=0.001).astype(float)),'honest finite nonconvergence status')
        pn,ps,cn,cs=[state[82+i*cap:82+(i+1)*cap].reshape(MI_SHAPE[:2],order='F') for i in range(4)]
        traces=state[82+4*cap:];require(np.all(traces==np.rint(traces)),'integer retained trace statuses')
        next_time=0.;updates=0
        while epoch>=next_time-64*np.finfo(float).eps*max(1.,abs(epoch),abs(next_time)):
            next_time+=10.;updates+=1
        mi=dict(last_step=int(state[4]),last_time=float(state[5]),residual=float(state[6]),raw_FAC_max_uA_m2=float(state[7]),
          CPCP_kV=np.array([np.ptp(pn),np.ptp(ps)])*0.001,PN=pn,PS=ps,CN=cn,CS=cs,
          native_SER=state[8:10].copy(),native_ITER=state[10:12].copy(),converged=state[12:14].copy(),
          state_bytes=state.tobytes(),state_sha256=hashlib.sha256(state.tobytes()).hexdigest(),
          next_time=next_time,cadence=10.,updates=updates,independent_state=state)
    else:require(count==0 and mi_crc==0,'CPC03 MI-off packet absent')
    p=manifest['physical'];sw=p['SW'];solar=np.asarray([sw['rho_g_cm3'],*np.asarray(sw['v_km_s'])*1e5,sw['c_km_s']*1e5,*np.asarray(sw['B_nT'])*1e-5])
    result=[];records=[]
    for j,line in enumerate(lines[5:]):
        gid,first,members,length,checksum=map(int,line.split())
        require([gid,first,members]==[j,j*g,min(g,ranks-j*g)],'CPC03 shard owner range')
        path=Path(str(stem)+f'-g{gid:06}.cp');raw=path.read_bytes()
        require(len(raw)==length==172+8*words*members,'CPC03 exact shard bytes')
        require(raw[:8]==b'LFMCPS03' and raw[8:72].decode()==fingerprint,'CPC03 shard schema/provenance')
        ss,tt,dt=struct.unpack_from('<qdd',raw,72);meta=list(struct.unpack_from('<15i',raw,96));nw=struct.unpack_from('<q',raw,156)[0]
        require(ss==step and tt==epoch and np.isfinite(dt) and dt>0,'CPC03 shared clocks')
        require(meta==[mode,ranks,g,gid,first,members,*local,*layout['global_shape'],*layout['dims']] and nw==words,'CPC03 field/layout description')
        stored=struct.unpack_from('<q',raw,len(raw)-8)[0]
        require(stored==checksum and crc64_bytes(raw[80:96]+raw[164:-8])==stored,'CPC03 shard CRC64 incl clock')
        data=np.frombuffer(raw,dtype='<f8',count=words*members,offset=164).reshape((words,members),order='F')
        require(np.isfinite(data).all(),'CPC03 finite every independent state word')
        records.append(dict(file=path.name,step=step,time_s=epoch,sha256=sha(path),bytes=length))
        for k in range(members):
            rank=first+k;v=data[:,k];q=v[:5*cells].reshape((*local,5),order='F');old=v[5*cells:10*cells].reshape((*local,5),order='F')
            for a in (q,old):require(a[...,0].min()>0 and a[...,4].min()>0,'CPC03 positive current/history rho/C')
            pos=10*cells;fluxes=[]
            for _ in range(2):
                flux=np.zeros((*tuple(n+1),3),order='F')
                for axis in range(3):
                    lo=[1,1,1];lo[axis]=0;shape=[int(n[d]+1-lo[d]) for d in range(3)];amount=int(np.prod(shape))
                    sl=tuple(slice(lo[d],n[d]+1) for d in range(3))+(axis,)
                    flux[sl]=v[pos:pos+amount].reshape(shape,order='F');pos+=amount
                div=flux[1:,1:,1:,0]-flux[:-1,1:,1:,0]+flux[1:,1:,1:,1]-flux[1:,:-1,1:,1]+flux[1:,1:,1:,2]-flux[1:,1:,:-1,2]
                scale=max(np.finfo(float).tiny,float(np.max(abs(flux))))
                require(np.max(abs(div))/scale<1e-11,'CPC03 current/history CT gate')
                fluxes.append(flux)
            require(pos==words,'CPC03 semantic field count')
            result.append(dict(file=path.name,rank=rank,step=step,time_s=epoch,dt=dt,rp=p['re_cm'],
              moment=p['dipole_equatorial_g']*p['re_cm']**3,sw=solar,mi=mi,state_bytes=v.tobytes(),
              independent_q=q,independent_history=old,phi=fluxes[0],phi_history=fluxes[1]))
    return result,records
