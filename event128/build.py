"""Build the 23 October 2002 event using source files from this checkout."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
MODULE = re.compile(r'^\s*module\s+(?!procedure\b|function\b|subroutine\b)(\w+)', re.M | re.I)
USE = re.compile(r'^\s*use\s*(?:::)?\s*(\w+)', re.M | re.I)
WRAPS = (
    'cube_yy_direct_reference_evaluate', 'cube_yy_magnetic_wall_create',
    'cube_yy_magnetic_wall_apply', 'cube_yy_magnetic_wall_destroy',
    'cube_yy_igrf_configure', 'cube_yy_igrf_field', 'cube_yy_footpoint_potential',
)


def run_build(out, jobs, environment):
    manifest = json.loads((HERE / 'native/units.json').read_text())
    mkl = Path(environment['MKLROOT']).resolve()
    hdf5 = Path(environment['HDF5_ROOT']).resolve()
    mpicc, mpifc = (environment.get(k, default) for k, default in
                     [('MPICC', 'mpicc'), ('MPIFC', 'mpifort')])
    units = manifest['units']
    out.mkdir(parents=True)
    (out / 'modules').mkdir()
    rules = ['.DELETE_ON_ERROR:', '.PHONY: all', 'all: cube_event128', '']
    archives = {}
    objects = {}
    providers = {}
    sources = {}
    for unit in units:
        key = unit['object']
        if key in objects:
            raise ValueError(f'Duplicate object {key}')
        unit_dir = out / 'units' / unit['group'] / key[:-2]
        unit_dir.mkdir(parents=True)
        source = REPO / unit['source']
        local = unit_dir / source.name
        local.symlink_to(source)
        for name, relative in unit['headers'].items():
            dest = unit_dir / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.symlink_to(REPO / relative)
        obj = str((unit_dir / key).relative_to(out))
        objects[key] = obj
        sources[key] = local
        if source.suffix.lower() != '.c':
            for module in MODULE.findall(source.read_text()):
                providers[module.lower()] = obj
        if unit['group'] != 'event':
            archives.setdefault(unit['group'], []).append(obj)

    for unit in units:
        key = unit['object']
        source = sources[key]
        c = source.suffix.lower() == '.c'
        flags = manifest['profiles'][unit['profile']]
        dependencies = [] if c else sorted({providers[m.lower()] for m in USE.findall(source.read_text())
                                            if m.lower() in providers and providers[m.lower()] != objects[key]})
        command = [mpicc if c else mpifc, *flags, '-I'+str(source.parent),
                   '-I'+str(HERE / 'event-headers'), '-I'+str(hdf5 / 'include'),
                   '-I'+str(mkl / 'include')]
        if not c:
            command += ['-Jmodules', '-Imodules']
        command += ['-c', str(source.relative_to(out)), '-o', objects[key]]
        rules += [objects[key]+': '+' '.join(dependencies), '\t'+shlex.join(command), '']

    libs = []
    for group, members in archives.items():
        lib = 'lib'+group+'.a'
        rules += [lib+': '+' '.join(members), '\t'+shlex.join(['ar','rcs',lib,*members]), '']
        libs.append(lib)
    original = json.loads((HERE/'COMPILE_LINK.json').read_text())['cube_event128']
    explicit = [objects[Path(arg).name] for arg in original if arg.endswith('.o')]
    mkl_lib = next((p for p in (mkl/'lib/intel64', mkl/'lib') if (p/'libmkl_rt.so').exists()), None)
    hdf_lib = next((p for p in (hdf5/'lib', hdf5/'lib64') if (p/'libhdf5.so').exists()), None)
    if not mkl_lib or not hdf_lib:
        raise ValueError('Set MKLROOT and HDF5_ROOT to the respective development installations')
    command = [mpifc,'-pthread','-O3',*explicit,'-Wl,--start-group',*libs,'-Wl,--end-group',
               *['-Wl,--wrap='+symbol for symbol in WRAPS],
               '-L'+str(hdf_lib),'-Wl,-rpath,'+str(hdf_lib),'-lhdf5',
               '-L'+str(mkl_lib),'-Wl,-rpath,'+str(mkl_lib),'-lmkl_rt','-ldl',
               '-fopenmp','-lm','-Wl,-Map,model.map','-o','cube_event128']
    rules += ['cube_event128: '+' '.join(explicit+libs), '\t'+shlex.join(command)]
    (out/'Makefile').write_text('\n'.join(rules)+'\n')
    subprocess.run(['make','-j',str(jobs),'all'], cwd=out, check=True)
    executable = out/'cube_event128'
    print(f'Built {executable} from {len(units)} source files. '
          f'SHA-256 {hashlib.sha256(executable.read_bytes()).hexdigest()}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path, help='new build directory')
    parser.add_argument('--jobs', type=int, default=8)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    run_build(args.directory.resolve(), args.jobs, os.environ)
