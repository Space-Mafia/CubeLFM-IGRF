"""Build functions used by the dated fine-grid event."""
from pathlib import Path

import hashlib

import json

import os

import shutil

import subprocess

import sys

import time

import traceback

O = Path(__file__).resolve().parent

EVENT_ITMAX = 500

def select_mi_source(source, build_dir):
    """Apply the fine event policy to build copies of the archived sources."""
    if source.name == 'lfm_cube_event_solve_cap.F':
        text = source.read_text().rstrip() + """

      SUBROUTINE CUBE_EVENT_MI_POLICY(LIMIT,STATUS)
#include "ion_bridge_common.inc"
      INTEGER LIMIT,STATUS
C     Use the original setup and tolerance, then install the event limit.
      CALL LFM_CUBE_YY_MI100_POLICY(STATUS)
      IF(STATUS.EQ.0) ION_ITMAX_BR=LIMIT
      END
"""
    elif source.name == 'lfm_cube_ionosphere.f90':
        text = source.read_text()
        edits = (
            (" include 'native_mi_shape.inc'", " include 'native_mi_shape.inc'\n"
             f" integer,parameter::event_itmax={EVENT_ITMAX}", 1),
            ('call lfm_cube_yy_mi100_policy(local)', 'call cube_event_mi_policy(event_itmax,local)', 1),
            ('solve_ITMAX(2)=100', 'solve_ITMAX(2)=event_itmax', 2),
            ('solve_ITMAX(1)=100', 'solve_ITMAX(1)=event_itmax', 1),
            ('h%solve_ITMAX=100;', 'h%solve_ITMAX=event_itmax;', 1),
            ('6._real64,100._real64,1.e-3_real64', '6._real64,real(event_itmax,real64),1.e-3_real64', 2),
        )
        for old, new, count in edits:
            if text.count(old) != count:
                raise ValueError(f'Unexpected event MI source at {old!r}')
            text = text.replace(old, new)
    else:
        return source
    selected = build_dir/source.name
    selected.write_text(text)
    return selected

def read(path):
    return json.loads(Path(path).read_bytes())

def put(name, record):
    with (O / name).open('x') as stream:
        json.dump(record, stream, indent=2, allow_nan=False)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()

def run(label, args, cwd=O, env=None):
    args = list(map(str,args))
    put(label+'.START.json',dict(argv=args,cwd=str(cwd)))
    with (O/(label+'.log')).open('xb') as log:
        start=time.monotonic()
        proc=subprocess.run(args,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
        wall=time.monotonic()-start
    put(label+'.receipt.json',dict(returncode=proc.returncode,wall_s=wall,argv=args))
    if proc.returncode:
        raise RuntimeError(label+' failed, see '+str(O/(label+'.log')))
    return wall

def symbols(path):
    out=subprocess.check_output(['nm','-S','--defined-only',str(path)],text=True)
    return {v[3]:(int(v[0],16),int(v[1],16)) for line in out.splitlines() if len(v:=line.split())==4}

def build(C,source_dir,*,global_solver=False,branch_wall=False):
    if branch_wall and not global_solver:
        raise ValueError("Branch-local return requires the global ionospheric solver")
    io=Path(__file__).resolve().parent
    if not (io/'select_output.py').exists():
        io=io.parents[3]/'cubesphere/io'
    sys.path.insert(0,str(io))
    from select_output import prepare_driver
    B=O/'build'; B.mkdir()
    modules=B/'modules'; modules.mkdir()
    fc='/public/software/compiler/gcc/12.3.0/bin/gfortran'
    includes=['-I'+str(modules),'-I'+str(source_dir),*['-I'+p for p in C['module_dirs']]]
    common=['-O3','-cpp','-ffree-line-length-none','-fcheck=all','-fbacktrace','-ffpe-trap=invalid,zero,overflow','-fcheck=all,no-bounds','-pthread','-J'+str(modules),'-DLFM_CUBE_NATIVE_MI','-DLFM_CUBE_BORIS_C_KM_S=1100','-mavx2','-mtune=znver2','-ffp-contract=off','-fno-fast-math']
    new={}
    reused=dict(C['objects'])
    extra_link=[]
    module,main=prepare_driver(source_dir/'cube_lfm.f90',B)
    sources=[source_dir/name for name in C['new_sources'] if name!='cube_lfm.f90']+[module,main]
    if global_solver:
        enhanced=('lfm_cube_yy_fac_map','lfm_cube_yy_fac_source','cube_event_aurora_cap')
        reused={p:digest for p,digest in reused.items() if Path(p).stem not in enhanced}
        sources=[source_dir/(name+'.f90') for name in (*enhanced,'cube_event_global')]+sources
        mkl=Path('/public/software/compiler/intel/onepai/2021.3.0/mkl/2021.3.0')
        # YYG's public MI header also imports its original grid declarations.
        reference64 = read(Path(C['reference_64']['root']) / 'CONFIG.json')
        tensor_includes = [source_dir, Path(reference64['donor_include']), mkl/'include']
        cc=Path(fc).with_name('gcc')
        c_sources=[source_dir/(n+'.c') for n in ('checked_tensor','global_tensor')]
        if branch_wall:
            from global_mi.select_wall_return import select
            adapter=Path(__file__).resolve().parent/'global_mi'
            donor=adapter/'cube_yy_released_volume_cache.c'
            if not donor.is_file():
                donor=adapter.parents[2]/'fixes/released_volume_return_v1/cube_yy_released_volume_cache.c'
            original, = [p for p in reused if Path(p).name=='volume_cache.o']
            del reused[original]
            c_sources.extend([adapter/'branch_wall.c',select(donor,B/'volume_cache.c')])
        for src in c_sources:
            name=src.stem
            obj=B/(name+'.o')
            extra_cflags=(['-fopenmp',*['-I'+p for p in C['module_dirs']]] if name=='volume_cache' else [])
            run('COMPILE_'+name,[cc,'-O3','-std=gnu11','-fPIC','-fno-fast-math','-DGAMERA_MI_SPARSE_DIRECT',
                                 *['-I'+str(p) for p in tensor_includes],*extra_cflags,'-c',src,'-o',obj],B)
            new[str(obj)]=sha(obj)
            put('OBJECT_'+name+'.json',dict(path=str(obj),source_sha256=sha(src),sha256=new[str(obj)]))
        extra_link=['-L'+str(mkl/'lib/intel64'),'-Wl,-rpath,'+str(mkl/'lib/intel64'),'-lmkl_rt','-ldl']
    from fixes.induction_reference.materialize import select_event_evolution
    for source in sources:
        name=source.name
        if not global_solver:
            source=select_mi_source(source,B)
        if name=='lfm_cube_evolution.f90':
            source=select_event_evolution(source,B)
            if global_solver and 'call induction_reference_edge(' not in source.read_text():
                raise ValueError('Selected event evolution source lacks the corrected induction reference')
        obj=B/(source.stem+'.o')
        flags=C['native128_fixed_flags'] if name.endswith('.F') else [*includes,*common]
        run('COMPILE_'+source.stem,[fc,*flags,'-c',source,'-o',obj],B)
        new[str(obj)]=sha(obj)
        put('OBJECT_'+source.stem+'.json',dict(path=str(obj),source_sha256=sha(source),sha256=new[str(obj)]))
    exe=B/'cube_event128'
    wraps=['cube_yy_direct_reference_evaluate','cube_yy_magnetic_wall_create','cube_yy_magnetic_wall_apply','cube_yy_magnetic_wall_destroy','cube_yy_igrf_configure','cube_yy_igrf_field']
    if global_solver:
        wraps.append('cube_yy_footpoint_potential')
    hlib=str(Path(sys.executable).resolve().parents[1]/'lib')
    run('LINK_EVENT128',[fc,*C['link_runtime_flags'],*new,*reused,'-Wl,--start-group',*C['libraries'],'-Wl,--end-group',*C['MPI_link_inputs'],*['-Wl,--wrap='+n for n in wraps],'-L'+hlib,'-Wl,-rpath,'+hlib,'-lhdf5',*extra_link,'-fopenmp','-lm','-Wl,-Map,'+str(B/'model.map'),'-o',exe],B)
    actual=symbols(exe)
    fine=symbols(C['reference_128']['exe'])
    for name in ['scrach_','scrch1_','lfm_rescue_mdims_','lfm_rescue_help_i_','lfm_rescue_help_r_','lfm_rescue_gdims_i_','lfm_rescue_gdims_l_']:
        assert actual[name][1]==fine[name][1] and actual[name][0]%32==0,name
    for name,size in {'lfm_rescue_ion_pn_':97*67*8,'lfm_rescue_ion_x_':4*97*129*8,'lfm_rescue_ion_rho_':3*96*128*8}.items():
        assert actual[name][1]==size,name
    for name in ['cube_event_psislv_','lfm_cube_event_solve_cap_','cube_event_dpb_step','cube_yy_aurora_point','cube_event_wind_sample']:
        assert name in actual,name
    mapping=(B/'model.map').read_text()
    assert 'liblfm_cube_contracts.a(lfm_cube_mpi.f90.o)' in mapping
    assert 'liblfm_cube_evolution.a(lfm_cube_evolution.f90.o)' not in mapping
    result=dict(exe=str(exe),sha256=sha(exe),new_objects=new,reused_objects=reused,libraries=C['libraries'])
    if branch_wall:
        if 'cube_project_branch_voltage' not in actual:
            raise ValueError('Missing branch-local return projection')
        result['wall_return']='branch-local-electrostatic'
    if global_solver:
        for name in ('gamera_mi_global_solve_tensor_checked','__wrap_cube_yy_footpoint_potential'):
            if name not in actual: raise ValueError('Missing global MI entry '+name)
        result.update(potential_solver='YYG global PARDISO',global_lon_lat=[192,97],
                      relative_tolerance=1.e-10,absolute_tolerance=1.e-8)
    else:
        result.update(ITMAX=EVENT_ITMAX,ERRMAX=C['ERRMAX'])
    put('BUILD_PASS.json',result)
    return result
