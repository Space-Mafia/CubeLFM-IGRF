"""Allocated-only first-version audit; frozen readers plus explicit MI100/v4 policy."""
from pathlib import Path
import json,sys,hashlib,os,numpy as np
from frozen_readers_igrf_v5 import read_frame,read_mi_snapshot,read_compact_epoch,csv,solve_indices,rank_tile
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def audit(case,contract_path):
 case=Path(case);C=json.loads(Path(contract_path).read_text())
 c={k:v for k,v in C.items()if k!='fingerprint'}
 assert hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()==C['fingerprint']
 assert C['ITMAX']==100 and C['ERRMAX']==.001 and C['finite_nonconvergence']=='warning_continue_honest_false'
 assert C['SigmaP']==5 and C['Hall']==0 and C['MI_cadence_s']==10 and C['MI_clip']is False
 assert C['rotation_period_s']==C['physical_seconds']==43200
 assert C['omega_rad_s']==2*np.pi/43200 and C['phase0_rad']==0
 assert C['wind']==dict(Bz_before_nT=5,Bz_after_nT=-5,switch_s=7200,time_reference='EARTH_x0_convected_native_real32_NINT')
 for n,h in C['inputs'].items():assert sha(case/n)==h,n
 clock=csv(case/'frame-steps.csv',['step','time_s','dt_s','cfl_dt_s','rho_min','c_min','CT_step','shared_face','divB_relative'])
 assert clock['step'][0]==clock['time_s'][0]==0 and clock['time_s'][-1]==C['physical_seconds']
 assert np.array_equal(clock['step'],np.arange(len(clock)))
 assert np.all(clock['dt_s'][1:]>0) and np.array_equal(clock['time_s'][:-1]+clock['dt_s'][1:],clock['time_s'][1:])
 assert np.array_equal(clock['dt_s'][1:-1],np.minimum(6.,clock['cfl_dt_s'][1:-1]))
 assert np.min(clock['rho_min'])>0 and np.min(clock['c_min'])>0
 for k in ['CT_step','shared_face','divB_relative']:assert np.max(clock[k])<1e-11,k
 mi=csv(case/'frame-mi-updates.csv',['step','time_s','SER_N','SER_S','ITER_N','ITER_S','ITMAX_N','ITMAX_S','converged_N','converged_S'])
 indices=solve_indices(clock,10.)
 assert len(indices)==4321
 assert np.array_equal(mi['step'],clock['step'][indices]) and np.array_equal(mi['time_s'],clock['time_s'][indices])
 assert np.array_equal(mi['updates'],np.arange(1,len(mi)+1))
 for h in ['N','S']:
  assert np.all(mi['ITMAX_'+h]==100) and np.all(mi['SER_'+h]>=0)
  assert np.array_equal(mi['converged_'+h],mi['SER_'+h]<=.001)
  assert np.all(mi['ITER_'+h]>=0) # Original nested ITER is deliberately not clamped to ITMAX.
 assert np.array_equal(mi['relative_residual'],np.maximum(mi['SER_N'],mi['SER_S']))
 layout=dict(mpi_count=96,dims=[4,2,2],global_shape=[64,24,24],explicit=True)
 owners={};frame_hashes={};maxdiv=0.;zero_frames=[]
 for p in sorted(case.glob('frame-*-r*.bin')):
  x=read_frame(p);rank=int(p.stem.rsplit('-r',1)[1]);face,first,tile=rank_tile(rank,layout)
  assert x['face']==face and x['first'].tolist()==first and x['n'].tolist()==tile and x['div']<1e-11
  if x['zero_flux']:zero_frames.append(dict(file=p.name,time_s=x['time'],rank=rank))
  maxdiv=max(maxdiv,x['div']);owners.setdefault(x['time'],[]).append(rank);frame_hashes[p.name]=sha(p)
 archive_indices=solve_indices(clock,float(C['output_cadence_s']))
 assert len(archive_indices)==361
 assert set(owners)==set(clock['time_s'][archive_indices]) and all(sorted(v)==list(range(96))for v in owners.values())
 snapshots=[]
 for p in sorted(case.glob('frame-mi-step*.bin')):
  x=read_mi_snapshot(p);matches=mi[mi['step']==x['step']];assert len(matches)==1
  row=matches[0];assert x['time']==x['last_time']==float(row['time_s'])
  assert np.allclose(x['CPCP_kV'],[row['CPCP_N_kV'],row['CPCP_S_kV']],rtol=0,atol=1e-12)
  snapshots.append(dict(file=p.name,time_s=x['time'],sha256=sha(p)))
 assert sorted(x['time_s']for x in snapshots)==sorted(owners)
 # This is a reader layout view, not a fabricated alternate physical manifest.
 manifest=dict(solver='cube_LFM',fingerprint=C['fingerprint'],grid=[64,24,24,6],execution_layout=layout,
  rotation=dict(omega_rad_s=C['omega_rad_s'],phase0_rad=C['phase0_rad']),
  physical=dict(mpi_ranks=96,re_cm=6.38e8,dipole_equatorial_g=.30,SW=dict(rho_g_cm3=8.3631096e-24,v_km_s=[-400.,0.,0.],c_km_s=40.,B_nT=[0.,0.,5.])))
 markers=sorted(case.glob('frame-step*.complete'));assert len(markers)==12
 cp_indices=solve_indices(clock,float(C['checkpoint_cadence_s']))[1:]
 assert len(cp_indices)==12
 assert sorted(float(m.read_text().splitlines()[1].split()[1])for m in markers)==sorted(clock['time_s'][cp_indices])
 checks=[]
 for marker in markers:
  states,records=read_compact_epoch(marker,manifest)
  assert len(states)==96 and {x['rank']for x in states}==set(range(96))
  assert len({x['mi']['state_bytes']for x in states})==1
  x=states[0];i=x['step'];assert x['time_s']==clock['time_s'][i] and x['dt']==clock['dt_s'][i]
  held=x['mi'];rows=mi[mi['step']==held['last_step']];assert len(rows)==1;row=rows[0]
  assert held['last_time']==row['time_s'] and held['residual']==row['relative_residual']
  assert np.array_equal(held['native_SER'],[row['SER_N'],row['SER_S']])
  assert np.array_equal(held['converged'],[row['converged_N'],row['converged_S']])
  checks.append(dict(marker=marker.name,step=i,time_s=x['time_s'],MI_last_time=held['last_time'],shards=records))
 result=dict(status='CUBE64_IGRF_SPIN12H_LONG12H_ARRAY_CLOCK_MI100_CHECKPOINT_PASS',steps=int(clock['step'][-1]),end_time_s=float(C['physical_seconds']),ranks=96,
  archived_frames=len(frame_hashes),frame_clocks=len(owners),natural_MI_solves=len(mi),checkpoint_clocks=len(markers),max_divB=maxdiv,
  max_CT_step=float(np.max(clock['CT_step'])),max_shared_face=float(np.max(clock['shared_face'])),
  MI_policy='LFM100/1e-3 finitewarning_continue_honest_converged',max_SER_N=float(np.max(mi['SER_N'])),max_SER_S=float(np.max(mi['SER_S'])),
  warned_solves_N=int(np.sum(mi['converged_N']==0)),warned_solves_S=int(np.sum(mi['converged_S']==0)),
  frame_sha256=frame_hashes,MI_snapshots=snapshots,checkpoints=checks,zero_flux_frames=zero_frames,full_YinYang_benchmark_complete=False)
 output=Path(os.environ.get('YY_AUDIT_RESULT_PATH',str(case/'LONG12H_NUMERICAL_AUDIT.json')))
 with output.open('x')as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
 return result

if __name__=='__main__':
 assert sys.argv[1]=='audit'
 r=audit(sys.argv[2],sys.argv[3])
 print(json.dumps({k:v for k,v in r.items()if k not in ['frame_sha256','MI_snapshots','checkpoints']}),flush=True)
