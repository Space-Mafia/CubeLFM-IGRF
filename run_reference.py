"""Fresh reference execution using the accepted Sugon ELF; no compilation."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,time
root=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
r=json.loads((root/'REFERENCE_RUNTIME.json').read_text());case=json.loads((root/'CASE_CONTRACT.json').read_text())
assert os.environ.get('SLURM_JOB_ID'),'Run inside a compute-node allocation'
assert os.environ.get('OMP_NUM_THREADS')=='1'and int(os.environ['SLURM_NTASKS'])==96
assert sha(r['binary'])==r['binary_sha256'],'Accepted executable differs; qualify a new build independently'
for name,h in r['inputs'].items():assert sha(Path(r['input_root'])/name)==h,name
parent=args.output.resolve().parent;v=os.statvfs(parent);assert v.f_bavail*v.f_frsize>=r['forecast_output_bytes'],'Measured available space does not fit own forecast outputs'
out=args.output.resolve();out.mkdir(exist_ok=False);(out/'rank_records').mkdir()
for name in r['inputs']:shutil.copy2(Path(r['input_root'])/name,out/name)
shutil.copy2(root/'rank_exec.py',out/'rank_exec.py');(out/'CASE_CONTRACT.json').write_text(json.dumps(case,indent=2)+'\n')
(out/'REPRODUCTION_ORIGIN.json').write_text(json.dumps(dict(source_tag='cubelfm-igrf-v0.2.0',new_job=os.environ['SLURM_JOB_ID'],reference_data_job='29433880',old_models_not_executed=True,binary_sha256=r['binary_sha256']),indent=2)+'\n')
cmd=['mpirun','--mca','plm','slurm','--mca','pml','ob1','--mca','btl','self,vader','--nooversubscribe','--map-by','ppr:12:numa','--bind-to','core','-np','96','python3',str(out/'rank_exec.py'),r['binary'],'shell_geometry.txt','frame','24','3','43200','1000000','120','-',r['case_fingerprint'],'3600','ionosphere_reference.txt','-','4','2','2','imf_driver.txt']
env=dict(os.environ);env['CUBE_SWEEP_CACHE']='1';tick=time.monotonic()
with(out/'MHD.log').open('x')as log:completed=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT)
wait=time.monotonic()-tick;(out/'IMMEDIATE_PROCESS_WAIT.json').write_text(json.dumps(dict(returncode=completed.returncode,full_wait_s=wait,init_output_exit_included=True,scientific_audit_pending=True),indent=2)+'\n')
raise SystemExit(completed.returncode)
