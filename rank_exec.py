from pathlib import Path
import json,os,socket,sys,time
q=Path(__file__).resolve().parent
rank=int(os.environ['OMPI_COMM_WORLD_RANK']);job=os.environ['SLURM_JOB_ID'];cores=set()
for cpu in os.sched_getaffinity(0):
 p=Path('/sys/devices/system/cpu')/('cpu'+str(cpu))/'topology';cores.add(((p/'physical_package_id').read_text().strip(),(p/'core_id').read_text().strip()))
assert len(cores)==1 and os.environ['OMP_NUM_THREADS']=='1'
numa=set()
for cpu in os.sched_getaffinity(0):
 for n in (Path('/sys/devices/system/cpu')/('cpu'+str(cpu))).glob('node[0-9]*'):numa.add(n.name)
assert len(numa)==1
x=dict(rank=rank,pid=os.getpid(),job_id=job,host=socket.gethostname().split('.')[0],cores=sorted(cores),cpus=sorted(os.sched_getaffinity(0)),numa=sorted(numa))
record=q/'rank_records'/('rank-%03d.json'%rank)
with record.with_suffix('.tmp').open('x')as f:json.dump(x,f)
os.rename(record.with_suffix('.tmp'),record)
ready=q/'rank_records'/'ALL_96_DISTINCT_PHYSICAL_CORES.json';t=time.monotonic()
if rank==0:
 while len(list((q/'rank_records').glob('rank-*.json')))<96:
  assert time.monotonic()-t<120;time.sleep(.05)
 rows=[json.loads(p.read_text())for p in sorted((q/'rank_records').glob('rank-*.json'))]
 assert {r['rank']for r in rows}==set(range(96));assert {r['job_id']for r in rows}=={job}
 assert len({(r['host'],*r['cores'][0])for r in rows})==96
 counts={}
 for r in rows:
  node=r['numa'][0];counts[node]=counts.get(node,0)+1
 assert len(counts)==8 and set(counts.values())=={12},counts
 with ready.open('x')as f:json.dump(dict(status='96_DISTINCT_PHYSICAL_CORES_8_NUMA_12_EACH_BEFORE_MHD',ranks=96,NUMA_counts=counts),f)
else:
 while not ready.exists():
  assert time.monotonic()-t<120;time.sleep(.05)
os.execv(sys.argv[1],sys.argv[1:])
