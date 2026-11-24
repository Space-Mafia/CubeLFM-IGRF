"""Compile the fine corrected global event, then launch with the completed event's inputs."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import deque
from pathlib import Path

event=Path(__file__).resolve().parents[1]

def run(label,argv,cwd,root):
    argv=list(map(str,argv))
    with (root/(label+'.log')).open('xb') as output:
        start=time.monotonic()
        p=subprocess.run(argv,cwd=cwd,stdout=output,stderr=subprocess.STDOUT)
        wall=time.monotonic()-start
    (root/(label+'.receipt.json')).write_text(json.dumps(dict(argv=argv,returncode=p.returncode,wall_s=wall)))
    p.check_returncode()
    return wall

def main(root,reference,source=None,branch_wall=False):
    original=json.loads((reference/'event12h.receipt.json').read_text())
    if original['returncode'] or not (reference/'PHYSICAL_FULL12H_CLOSED.json').is_file():
        raise ValueError('The baseline MHD process has not closed successfully')
    argv=original['argv'];inputs=argv[-16:];mpi=argv[:argv.index('-np')+2]
    for i in (0,10,15):inputs[i]=str((reference/'event12h'/inputs[i]).resolve(strict=True))
    build_args=[sys.executable,event/'global_mi/build.py',reference,source or reference]
    if branch_wall:
        build_args.append('--branch-local-wall')
    run('build',build_args,root,root)
    build=json.loads((root/'BUILD_PASS.json').read_text());exe=build['exe']
    inputs[8]=hashlib.sha256((inputs[8]+build['sha256']).encode()).hexdigest()
    os.environ.update(CUBE_SWEEP_CACHE='1',CUBELFM_OUTPUT_FORMAT='compact',CUBELFM_IO_PYTHON=sys.executable,
                      CUBELFM_IO_TOOL=str(root/'build/compact_io.py'),MKL_NUM_THREADS='1',MKL_THREADING_LAYER='GNU')
    case=root/'event12h';case.mkdir()
    wall=run('event12h',mpi+[exe]+inputs,case,root)
    with (case/'frame-steps.csv').open('rb') as f:
        f.seek(max(0,os.fstat(f.fileno()).st_size-10000))
        last=deque(f,maxlen=1)[0].decode().split(',')
    if float(last[1])!=43200:raise ValueError('Fine global run ended before 16 UT')
    (root/'PHYSICAL_CLOSED.json').write_text(json.dumps(dict(seconds=43200,full_wait_s=wall,rt=43200/wall,
        data_job=os.environ['SLURM_JOB_ID'],step=int(last[0]),reference=str(reference))))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference',type=Path)
    parser.add_argument('--source',type=Path)
    parser.add_argument('--branch-local-wall',action='store_true')
    args=parser.parse_args()
    main(Path.cwd(),args.reference,args.source,args.branch_local_wall)
