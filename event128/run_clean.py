"""Cold-start the dated 128-grid case with the freshly built executable."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('inputs', type=Path)
    parser.add_argument('directory', type=Path, help='new output directory')
    parser.add_argument('--seconds', type=int, default=43200, help='seconds after 04:00 UT')
    args = parser.parse_args()
    if not 0 < args.seconds <= 43200:
        parser.error('--seconds must be in the interval 1 to 43200')
    exe = (args.build/'cube_event128').resolve(strict=True)
    wind = (args.inputs/'omni.h5').resolve(strict=True)
    mesh = HERE/'inputs/geometry/shell_geometry.txt'
    ion = HERE/'inputs/geometry/ionosphere_reference.txt'
    fingerprint = hashlib.sha256()
    for path in (exe, wind, mesh, ion, HERE/'CASE.json'):
        fingerprint.update(hashlib.sha256(path.read_bytes()).digest())
    args.directory.mkdir(parents=True)
    mpi = shlex.split(os.environ.get('CUBELFM_MPI_LAUNCH', 'mpirun --bind-to core'))
    cmd = [*mpi,'-np','192',str(exe),str(mesh),'frame','48','3',str(args.seconds),
           '10000000','120','-',fingerprint.hexdigest(),'3600',str(ion),'-','8','2','2',str(wind)]
    env = dict(os.environ, CUBE_SWEEP_CACHE='1', CUBELFM_OUTPUT_FORMAT='compact',
               CUBELFM_IO_PYTHON=sys.executable, CUBELFM_IO_TOOL=str(HERE/'runtime/compact_io.py'),
               OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', MKL_THREADING_LAYER='GNU')
    start = time.monotonic()
    result = subprocess.run(cmd, cwd=args.directory, env=env)
    elapsed = time.monotonic()-start
    (args.directory/'process.json').write_text(json.dumps(
        dict(exit=result.returncode,seconds=args.seconds,full_process_wall_s=elapsed,
             mpi_ranks=192,fingerprint=fingerprint.hexdigest()),indent=2)+'\n')
    result.check_returncode()


if __name__ == '__main__':
    main()
