"""Inspect the final compact frame and checkpoint of a completed 128-grid event."""
import argparse
import csv
import json
from pathlib import Path
import h5py
import numpy as np


def check(directory):
    process = json.loads((directory/'process.json').read_text())
    if process['exit'] != 0:
        raise ValueError('The MPI process failed')
    with (directory/'frame-steps.csv').open() as stream:
        step, time_s, dt, *_ = list(csv.reader(stream))[-1]
    step, time_s, dt = int(step), float(time_s), float(dt)
    if step <= 0 or time_s != process['seconds'] or not np.isfinite(dt) or dt <= 0:
        raise ValueError('The model did not advance to the requested endpoint')
    with (directory/'frame-mi-updates.csv').open() as stream:
        mi = list(csv.reader(stream))
    if len(mi)<3:
        raise ValueError('Expected at least one natural MI update after initialization')
    with h5py.File(directory/f'frame-{step:07d}.h5') as frame:
        if frame.attrs['time_s'] != time_s or frame.attrs['ranks'] != 192:
            raise ValueError('Incorrect frame time or MPI size')
        for rank in range(192):
            tile = frame[str(rank)]
            state = tile['primitive'][:]
            if state.shape!=(8,24,24,16) or not np.isfinite(state).all() or state[0].min()<=0 or state[4].min()<=0:
                raise ValueError(f'Invalid MHD fields on rank {rank}')
            if not np.isfinite(tile['ct_flux'][:]).all():
                raise ValueError(f'Invalid CT flux on rank {rank}')
        global_mi=frame['global_mi']
        for name in ('potential_V','outward_fac_A_m2','pedersen_S','hall_S'):
            field=global_mi[name][:]
            if field.shape!=(97,192) or not np.isfinite(field).all():
                raise ValueError(f'Invalid global MI field {name}')
        residual=float(global_mi.attrs['final_residual'])
        held=float(global_mi.attrs['held_s'])
        if not np.isfinite(residual) or held!=float(mi[-1][1]):
            raise ValueError('Invalid global residual or inconsistent held clock')
    with h5py.File(directory/f'frame-step{step:07d}.h5') as cp:
        if cp.attrs['format']!='CubeLFM-checkpoint-1' or 'manifest' not in cp:
            raise ValueError('Missing compact checkpoint manifest')
    return dict(physical_s=time_s,steps=step,final_dt_s=dt,mpi_ranks=192,mi_rows=len(mi)-1,
                held_s=held,global_residual=residual,full_process_wall_s=process['full_process_wall_s'])


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    result=check(parser.parse_args().directory)
    print(json.dumps(result,indent=2,allow_nan=False))
