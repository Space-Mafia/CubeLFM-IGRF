"""Prepare the dated solar wind for the 128-grid reference event."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--raw-csv', type=Path, help='NASA HAPI CSV for an offline compute node')
    args = parser.parse_args()
    data = args.directory.resolve()
    data.mkdir(parents=True, exist_ok=True)
    case = json.loads((HERE/'CASE.json').read_text())
    model = data/'event_model.json'
    model.write_text(json.dumps(case['event_model'], indent=2)+'\n')
    command = [sys.executable, str(HERE/'driver/fetch_nasa_omni_event.py'),
               '--simulation-start', case['start_utc'], '--simulation-stop', case['stop_utc'],
               '--padding-minutes', '160', '--event-model', str(model),
               '--domain-radius-re', str(3515454736.38853/6371000.0),
               '--output', str(data/'omni.h5'), '--no-qa-plot']
    if args.raw_csv:
        command += ['--raw-csv', str(args.raw_csv.resolve(strict=True))]
    subprocess.run(command, check=True)


if __name__ == '__main__':
    main()
