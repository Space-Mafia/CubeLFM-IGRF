"""Reproduce the fine event's compiler inputs without a compiler, simulation or network."""
import ast
import hashlib
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parent
EVENT = ROOT / 'ports/yinyang_igrf/events/20021023'


def read(name):
    return json.loads((ROOT / name).read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(path, expected):
    if digest(path) != expected:
        raise ValueError(f'SHA-256 mismatch: {path}')


def main():
    sys.dont_write_bytecode = True
    sys.path[:0] = [str(EVENT / 'global_mi'), str(EVENT), str(ROOT / 'cubesphere/io')]
    from materialize import materialize
    from build import install_global_writer
    from fixes.induction_reference.materialize import select_event_evolution
    from global_mi.select_wall_return import select
    from select_output import prepare_driver

    for name, value in read('FILES.json').items():
        check(ROOT / name, value)
    build = read('BUILD.json')
    for record in build['inherited'].values():
        check(ROOT.parent / record['repo_source'], record['source_sha256'])
    compiler = (EVENT / 'run_event128.py').read_text()
    routine = next(n for n in ast.parse(compiler).body if getattr(n, 'name', '') == 'build')
    actual = hashlib.sha256(ast.get_source_segment(compiler, routine).encode()).hexdigest()
    if actual != build['build_function_sha256']:
        raise ValueError('The 128-grid build routine differs from the recorded version')

    with TemporaryDirectory() as tmp:
        selected = Path(tmp) / 'source'
        output = Path(tmp) / 'build'
        output.mkdir()
        materialize(ROOT / 'reference-source', selected)
        module, driver = prepare_driver(selected / 'cube_lfm.f90', output)
        evolution = select_event_evolution(selected / 'lfm_cube_evolution.f90', output)
        donor = ROOT / 'ports/yinyang_igrf/fixes/released_volume_return_v1/cube_yy_released_volume_cache.c'
        cache = select(donor, output / 'volume_cache.c')
        sources = {p.name: p for p in (module, driver, evolution, cache, EVENT / 'global_mi/branch_wall.c')}
        for name, expected in build['compile_sources'].items():
            check(sources.get(name, selected / name), expected)
            check(ROOT / 'compiled' / name, expected)
        install_global_writer(output)
        for name, expected in build['output_writers'].items():
            check(output / name, expected)
            check(ROOT / 'runtime' / name, expected)
    print(f"Verified {len(build['compile_sources'])} regenerated compiler inputs, "
          f"{len(build['inherited'])} inherited source mappings, and both output writers.")


if __name__ == '__main__':
    main()
