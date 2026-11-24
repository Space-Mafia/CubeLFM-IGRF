"""Build the global-MI event with the reference case's MHD resolution."""
import json
import sys
from pathlib import Path

from materialize import materialize, once

here = Path(__file__).resolve().parent
event = here.parent
repo = here.parents[4]
sys.path[:0] = [str(event), str(repo/'cubesphere/io')]


def main(reference, source, equipotential=False, branch_wall=False):
    case = Path.cwd()
    selected = case/'source'
    materialize(source, selected)
    cfg=json.loads((reference/'CONFIG.json').read_text())
    if cfg.get('grid') != [128, 48, 48, 6] or equipotential or not branch_wall:
        raise ValueError('This release builds the 128-grid global event with the branch-local wall')
    import run_event128
    run_event128.O = case
    run_event128.build(cfg, selected, global_solver=True, branch_wall=True)

    install_global_writer(case / 'build')


def install_global_writer(directory):
    writer = directory / 'compact_io.py'
    text = once(writer.read_text(), '            write_mi(out,mi,step,t)',
                "            write_mi(out,mi,step,t)\n            from archive_global import write_global\n            write_global(out,stage/'global-state')")
    writer.write_text(once(text, '        mi.unlink()', "        mi.unlink()\n        (stage/'global-state').unlink()"))
    archive_src=(here/'archive_global.py').read_text()
    writer.with_name('archive_global.py').write_text(archive_src)


if __name__ == '__main__':
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), '--equipotential' in sys.argv[3:], '--branch-local-wall' in sys.argv[3:])
