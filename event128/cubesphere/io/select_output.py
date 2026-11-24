"""Wire the CubeLFM output interface into a selected (possibly archived) driver."""
from pathlib import Path


def with_compact_output(text):
    if 'use cube_compact_output' in text:
        return text
    replacements = [
        (' use mpi_f08', ' use cube_compact_output,only:cube_output\n use mpi_f08', 1),
        (' type(cube_model)::c',' type(cube_output)::archive\n type(cube_model)::c', 1),
        (' if(len_trim(restart)>0.and.trim(restart)/="-") call c%restore(trim(restart),fingerprint)',
         ' call archive%open(c,trim(path))\n if(len_trim(restart)>0.and.trim(restart)/="-") call archive%restore(c,trim(restart),fingerprint)', 1),
        ('call c%dump(trim(path));call c%dump_mi(trim(path))','call archive%frame(c)',2),
        ('call c%checkpoint(trim(path),fingerprint)','call archive%checkpoint(c,fingerprint)',1),
        (' call MPI_Finalize(st)',' call archive%close(c)\n call MPI_Finalize(st)',1),
    ]
    for old,new,count in replacements:
        if text.count(old)!=count:
            raise ValueError('Unsupported Cube driver at '+old)
        text=text.replace(old,new)
    return text


def prepare_driver(source,build):
    io=Path(__file__).resolve().parent
    main=build/'cube_lfm.f90'
    main.write_text(with_compact_output(Path(source).read_text()))
    helper=build/'compact_io.py'
    helper.write_bytes((io/'compact_io.py').read_bytes())
    return io/'cube_compact_output.f90',main
