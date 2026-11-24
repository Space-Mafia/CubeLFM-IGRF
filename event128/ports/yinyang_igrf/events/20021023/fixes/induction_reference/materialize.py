"""Select the induction reference independently of the physical cold-start flux."""
from pathlib import Path
import sys
import hashlib


def apply_to_text(text):
    start = text.index(' subroutine prepare_hydro_time_reference(c)')
    stop = text.index(' subroutine prepare_zero_reference_records(c)', start)
    routine = text[start:stop]
    begin = routine.index('   do d=1,3\n    do k=0,n(3)')
    end = routine.index('   deallocate(reference_faces)', begin) + len('   deallocate(reference_faces)')
    reconstruction = routine[begin:end]
    s = reconstruction.index('       line=(b(1)-a(1))*gint')
    e = reconstruction.index('\n      enddo',s)
    reconstruction = reconstruction[:s] + """       ! A straight segment with both endpoints inside the ball stays inside.
       ! Its existing field and line integral already use unit taper.
       if(max(norm2(a),norm2(b))<=16._real64*c%rp)cycle
       call induction_reference_edge(a,b,t,ap(i,j,k,d),edge_vector(:,i,j,k,d))""" + reconstruction[e:]
    s = reconstruction.index('   if(c%mode==3.and.sample==1)')
    e = reconstruction.index('   delta_phi=phi_ref-c%background_flux', s)
    reconstruction = reconstruction[:s] + """   call c%faces%sync(phi_ref,st);call demand(st==0,'uncut induction reference face ownership')
""" + reconstruction[e:]
    reconstruction = reconstruction.replace('      yy_hydro_cell_coef(:,sample,i,j,k)=reference_channels(i,j,k,6:8)\n','')
    reconstruction = reconstruction.replace('      ! Original prescribed planetary taper is identically zero here;\n      ! nearest active planetary reference also gives the inflow value0.', '      ! Extend the numerical reference with the original radial ghost convention.')
    # The physical, tapered solve has already populated all HYDRO, FAC and initial-flux arrays.
    # Reuse its temporary edge samples and the same owner/CT/ghost reconstruction procedure.
    reconstruction = reconstruction.replace('   call cube_curl_update(', '   phi_ref=0\n   call cube_curl_update(')
    assert 'yy_hydro_' not in reconstruction and 'yy_planetary_initial_phi' not in reconstruction
    key = "   if(c%rank==0)write(*,*)'INIT_CT_MATCHED_HYDRO_REFERENCE_PHASE',sample,ns"
    assert routine.count(key)==1
    routine = routine.replace(key, '   if(c%mode>=2)then\n' + reconstruction + '\n   endif\n' + key)
    edge = Path(__file__).with_name('bare_edge.f90.inc').read_text()
    return text[:start]+edge+routine+text[stop:]


def select_event_evolution(source, build):
    """Apply the accepted induction fix to the two known pre-fix event sources."""
    source = Path(source)
    data = source.read_bytes()
    original = {
        'a32608e2457a11ce390788a17c0160b97c6469b8d594f8624c241b18c57e26d2',
        '8f79eee5c9e97310d3d90568fbe4a95d302bd4e13f86a30d05ee608a26e6911e',
    }
    if hashlib.sha256(data).hexdigest() not in original:
        return source
    selected = build / source.name
    selected.write_text(apply_to_text(data.decode()))
    return selected


if __name__=='__main__':
    Path(sys.argv[2]).write_text(apply_to_text(Path(sys.argv[1]).read_text()))
