"""Select the global MI adapter in an isolated copy of the dated event sources."""
from pathlib import Path
import shutil
here=Path(__file__).resolve().parent

def once(s,old,new):
    if s.count(old)!=1:
        raise ValueError(f"Unexpected source at {old[:90]!r}")
    return s.replace(old,new)

def materialize(source,target):
    target.mkdir()
    for f in source.iterdir():
        if f.is_file() and f.suffix.lower() in ('.f','.f90','.c','.h','.inc'):
            (target/f.name).symlink_to(f.resolve())
    for name in ('lfm_cube_yy_fac_map.f90','lfm_cube_yy_fac_source.f90','cube_event_aurora_cap.f90'):
        if not (target/name).exists():
            (target/name).symlink_to((here.parent/'native'/name).resolve(strict=True))
    def edit(name,fn):
        src=target/name
        new=fn(src.read_text())
        if src.is_symlink():src.unlink()
        src.write_text(new)
    for f in [*sorted((here/'yyg').iterdir()),here/'global_tensor.c',here/'cube_event_global.f90']:
        if f.is_file():
            dest=target/f.name
            if dest.is_symlink(): dest.unlink()
            shutil.copy2(f,dest)
    edit('cube_event_aurora_cap.f90',lambda s:once(s,' public::event_cap_update',' public::auroral_point,feet,advance_dpb,masks,precip\n public::event_cap_update'))
    edit('lfm_cube_yy_fac_map.f90',lambda s:once(once(s,' public::yy_released_source_map',' public::yy_map_global_points,yy_released_source_map'),'end module', (here/'global_map.inc').read_text()+'end module'))
    def fac(s):
        s=once(s,'only:yy_released_source_map','only:yy_released_source_map,yy_map_global_points')
        s=once(s,'public::yy_source_prepare','public::yy_source_global_moments,yy_source_prepare')
        return once(s,'end module', ''' subroutine yy_source_global_moments(t,theta,radii,rion,v,samples,status)
  real(real64),intent(in)::t,theta(:),radii(:),rion,v(:,:)
  real(real64),intent(out)::samples(:,:)
  integer,intent(out)::status
  status=1
  if(.not.source_ready.or.t/=source_epoch)return
  call yy_map_global_points(mesh,length_m,global_parallel,t,theta,radii,rion,v,samples,status)
 end subroutine
end module''')
    edit('lfm_cube_yy_fac_source.f90',fac)
    def ion(s):
        s=once(s,' use cube_event_aurora_cap',' use cube_event_global\n use cube_event_aurora_cap')
        s=once(s,'integer::p,first,last,total,ranks','integer::p,first,last,total,ranks,iterations')
        loc='  if(h%rank==0) then\n   rn=reshape'
        s=once(s,loc,'''  call global_update(h%comm,time,h%source_theta,h%source_radii,rion,h%event_ready,h%event_radius,h%event_edge,iterations,status)
  if(status/=0)return
'''+loc)
        start=s.index('   oldpot(:,:,1)='); end=s.index('    h%worst_residual=',start)
        s=s[:start]+'''   call global_caps(time,h%cap_x,h%cap_y,h%pn,h%ps,h%event_p,h%event_h,h%event_diag)
   if(local==0)then
    h%solve_SER=global_metrics(1);h%solve_ITER=iterations;h%solve_ITMAX=20000;h%solve_converged=1
'''+s[end:]
        s=once(s,'buffer=[6._real64,100._real64,1.e-3_real64,1._real64','buffer=[7._real64,20000._real64,1.e-10_real64,1.e-8_real64')
        s=once(s,'reshape(h%event_diag,[2*event_channels*capj*capk])]','reshape(h%event_diag,[2*event_channels*capj*capk]),reshape(global_fields,[global_size]),global_metrics]')
        s=once(s,'(10+2*event_channels)*capj*capk) return','(10+2*event_channels)*capj*capk+global_size+4) return')
        s=once(s,'[6._real64,100._real64,1.e-3_real64,1._real64]','[7._real64,20000._real64,1.e-10_real64,1.e-8_real64]')
        s=once(s,'any(buffer(13:14)/=real(merge(1,0,buffer(9:10)<=1.e-3_real64),real64))','any(buffer(13:14)/=1._real64)')
        s=once(s,'h%solve_ITMAX=100;n=82','h%solve_ITMAX=20000;n=82')
        s=once(s,'h%event_diag=reshape(buffer(n+1:),shape(h%event_diag))','''h%event_diag=reshape(buffer(n+1:n+size(h%event_diag)),shape(h%event_diag));n=n+size(h%event_diag)
  global_fields=reshape(buffer(n+1:n+global_size),shape(global_fields))
  global_metrics=buffer(n+global_size+1:);global_held=h%last_time
  if(any(global_fields(:,:,2)<=0).or.any(global_fields(:,:,3)<0))return''')
        return s
    edit('lfm_cube_ionosphere.f90',ion)
    def main(s):
        s=once(s,' use mpi_f08',' use cube_event_global,only:global_write\n use mpi_f08')
        # Both old and compact event drivers pass through the regular compact IO selector.
        from select_output import with_compact_output
        s=with_compact_output(s)
        if s.count('call archive%frame(c)')!=2: raise ValueError('Unexpected output paths')
        return s.replace('call archive%frame(c)','if(c%rank==0)call global_write(trim(path),c%step,c%time)\n   call archive%frame(c)').replace('relative_residual,','global_final_residual,').replace('SER_N,SER_S,ITER_N,ITER_S,ITMAX_N,ITMAX_S,converged_N,converged_S','global_residual_copy_N,global_residual_copy_S,refinements_copy_N,refinements_copy_S,budget_copy_N,budget_copy_S,converged_copy_N,converged_copy_S')
    edit('cube_lfm.f90',main)
