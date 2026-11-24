 ! Native GINT line averages of the untapered IGRF and its vector potential.
 ! Only the induction reconstruction uses them. CT state is initialized separately.
 subroutine induction_reference_edge(left,right,time,line,field)
  use lfm_cube_yy_background,only:cube_yy_evaluate_cgs
  real(real64),intent(in)::left(3),right(3),time
  real(real64),intent(out)::line,field(3)
  real(real64),parameter::nodes(6)=[0.1252334085_real64,0.3678314989_real64,0.5873179542_real64, &
      0.7699026741_real64,0.9041172563_real64,0.9815606342_real64]
  real(real64),parameter::weights(6)=[0.2491470458_real64,0.2334925365_real64,0.2031674267_real64, &
      0.1600783285_real64,0.1069393259_real64,0.0471753363_real64]
  real(real64)::mid(3),half(3),pa(3),ma(3),pb(3),mb(3),potential(3)
  integer::m,st
  mid=(left+right)/2;half=(right-left)/2
  field=0;potential=0
  do m=1,6
   call cube_yy_evaluate_cgs(mid+nodes(m)*half,time,0,0,pb,pa,st)
   call demand(st==0,'positive GINT abscissa in uncut induction reference')
   call cube_yy_evaluate_cgs(mid-nodes(m)*half,time,0,0,mb,ma,st)
   call demand(st==0,'negative GINT abscissa in uncut induction reference')
   field=field+weights(m)*(pb+mb)
   potential=potential+weights(m)*(pa+ma)
  enddo
  field=field/2
  line=dot_product(half,potential)
 end subroutine
