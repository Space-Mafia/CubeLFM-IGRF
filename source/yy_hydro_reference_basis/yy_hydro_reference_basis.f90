module yy_hydro_reference_basis
 use iso_fortran_env,only:real64
 implicit none
 private
 public::yy_reference_coefficients,yy_reference_at
 contains
 subroutine yy_reference_coefficients(samples,coef)
 real(real64),intent(in)::samples(:,:)
 real(real64),intent(out)::coef(:,:)
 real(real64)::phi,pi
 integer::n,m,k
 n=size(samples,2);pi=acos(-1._real64)
 if(mod(n,2)/=1.or.any(shape(coef)/=shape(samples)))error stop 'odd exact harmonic sample count'
 coef=0
 do k=1,n
 phi=2*pi*(k-1)/n
 coef(:,1)=coef(:,1)+samples(:,k)/n
 do m=1,(n-1)/2
 coef(:,2*m)=coef(:,2*m)+2*samples(:,k)*cos(m*phi)/n
 coef(:,2*m+1)=coef(:,2*m+1)+2*samples(:,k)*sin(m*phi)/n
 enddo
 enddo
 end subroutine
 subroutine yy_reference_at(coef,phase,values)
 real(real64),intent(in)::coef(:,:),phase
 real(real64),intent(out)::values(:)
 integer::m
 if(size(values)/=size(coef,1).or.mod(size(coef,2),2)/=1)error stop 'reference basis shape'
 values=coef(:,1)
 do m=1,(size(coef,2)-1)/2
 values=values+coef(:,2*m)*cos(m*phase)+coef(:,2*m+1)*sin(m*phase)
 enddo
 end subroutine
end module
