module Land_Eck
!Purpose: 
! 1) boost energy-momentum tensor Tmn to Landau frame and to Eckart frame
! 2) find the 4-velocity of such boosts
! 3) implement test procedure, which reads list of particle momenta from file,
!    builds energy-momentum tensor Tmn=\sum p^mu p^nu/p^0 and boosts it to Landau and Eckart frames
!    and prints results
!
! Uses:
!  subroutine rg from EISPACK free library to solve eigenvalue problem
!
!Author:
!  Dmytro Oliinychenko (oliiny@fias.uni-frankfurt.de)

contains

subroutine TestLandEck(ftest)
double precision p(4),Tmn(4,4), umu(4), Tmn_b(4,4), j(4)
integer npart, i, mu, nu
 character(LEN = *) ftest

Tmn=0.d0; j=0.d0;

!---Read-test-set-of-particles-and-form-tensor---
print *,"RUNNING LANDAU ECKART TESTING ROUTINE"
print *,"Reading particles from file ", ftest
print *

 open(unit=8, file=trim(adjustl(ftest)))
  read(8,*)npart
  do i=1,npart
   read(8,*)(p(mu),mu=1,4)
   do mu=1,4; do nu=1,4; Tmn(mu,nu) = Tmn(mu,nu) + p(mu)*p(nu)/p(1); end do; end do
   do mu=1,4; j(mu) = j(mu) + p(mu)/p(1); end do
  end do
 close(8)

!----Print--tensor-------------------------------
print *,"Energy-momentum tensor"
do mu=1,4; write(*,'(4f15.8)')(Tmn(mu,nu),nu=1,4); end do
print *,"tensor determinant: ",m44det(Tmn)
print *

!---Landau--frame--------------------------------
call FindLandau(Tmn,Tmn_b,umu)
print *,"Tensor in Landau frame"
do mu=1,4; write(*,'(4f15.8)')(Tmn_b(mu,nu),nu=1,4); end do
print *,"tensor determinant: ",m44det(Tmn_b)
print *,"Frame u_{mu}: ",(umu(mu),mu=1,4)
print *


!---Eckart--frame--------------------------------
call FindEckart(Tmn,j,Tmn_b,umu)
print *,"Tensor in Eckart frame"
do mu=1,4; write(*,'(4f15.8)')(Tmn_b(mu,nu),nu=1,4); end do
print *,"tensor determinant: ",m44det(Tmn_b)
print *,"Frame u_{mu}: ",(umu(mu),mu=1,4)
print *


end subroutine


!==============FindEckart===============================================
subroutine FindEckart(T,j,T_b,u)
!Given energy-momentum tensor T and conserved charge flow j subroutine finds:
! j is supposed to be with upper index
! 1) T_b - T boosted to Eckart frame
! 2) u - boost 4-velocity to Eckart frame, u^mu=j^mu/sqrt(j^nu j_nu)
implicit none
 double precision, intent(in) :: T(4,4), j(4)
 double precision, intent(out) :: T_b(4,4),u(4)

 call NormalizeTo4Vel(j,u, "FindEckart")
 u(2:4) = - u(2:4) !lower the index
 call BoostTensor(u,T,T_b) 
end subroutine

!==================FindLandau===========================================
subroutine FindLandau(T,T_b,u)
!Given energy-momentum tensor T subroutine finds:
! 1) T_b - T boosted to Landau frame
! 2) u - boost 4-velocity to make T^{0i} = 0 <=> boost to Landau frame 
implicit none

 double precision, intent(in) :: T(4,4)
 double precision, intent(out) :: T_b(4,4),u(4)
 double precision gmn(4,4),hlp(4,4),eigRe(4),eigIm(4),eiv(4,4),maxv
 integer mu, ierr, nu, nmax
 logical low_acc

 gmn=0.d0; gmn(1,1)=1.d0; do mu=2,4; gmn(mu,mu)=-1.d0; end do
 hlp=matmul(gmn,T)

 call rg(4, hlp, eigRe, eigIm, 1, eiv, ierr)
 
 !Catch errors
 if (ierr.ne.0) then; print *,"FindLandau: ierr: ",ierr; call exit(); endif
 do mu=1,4; if (abs(eigIm(mu))>1.d-14) then
   print *,"FindLandau: imaginary eigenvalue"; 
   print *,"eigRe: ",eigRe(mu)
   print *,"eigIm: ",eigIm(mu)
   print *,"eigen vector:";  print *,(eiv(nu,mu),nu=1,4)
   print *,"Tmn:"
   print *,T(1:4,1)
   print *,T(1:4,2)
   print *,T(1:4,3)
   print *,T(1:4,4)
   call exit(); 
 endif; end do

 !The largest eigenvalue (and the only positive) is energy density
 !Corresponding eigenvector is parallel to boost 4-velocity
 maxv=-1.d-12; do mu=1,4; if (eigRe(mu)>maxv) then; maxv=eigRe(mu); nmax=mu; endif; end do
 call NormalizeTo4Vel(eiv(1:4,nmax),u,"FindLandau")
 call BoostTensor(u,T,T_b)

 !Check that after boost T_b 0i components are computationally zero, i.e.
 !12 orders of magnitude less than max eigenvalue
 low_acc=.false.; do mu = 2,4;  if (abs(T_b(1,mu)/maxv) > 1.d-10) then; low_acc=.true.; endif; end do
 if (low_acc) then
   print *,"FindLandau: accuracy less than 10^-10";
   print *,"Tensor to transform:"
   do mu=1,4; print *,(T(mu,nu),nu=1,4); end do
   print *,"Tensor in Landau frame:"
   do mu=1,4; print *,(T_b(mu,nu),nu=1,4); end do
 endif

end subroutine


!===========================================================================
subroutine BoostTensor(v,T,T_b)
!boosts tensor T to 4-vector v(with lower index), returns boosted tensor T_b
!g is diagonal of metric tensor
implicit none

 double precision, intent(in) :: v(0:3)
 double precision, intent(in) :: T(0:3,0:3)
 double precision, intent(out) :: T_b(0:3,0:3)

 double precision u(0:3), M(0:3,0:3)

 call NormalizeTo4Vel(v,u, "BoostTensor")
 call GetBoostMatrix(u, M)
 T_b=matmul(M,matmul(T,M))

end subroutine

!===========================================================================
subroutine GetBoostMatrix(u_boost, M_boost)
!Gets boost matrix to 4-vector u_boost (lower index)
implicit none
double precision, intent(in) :: u_boost(0:3)
double precision, intent(out) :: M_boost(0:3,0:3)
integer ic,jc

 M_boost(0,0) = u_boost(0) !gamma

 do ic=1,3; M_boost(0,ic)=u_boost(ic); M_boost(ic,0)=M_boost(0,ic); end do
 do ic=1,3; do jc=1,3
  M_boost(ic,jc) = u_boost(ic)*u_boost(jc)/(u_boost(0)+1.d0);
  if (ic==jc) then;  M_boost(ic,jc) = M_boost(ic,jc)+1.d0; endif
 end do; end do

end subroutine

!===========================================================================
subroutine NormalizeTo4Vel(v,u, my_err_msg)
!Creates 4-vector u: u||v, u_mu*u^mu=1
implicit none
 double precision, dimension(0:3), intent(in) :: v
 double precision, dimension(0:3), intent(out) :: u
 double precision norm
 integer mu
 character(LEN=*) my_err_msg

 norm=LorentzProduct(v,v)
 if (norm.le.0.d0) then
  print *,"Error in ", my_err_msg
  print *,"NormalizeTo4vel: norm = ",norm," v: ",(v(mu),mu=0,3) 
  call exit(); 
 endif
 norm=sqrt(norm)
 u=v/norm
 if (u(0)<0) then; u=-u; endif

end subroutine

!==============LORENTZ=PRODUCT==============================================
double precision function LorentzProduct(v1,v2) result (prod)
!Returns Lorentz product of two 4-vectors
implicit none
double precision, dimension(0:3), intent(in) :: v1,v2
prod = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
return
end function LorentzProduct

!==============EUCLID=PRODUCT===============================================
double precision function EuclidProduct(v1,v2) result (prod)
!Returns Euclid product of two 4-component vectors
implicit none
double precision, dimension(0:3), intent(in) :: v1,v2
prod = v1(0)*v2(0)+v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
return
end function EuclidProduct

!================4x4==MATRIX=DETERMINANT====================================
double precision function m44det (A) result (det)
!Retutns determinant of 4x4 matrix A
implicit none
double precision, dimension(4,4), intent(in)  :: A

det =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)- &
       A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
       A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)- &
       A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
       A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

return
end function m44det
!==========================================================================
end module Land_Eck
