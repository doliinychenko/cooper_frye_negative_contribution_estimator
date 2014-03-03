program mvtk
!Prints out vtk file for energy density in the local Landau rest frame of each cell
!Written by Dmytro Oliinychenko (oliiny@fias.uni-frankfurt.de)
use read_f14f15, only: gen_part
use Land_Eck, only: FindLandau

implicit none

 !files and strings
 double precision, parameter :: E_collision = 40.d0
 double precision, parameter :: b_collision = 0.d0
 character(LEN=*), parameter :: input_prename = "/tmp/Tmn_proj_data/TrajLines/"

 !Grid dimensions
 integer, parameter :: nx = 10
 integer, parameter :: ny = 10
 integer, parameter :: nz = 80
 integer tstart_step, tend_step

 !Cell spacings
 double precision, parameter :: dx = 1.0d0
 double precision, parameter :: dy = 1.0d0
 double precision, parameter :: dz = 0.3d0
 double precision, parameter :: dt = 0.3d0

 !Particles and events
 type (gen_part) :: part
 integer Npart_evt, Npart_file, nevt


 double precision, dimension(:,:,:,:,:), allocatable :: Tmn
 double precision, allocatable :: EdensL(:,:,:,:) !Energy density in Landau frame


 integer it,ix,iy,iz, io, mu, nu, i, fnum
 character*100 s_hlp, t_option, input_fname
 double precision Tmn_hlp(0:3,0:3), TmnL_hlp(0:3, 0:3), umu(0:3), r1(0:3)
 logical op

 call getarg(1, t_option)
 read(t_option,*)tstart_step
 print *,"Start timestep: ",tstart_step, "      Start time:", tstart_step*dt, " fm/c"
 call getarg(2, t_option)
 read(t_option,*)tend_step
 print *,"End timestep: ",tend_step, "      End time:", tend_step*dt, " fm/c"


 allocate(Tmn(0:9, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(EdensL(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))


 fnum = 0
 nevt = 0

 do !Input files cycle 

   Npart_file = 0
   call make_in_name(input_prename,E_collision,b_collision,fnum,input_fname)
   inquire(file=input_fname, exist = op)
   if (.not. op) then; exit; endif
   
   print *,"Reading file: ",input_fname
   open(unit=8, form="unformatted", access="sequential", file = input_fname)
   do !event cycle
     read(8, iostat = io)Npart_evt
     if (io .ne. 0) then; exit; endif
     nevt = nevt + 1
     Npart_file = Npart_file + Npart_evt

     do i=1,Npart_evt
       read(8)part

       !Add particle to Tmn, jB, jS, etc
       do it = tstart_step, tend_step + 1
          if (part%ri(0) > it*dt .or. part%rf(0) < it*dt) then; cycle; endif

          call GetTrajPoint(it*dt, r1, part%ri, part%rf)
          ix = nint(r1(1)/dx); iy = nint(r1(2)/dy); iz = nint(r1(3)/dz)
          if (abs(ix) > nx .or. abs(iy) > ny .or. abs(iz) > nz) then; cycle; endif

          do mu=0,3
           do nu=mu,3
            Tmn(nu + 3*mu - mu*(mu-1)/2 ,it,ix,iy,iz) = Tmn(nu + 3*mu - mu*(mu-1)/2 ,it,ix,iy,iz) + &
                                                        part%p(mu)*part%p(nu)/part%p(0)
           end do
          end do


      
       end do 


     end do !particles in event cycle
   end do ! end event cycle
   close(8)
   print *,Npart_file, " lines read from file ", input_fname
   !print *,"ri(Npart), rf(Npart): ", part(395)%ri(0:3), part(395)%rf(0:3)
   fnum = fnum + 1  

 end do !End input files cycle

 print *,"Total ", nevt, " events read from ", fnum, " files."

 Tmn = Tmn/(dx*dy*dz*nevt)

!======================================================================================================
! Get Landau rest-frame energy density on the grid
!======================================================================================================
 print *,"Getting Landau rest frame energy densities."

 do it = tstart_step, tend_step + 1
  do ix=-nx,nx; do iy = -ny,ny; do iz = -nz,nz

    if (Tmn(0,it,ix,iy,iz) < 1.d-2) then; cycle; endif
    Tmn_hlp(0:3,0) = Tmn(0:3,it,ix,iy,iz); Tmn_hlp(0,0:3) = Tmn(0:3,it,ix,iy,iz)
    Tmn_hlp(1:3,1) = Tmn(4:6,it,ix,iy,iz); Tmn_hlp(1,1:3) = Tmn(4:6,it,ix,iy,iz)
    Tmn_hlp(2:3,2) = Tmn(7:8,it,ix,iy,iz); Tmn_hlp(2,2:3) = Tmn(7:8,it,ix,iy,iz)
    Tmn_hlp(3,3)   = Tmn(  9,it,ix,iy,iz)

    call FindLandau(Tmn_hlp,TmnL_hlp, umu(0:3))
    EdensL(it,ix,iy,iz) = TmnL_hlp(0,0)
  end do; end do; end do
 end do


!--------------------WRITE VTK FILE---------------------------------------------------------------
print *,"Writing vtk."
do it=tstart_step, tend_step + 1
 write(s_hlp,*)it; 
 open(unit=8, file="output/E40.0_b0.0prod/vtk/edens.vtk."//trim(adjustl(s_hlp)) )
  write(8,'(A)')"# vtk DataFile Version 2.0"
  write(8,'(A)')"energy density"
  write(8,'(A)')"ASCII"
  write(8,'(A)')"DATASET STRUCTURED_POINTS"
  write(8,'(A11,3I5)')"DIMENSIONS ", 2*nx+1, 2*ny+1, 2*nz+1
  write(8,'(A8,3I5)')"SPACING ", 1, 1, 1
  write(8,'(A7,3I5)')"ORIGIN ",-nx,-ny,-nz
  write(8,*)"POINT_DATA", (2*nx+1)*(2*ny+1)*(2*nz+1)
  write(8,*)"SCALARS energy_density float 1"
  write(8,*)"LOOKUP_TABLE default"

  do iz=-nz,nz; do iy=-ny,ny;
   write(8,'(200f10.3)')(EdensL(it,ix,iy,iz), ix=-nx,nx)
  end do; end do
 
 close(8)
end do

contains


!===============================================================================================
subroutine GetTrajPoint(t, r, ri,rf)
! ri(0:3) - initial point of a straight line, rf(0:3) - final point of straight line
!Subroutine finds point of straight line in moment t, ri(0)<t<rf(0)
implicit none
 double precision, intent(in) :: t, ri(0:3), rf(0:3)
 double precision, intent(out) :: r(0:3)
 double precision v(1:3)
 integer i
 
 if (ri(0) > t .or. rf(0) < t) then
   print *," - GetTrajPoint: t out of line range. ti,tf,t: ", ri(0), rf(0), t 
   call exit()
 endif
 do i=1,3;  v(i) = (rf(i) - ri(i))/(rf(0) - ri(0));  end do
 r(0) = t
 do i=1,3; r(i) = ri(i) + v(i) * (r(0) - ri(0)); end do

end subroutine


!===============================================================================================
subroutine make_in_name(prename, E, b, cnt, in_name)
implicit none
 character(LEN=*), intent(in) :: prename 
 double precision, intent(in) :: E, b     !Collision energy[GeV], centrality[fm]
 integer, intent(in) :: cnt                !file marker
 character(LEN=*), intent(out) :: in_name 
 character(LEN=10) :: Estr, bstr, nstr

 write(Estr,'(F10.1)')E
 write(bstr,'(F10.1)')b
 write(nstr,'(I10)')cnt

 in_name = prename//"E"//trim(adjustl(Estr))//"_b"//trim(adjustl(bstr))//"_n"//trim(adjustl(nstr))//".dat"

end subroutine

end program mvtk
