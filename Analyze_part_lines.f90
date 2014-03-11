program analyze_part_lines

 use read_f14f15, only: gen_part
 use Land_Eck, only: FindLandau, LorentzProduct, EuclidProduct
 use Histogram_module, only: histo, HistNew, HistAdd, BinToValue, HistPrint, HistBinSize
 use T_from_EoS, only: read_eos_tables, GetTemp, GetMuB, GetMuS, E_hadgas, B_hadgas, S_hadgas, pressure_hadgas

 implicit none

 type part_equil
  sequence
  character(len=10)pname                !name
  !In UrQMD particle is identified by two numbers: ityp and i3
  integer ityp
  integer Ni3                            !number of i3's 
  integer, dimension(10) :: i3          !list of 2*i3 
  double precision m                    !mass in GeV
  integer g                              !degeneracy factor
  integer B                              !baryon charge
  integer S                              !strangeness
 end type part_equil

 !Particles of interest arrays: y and pt spectra will be calculated for these particles
 integer, parameter :: max_plot_part = 8
 type (part_equil), dimension(max_plot_part) :: pequil
 integer Npequil

 !Grid dimensions
 integer, parameter :: nx = 20
 integer, parameter :: ny = 20
 integer, parameter :: nz = 40
 integer tstart_step, tend_step

 !Cell spacings
 double precision, parameter :: dx = 1.0d0
 double precision, parameter :: dy = 1.0d0
 double precision, parameter :: dz = 1.0d0
 double precision, parameter :: dt = 1.0d0

 !Pressure isotropy criterion
 double precision, parameter :: p_iso_crit = 0.1d0 ! 2*|p_a - p_b|/(p_a+p_b) < p_iso_crit - isotropy

 !files and strings
 double precision, parameter :: E_collision = 40.d0
 double precision, parameter :: b_collision = 0.d0
 character(LEN=*), parameter :: input_prename = "/tmp/Tmn_proj_data/TrajLines/"
 character(LEN=*), parameter :: dir_stamp = "thermalization" !stamp for output directory
 character(LEN=100) :: out_dir
 character(LEN=100) :: s_hlp, input_fname, t_option
 logical op

 !Particles and events
 type (gen_part) :: part
 integer Npart_evt, Npart_file, nevt

 !Arrays on grid
 double precision, dimension(:,:,:,:,:), allocatable :: Tmn, TmnL
 double precision, dimension(:,:,:,:), allocatable :: Temp, MuB, MuS
 logical, dimension(:,:,:,:), allocatable :: pres_isotr
 double precision, dimension(:,:), allocatable :: mom_tot
 double precision, dimension(:), allocatable :: B_tot, S_tot, I3_tot

 
 integer fnum, mu, nu, io
 integer i, it,ix,iy,iz
 double precision r1(0:3), Tmn_hlp(0:3,0:3), TmnL_hlp(0:3, 0:3), umu_hlp(0:3), presx, presy, presz, pcrit


 !===========================================================================================
 
 !Get tstart_step and tend_step from command line. Allocate arrays.

 call getarg(1, t_option)
 read(t_option,*)tstart_step
 print *,"Start timestep: ",tstart_step, "      Start time:", tstart_step*dt, " fm/c"
 call getarg(2, t_option)
 read(t_option,*)tend_step
 print *,"End timestep: ",tend_step, "      End time:", tend_step*dt, " fm/c"

 allocate(Tmn(0:9, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(TmnL(0:9, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))

 allocate(Temp(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(MuB(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(MuS(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))

 allocate(pres_isotr(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))

 allocate(mom_tot(0:3,tstart_step:tend_step+1))
 allocate(B_tot(tstart_step:tend_step+1))
 allocate(S_tot(tstart_step:tend_step+1))
 allocate(I3_tot(tstart_step:tend_step+1))

 Tmn      = 0.d0
 TmnL     = 0.d0
 pres_isotr = .FALSE.

 call create_output_folder_structure(E_collision, b_collision, dir_stamp, out_dir)
 print *,"Output directory is: ", trim(out_dir)

 call allhist_init()

 !Read input particle files and find Tmn 

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

          mom_tot(0:3,it) = mom_tot(0:3,it) + part%p(0:3)
          B_tot(it) = B_tot(it) + part%B*1.d0
          S_tot(it) = S_tot(it) + part%S*1.d0
          I3_tot(it) = I3_tot(it) + part%I3*0.5d0

          call GetTrajPoint(it*dt, r1, part%ri, part%rf)
          ix = nint(r1(1)/dx); iy = nint(r1(2)/dy); iz = nint(r1(3)/dz)
          if (abs(ix) > nx .or. abs(iy) > ny .or. abs(iz) > nz) then; cycle; endif
          
          !Tmn: 00, 01, 02, 03, 11, 12, 13, 22, 23, 33
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

 mom_tot = mom_tot/nevt
 B_tot   =   B_tot/nevt
 S_tot   =   S_tot/nevt
 I3_tot  =  I3_tot/nevt

 open(unit=8,file="output/"//trim(out_dir)//"/total_cons_quantities.txt")
   do it= tstart_step, tend_step + 1
    write(8,*)it*dt, mom_tot(0:3,it), B_tot(it), S_tot(it), I3_tot(it)
   end do
 close(8)

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

    call FindLandau(Tmn_hlp,TmnL_hlp, umu_hlp(0:3))
    TmnL(0:3,it,ix,iy,iz) = TmnL_hlp(0,0:3)
    TmnL(4:6,it,ix,iy,iz) = TmnL_hlp(1,1:3)
    TmnL(7:8,it,ix,iy,iz) = TmnL_hlp(2,2:3)
    TmnL(9:9,it,ix,iy,iz) = TmnL_hlp(3,3:3)

    presx = TmnL_hlp(1,1)
    presy = TmnL_hlp(2,2)
    presz = TmnL_hlp(3,3)

    if ( abs(presx-presy)*2.d0/(presx+presy) < p_iso_crit .AND. &
         abs(presx-presz)*2.d0/(presx+presz) < p_iso_crit .AND. &
         abs(presy-presz)*2.d0/(presy+presz) < p_iso_crit ) then
     pres_isotr(it,ix,iy,iz) = .TRUE.
    endif

  end do; end do; end do
 end do

!======================================================================================================
!  Analysis of thermalization
!======================================================================================================

 
!======================================================================================================
!  Finalize and output
!======================================================================================================

 open(unit=8, file = "output/"//trim(out_dir)//"/isotropy_cells.txt")
  do it = tstart_step, tend_step;
   do ix=-nx,nx; do iy = -ny,ny; do iz = -nz,nz
    if (pres_isotr(it,ix,iy,iz)) then
     write(8,*)it*dt, ix*dx, iy*dy, iz*dz
    endif
   end do; end do; end do
  end do
 close(8)

 print *,"Writing vtk: isotropic/not isotropic pressure in cells"
 do it=tstart_step, tend_step + 1
  write(s_hlp,*)it; 
  open(unit=8, file="output/"//trim(out_dir)//"/vtk/p_iso.vtk."//trim(adjustl(s_hlp)) )
   write(8,'(A)')"# vtk DataFile Version 2.0"
   write(8,'(A)')"pressure isotropy"
   write(8,'(A)')"ASCII"
   write(8,'(A)')"DATASET STRUCTURED_POINTS"
   write(8,'(A11,3I5)')"DIMENSIONS ", 2*nx+1, 2*ny+1, 2*nz+1
   write(8,'(A8,3I5)')"SPACING ", 1, 1, 1
   write(8,'(A7,3I5)')"ORIGIN ",-nx,-ny,-nz
   write(8,'(A,X,I10)')"POINT_DATA", (2*nx+1)*(2*ny+1)*(2*nz+1)
   write(8,'(A)')"SCALARS pressure_isotropy float 1"
   write(8,'(A)')"LOOKUP_TABLE default" 

   do iz=-nz,nz; do iy=-ny,ny; do ix=-nx,nx
    presx = TmnL(4,it,ix,iy,iz)
    presy = TmnL(7,it,ix,iy,iz)
    presz = TmnL(9,it,ix,iy,iz)
    if (presx > 0.d0 .AND. presy > 0.d0 .AND. presz>0.d0) then
     pcrit = (abs(presx-presy)+abs(presy-presz)+abs(presz-presx))/(presx+presy+presz)
    else
     pcrit = 2.d0
    endif

    write(8,'(f10.3)') pcrit
    
   end do; end do; end do
  close(8)
 end do


 print *,"Writing vtk: energy density in cells"
 do it=tstart_step, tend_step + 1
  write(s_hlp,*)it; 
  open(unit=8, file="output/"//trim(out_dir)//"/vtk/edens.vtk."//trim(adjustl(s_hlp)) )
   write(8,'(A)')"# vtk DataFile Version 2.0"
   write(8,'(A)')"energy density"
   write(8,'(A)')"ASCII"
   write(8,'(A)')"DATASET STRUCTURED_POINTS"
   write(8,'(A11,3I5)')"DIMENSIONS ", 2*nx+1, 2*ny+1, 2*nz+1
   write(8,'(A8,3I5)')"SPACING ", 1, 1, 1
   write(8,'(A7,3I5)')"ORIGIN ",-nx,-ny,-nz
   write(8,'(A,X,I10)')"POINT_DATA", (2*nx+1)*(2*ny+1)*(2*nz+1)
   write(8,'(A)')"SCALARS energy_density float 1"
   write(8,'(A)')"LOOKUP_TABLE default" 

  do iz=-nz,nz; do iy=-ny,ny;
   write(8,'(200f10.3)')(TmnL(0,it,ix,iy,iz), ix=-nx,nx)
  end do; end do
  close(8)
 end do


 call allhist_finalize_and_output(out_dir)


 deallocate(Tmn)
 deallocate(TmnL)
 deallocate(Temp)
 deallocate(MuB)
 deallocate(MuS)
 deallocate(pres_isotr)
 deallocate(mom_tot)
 deallocate(B_tot)
 deallocate(S_tot)
 deallocate(I3_tot)

!======================================END=============================================================

contains
!======================================================================================================
subroutine create_output_folder_structure(E,b,stamp, run_folder)
 implicit none
 double precision E,b
 character(len=*) stamp
 character(len=*) run_folder
 character(len=10) :: Estr, bstr

 write(Estr,'(F10.1)')E
 write(bstr,'(F10.1)')b

 run_folder = "E"//trim(adjustl(Estr))//"_b"//trim(adjustl(bstr))//trim(adjustl(stamp))
 call system('mkdir -p ./output')
 call system('mkdir -p ./output/'//trim(run_folder))
 call system('mkdir -p ./output/'//trim(run_folder)//"/vtk")
 
end subroutine


!===========================================================================================================
subroutine allhist_init()
implicit none

end subroutine


!===============================================================================================
double precision function GetY(p) result (y)
implicit none
double precision p(0:3)
 y = 0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
end function GetY

!===============================================================================================
double precision function GetPt(p) result (pt)
implicit none
double precision p(0:3)
 pt = sqrt(p(1)**2 + p(2)**2)
end function GetPt

!===========================================================================================================
subroutine allhist_finalize_and_output(dname)
implicit none
 !integer i,j
 character(len=*)dname

end subroutine



!===========================================================================================================
subroutine pequil_init()
implicit none
  integer i

  Npequil = 6
  pequil(1)%pname = "pi";  pequil(1)%m = 0.138d0; pequil(1)%g =  3; pequil(1)%B = 0; pequil(1)%S =  0
  pequil(1)%ityp = 101; pequil(1)%Ni3 = 3; pequil(1)%i3(1:3) = (/-2,0,2/)
  pequil(2)%pname = "N";   pequil(2)%m = 0.938d0; pequil(2)%g =  4; pequil(2)%B = 1; pequil(2)%S =  0
  pequil(2)%ityp =   1; pequil(2)%Ni3 = 2; pequil(2)%i3(1:2) = (/-1,1/)
  pequil(3)%pname = "De";  pequil(3)%m = 1.232d0; pequil(3)%g = 16; pequil(3)%B = 1; pequil(3)%S =  0
  pequil(3)%ityp =  17; pequil(3)%Ni3 = 4; pequil(3)%i3(1:4) = (/-3,-1,1,3/)
  pequil(4)%pname = "La";  pequil(4)%m = 1.115d0; pequil(4)%g =  1; pequil(4)%B = 1; pequil(4)%S = -1
  pequil(4)%ityp =  27; pequil(4)%Ni3 = 1; pequil(4)%i3(1:1) = (/0/)
  pequil(5)%pname = "Kmi"; pequil(5)%m = 0.495d0; pequil(5)%g =  1; pequil(5)%B = 0; pequil(5)%S = -1
  pequil(5)%ityp =-106; pequil(5)%Ni3 = 1; pequil(5)%i3(1:1) = (/-1/)
  pequil(6)%pname = "Kpl"; pequil(6)%m = 0.495d0; pequil(6)%g =  1; pequil(6)%B = 0; pequil(6)%S =  1
  pequil(6)%ityp = 106; pequil(6)%Ni3 = 1; pequil(6)%i3(1:1) = (/1/)

  print *,"=============================================================="
  print *,"Initializing particles for comparing N(T,mu) and N(UrQMD):    "
  print *,"=============================================================="

  write(*,'(A4,2X,A5,2X,3A3,2X,A5,2X,A)')"name"," mass","  g", "  B", "  S", " ityp"," 2*i3 list"
  do i=1, Npequil
   write(*,'(A4,2X,F5.3,2X,3I3,2X,I5,2X,6I3)')pequil(i)%pname, pequil(i)%m, pequil(i)%g, pequil(i)%B, pequil(i)%S,&
                                              pequil(i)%ityp, pequil(i)%i3( 1:pequil(i)%Ni3 )
  end do

end subroutine


!===========================================================================================================
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


end program analyze_part_lines
