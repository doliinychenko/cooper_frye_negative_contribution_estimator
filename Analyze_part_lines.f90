program read_sorted_lines

 use read_f14f15, only: gen_part
 use Land_Eck, only: FindLandau, LorentzProduct, EuclidProduct

 implicit none

 !Grid dimensions
 integer, parameter :: nx = 35 
 integer, parameter :: ny = 35
 integer, parameter :: nz = 60
 integer tstart_step, tend_step

 !Cell spacings
 double precision, parameter :: dx = 0.3d0
 double precision, parameter :: dy = 0.3d0
 double precision, parameter :: dz = 0.3d0
 double precision, parameter :: dt = 0.3d0

 !files and strings
 character(LEN=*), parameter :: urqmdf14_input = "/tmp/Tmn_proj_data/UrQMD_f14f15/E40b0ev4500.f14"
 character(LEN=*), parameter :: urqmdf15_input = "/tmp/Tmn_proj_data/UrQMD_f14f15/E40b0ev4500.f15"
 character(LEN=*), parameter :: input_prename = "/tmp/Tmn_proj_data/TrajLines/"
 character(LEN=100) :: input_fname, t_option
 logical op

 !For hypersurface finding via Cornelius
 double precision, parameter :: e0 = 0.15d0 !GeV
 double precision,dimension(0:1,0:1,0:1,0:1) :: HyperCube      
 double precision,dimension(0:3,8) :: dSigma
 integer        :: Nsurf
 double precision,dimension(0:3,8) :: Vmid
 integer      :: Nambi, Ndisc
 
 !Particles and events
 type (gen_part) :: part
 integer Npart_evt, Npart_file, nevt

 !Arrays on grid
 double precision, dimension(:,:,:,:,:), allocatable :: Tmn, jB, jS, umu
 double precision, dimension(:,:,:,:), allocatable :: EdensL
 double precision, dimension(:,:), allocatable :: mom_tot
 double precision, dimension(:), allocatable :: B_tot, S_tot, I3_tot
 double precision, dimension(:), allocatable :: E_fl
 double precision, dimension(:), allocatable :: E_in

 
 integer fnum, i, mu, nu, io
 integer it,ix,iy,iz
 double precision r1(0:3), r2(0:3), Tmn_hlp(0:3,0:3), TmnL_hlp(0:3, 0:3), ar


 !===========================================================================================
 
 !Get tstart_step and tend_step from command line. Allocate arrays.

 call getarg(1, t_option)
 read(t_option,*)tstart_step
 print *,"Start timestep: ",tstart_step
 call getarg(2, t_option)
 read(t_option,*)tend_step
 print *,"End timestep: ",tend_step

 allocate(Tmn(0:9, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(jB(0:3, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(jS(0:3, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(umu(0:3, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(EdensL(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(mom_tot(0:3,tstart_step:tend_step+1))
 allocate(B_tot(tstart_step:tend_step+1))
 allocate(E_fl(tstart_step:tend_step+1))
 allocate(E_in(tstart_step:tend_step+1))
 allocate(S_tot(tstart_step:tend_step+1))
 allocate(I3_tot(tstart_step:tend_step+1))
 Tmn     = 0.d0
 jB      = 0.d0
 jS      = 0.d0
 mom_tot = 0.d0
 B_tot   = 0.d0 
 S_tot   = 0.d0 
 I3_tot  = 0.d0 
 EdensL  = 0.d0
 umu(0,:,:,:,:) = 1.d0;  umu(1:3,:,:,:,:) = 0.d0


 !Read input particle files and find Tmn 

 fnum = 0
 nevt = 0

 do !Input files cycle 

   Npart_file = 0
   call make_in_name(input_prename,40.d0,0.d0,fnum,input_fname)
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

          do mu=0,3
           do nu=mu,3
            Tmn(nu + 3*mu - mu*(mu-1)/2 ,it,ix,iy,iz) = Tmn(nu + 3*mu - mu*(mu-1)/2 ,it,ix,iy,iz) + &
                                                        part%p(mu)*part%p(nu)/part%p(0)
           end do
          end do


          jB(0:3,it,ix,iy,iz) =  jB(0:3,it,ix,iy,iz) + part%p(0:3)/part%p(0)*part%B
          jS(0:3,it,ix,iy,iz) =  jS(0:3,it,ix,iy,iz) + part%p(0:3)/part%p(0)*part%S
      
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
 jB  =  jB/(dx*dy*dz*nevt)
 jS  =  jS/(dx*dy*dz*nevt)

 mom_tot = mom_tot/nevt
 B_tot   =   B_tot/nevt
 S_tot   =   S_tot/nevt
 I3_tot  =  I3_tot/nevt

 open(unit=8,file="output/total_quantities_t1_t2.txt")
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

    call FindLandau(Tmn_hlp,TmnL_hlp, umu(0:3,it,ix,iy,iz))
    EdensL(it,ix,iy,iz) = TmnL_hlp(0,0)
  end do; end do; end do
 end do


!======================================================================================================
!  Get hypersurface using Cornelius
!======================================================================================================
 print *,"HYDRO STYLE"

 E_fl = 0.d0

 do it = tstart_step, tend_step
  do ix=-nx,nx-1; do iy = -ny,ny-1; do iz = -nz,nz-1


   HyperCube=EdensL(it:it+1,ix:ix+1,iy:iy+1,iz:iz+1)
   
   if (minval(HyperCube) <= e0 .AND. maxval(HyperCube) >= e0)    then 

      call Cornelius(e0, HyperCube, dSigma, Nsurf, Vmid, dt, dx, dy, dz, Nambi, Ndisc) ! Let Cornelius do all the work
 
      do i=1,Nsurf
        
       do mu=0,3
        Tmn_hlp(mu,0) = Cube4Intpl(Vmid(0:3,i), Tmn(mu,it:it+1,ix:ix+1,iy:iy+1,iz:iz+1))
       end do    

       E_fl(it+1) = E_fl(it+1) + EuclidProduct(Tmn_hlp(0:3,0), dSigma(0:3,i))

      end do

   endif


  end do; end do; end do
 end do


!======================================================================================================
!  Get particle crossings by explicit counting
!======================================================================================================
 print *,"BY PARTICLES"
 E_in = 0.d0

 fnum = 0
 nevt = 0

 do !Input files cycle 

   Npart_file = 0
   call make_in_name(input_prename,40.d0,0.d0,fnum,input_fname)
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
       !Analyse here

       do it = tstart_step, tend_step
        if (part%ri(0) > it*dt .or. part%rf(0) < it*dt) then; cycle; endif
        call GetTrajPoint(it*dt, r1, part%ri, part%rf)
        ix = floor(r1(1)/dx)
        iy = floor(r1(2)/dy)
        iz = floor(r1(3)/dz)
        r1(0) = 0.d0;
        r1(1) = r1(1)/dx - ix*1d0
        r1(2) = r1(2)/dy - iy*1d0
        r1(3) = r1(3)/dz - iz*1d0

        if (ix < -nx .or. ix >= nx .or. iy < -ny .or. iy >= ny .or. iz < -nz .or. iz >= nz ) then
         ar=0.d0
        else
         ar = Cube4Intpl(r1,EdensL(it:it+1,ix:ix+1, iy:iy+1, iz:iz+1))
        endif

        if (ar .ge. e0) then
          E_in(it) = E_in(it) + part%p(0)
        endif
       end do 

     end do !particles in event cycle
   end do ! end event cycle
   close(8)
   print *,Npart_file, " lines read from file ", input_fname
   !print *,"ri(Npart), rf(Npart): ", part(395)%ri(0:3), part(395)%rf(0:3)
   fnum = fnum + 1  

 end do !End input files cycle
 print *,"Total ", nevt, " events read from ", fnum, " files."

 E_in = E_in/nevt

!======================================================================================================
!  Finalize and output
!======================================================================================================

 open(unit=8, file = "output/Econs.txt")
  do it = tstart_step, tend_step
   write(8,*)it*dt, E_fl(it+1), E_in(it) - E_in(it+1)
  end do
 close(8)


 deallocate(Tmn)
 deallocate(jB)
 deallocate(jS)
 deallocate(umu)
 deallocate(EdensL)
 deallocate(mom_tot)
 deallocate(B_tot)
 deallocate(S_tot)
 deallocate(I3_tot)

contains
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

!===============================================================================================
double precision function Cube4Intpl(v,CUBE) result (res)
!Given values on cubic 4d grid [0,1]x[0,1]x[0,1]x[0,1] linearly interpolates to v(0),v(1),v(2),v(3)
implicit none
 double precision, parameter :: u = 1.d0
 double precision, dimension(0:3), intent(in) :: v !each component in [0,u]
 double precision, dimension(0:1,0:1,0:1,0:1), intent(in) :: CUBE
 double precision c_weight
 integer i,j,k,l,mu
 
 res=0.d0
 do mu = 0,3
  if (v(mu) > 1.d0 .OR. v(mu) < 0.d0) then
   print *,"Warning: trying to interpolate to point: ", v(0:3)
  endif
 end do
 do i=0,1; do j=0,1; do k=0,1; do l=0,1
   c_weight = abs(  (u - i*u - v(0)) * (u - j*u - v(1)) * ( u - k*u - v(2)) * (u - l*u - v(3)) )
   res=res+c_weight*CUBE(i,j,k,l)
 end do; end do; end do; end do
 
end function Cube4Intpl

end program read_sorted_lines
