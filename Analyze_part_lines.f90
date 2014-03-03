program analyze_part_lines

 use read_f14f15, only: gen_part
 use Land_Eck, only: FindLandau, LorentzProduct, EuclidProduct
 use Histogram_module, only: histo, HistNew, HistAdd, BinToValue, HistPrint, HistBinSize
 use T_from_EoS, only: read_eos_tables, GetTemp, GetMuB, GetMuS, E_hadgas, B_hadgas, S_hadgas, pressure_hadgas
 use CF_int, only: dNdy, dNptdpt, init_CF_int

 implicit none

 type part_hyd
  sequence
  character(len=10)pname    !name
  double precision m        !mass in GeV
  integer g                  !degeneracy factor
  integer B                  !baryon charge
  integer S                  !strangeness
 end type part_hyd

 type part_prt
  sequence
  character(len=10)pname                !name
  !In UrQMD particle is identified by two numbers: ityp and i3
  integer ityp
  integer Ni3                            !number of i3's 
  integer, dimension(10) :: i3          !list of 2*i3 
 end type part_prt

 !Particles of interest arrays: y and pt spectra will be calculated for these particles
 integer, parameter :: max_plot_part = 8
 type (part_hyd), dimension(max_plot_part) :: phyd
 integer Nphyd
 type (part_prt), dimension(max_plot_part) :: pprt
 integer Npprt

 !Options
 logical, parameter :: calculate_hyd_spectra = .FALSE.
 logical, parameter :: calculate_prt_spectra = .FALSE.

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

 !Histograms
 type (histo) hyd_surf_sig(0:3), hyd_surf_Vmid(0:3), hyd_surf_e, hyd_surf_nb
 type (histo) hyd_surf_T, hyd_surf_mub, hyd_surf_mus, hyd_surf_v(1:3), hyd_surf_gam
 type (histo), dimension(max_plot_part) :: hyd_spect_y, hyd_spect_pt, hyd_neg_y, hyd_neg_pt
 type (histo), dimension(max_plot_part) :: prt_spect_y, prt_spect_pt, prt_neg_y, prt_neg_pt
 type (histo) pres_CF, pres_x, pres_y, pres_z 

 !files and strings
 double precision, parameter :: E_collision = 10.d0
 double precision, parameter :: b_collision = 0.d0
 character(LEN=*), parameter :: input_prename = "/tmp/Tmn_proj_data/TrajLines/"
 character(LEN=*), parameter :: dir_stamp = "testing" !stamp for output directory: testing, prod or "" are suggested
 character(LEN=100) :: out_dir
 character(LEN=100) :: input_fname, t_option
 logical op

 !For hypersurface finding via Cornelius
 double precision, parameter :: e0 = 0.15d0 !GeV
 double precision,dimension(0:1,0:1,0:1,0:1) :: e_HC, nb_HC, T_HC, mub_HC, mus_HC
 double precision,dimension(0:3,8) :: dSigma
 integer        :: Nsurf
 double precision,dimension(0:3,8) :: Vmid
 integer      :: Nambi, Ndisc

 !Interpolated values on the surface
 double precision txyz_intpl(0:3), e_intpl, T_intpl, nb_intpl, muB_intpl, muS_intpl, v_intpl(1:3), u_intpl(0:3)
 double precision presx_intpl, presy_intpl, presz_intpl
 
 !Particles and events
 type (gen_part) :: part
 integer Npart_evt, Npart_file, nevt

 !Arrays on grid
 double precision, dimension(:,:,:,:,:), allocatable :: Tmn, jB, jS, umu
 double precision, dimension(:,:,:,:), allocatable :: EdensL, TL11, TL22, TL33
 double precision, dimension(:,:), allocatable :: mom_tot
 double precision, dimension(:), allocatable :: B_tot, S_tot, I3_tot
 double precision, dimension(:), allocatable :: E_fl, B_fl, S_fl, E_in, B_in, S_in, E_CF, B_CF, S_CF

 
 integer fnum, mu, nu, io, a,b,c,d
 integer i,j,k, it,ix,iy,iz
 double precision r1(0:3), r2(0:3), Tmn_hlp(0:3,0:3), TmnL_hlp(0:3, 0:3), jB_hlp(0:3), jS_hlp(0:3), ar, br
 double precision dt_ptraj, t_ptraj, umu_sigmu, bin_val


 !===========================================================================================
 
 !Get tstart_step and tend_step from command line. Allocate arrays.

 call getarg(1, t_option)
 read(t_option,*)tstart_step
 print *,"Start timestep: ",tstart_step, "      Start time:", tstart_step*dt, " fm/c"
 call getarg(2, t_option)
 read(t_option,*)tend_step
 print *,"End timestep: ",tend_step, "      End time:", tend_step*dt, " fm/c"

 allocate(Tmn(0:9, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(jB(0:3, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(jS(0:3, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(umu(0:3, tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(EdensL(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(TL11(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(TL22(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(TL33(tstart_step:tend_step+1, -nx:nx, -ny:ny, -nz:nz))
 allocate(E_fl(tstart_step:tend_step+1))
 allocate(E_in(tstart_step:tend_step+1))
 allocate(E_CF(tstart_step:tend_step+1))
 allocate(B_fl(tstart_step:tend_step+1))
 allocate(B_in(tstart_step:tend_step+1))
 allocate(B_CF(tstart_step:tend_step+1))
 allocate(S_fl(tstart_step:tend_step+1))
 allocate(S_in(tstart_step:tend_step+1))
 allocate(S_CF(tstart_step:tend_step+1))

 allocate(mom_tot(0:3,tstart_step:tend_step+1))
 allocate(B_tot(tstart_step:tend_step+1))
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

 call create_output_folder_structure(E_collision, b_collision, dir_stamp, out_dir)
 print *,"Output directory is: ", trim(out_dir)

 print *,"Hyper-surface criterion E0(GeV/fm^3) is: ", e0

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

    call FindLandau(Tmn_hlp,TmnL_hlp, umu(0:3,it,ix,iy,iz))
    EdensL(it,ix,iy,iz) = TmnL_hlp(0,0)
    TL11(it,ix,iy,iz) = Tmn_hlp(1,1)
    TL22(it,ix,iy,iz) = Tmn_hlp(2,2)
    TL33(it,ix,iy,iz) = Tmn_hlp(3,3)
  end do; end do; end do
 end do


!======================================================================================================
!  Get hypersurface using Cornelius
!======================================================================================================
 print *,"HYDRO STYLE"
 call phyd_init()
 print *,"Reading EoS tables"
 call read_eos_tables()
 print *,"Calculating coefficients and points for Gauss-Legendre integration."
 call init_CF_int()
 print *,"Starting hydro-style calculations: cycle over hypersurface."

 
 open(unit=8, file="output/"//trim(out_dir)//"/hyd/surfTX.txt")
 open(unit=9, file="output/"//trim(out_dir)//"/hyd/surfTZ.txt")

 E_fl = 0.d0; B_fl = 0.d0; S_fl = 0.d0
 E_CF = 0.d0; B_CF = 0.d0; S_CF = 0.d0

 do it = tstart_step, tend_step
  write(*,'(A8,I5,A)',advance="no")"timestep: ",it,char(13);
  do ix=-nx,nx-1; do iy = -ny,ny-1; do iz = -nz,nz-1


   e_HC=EdensL(it:it+1,ix:ix+1,iy:iy+1,iz:iz+1)
   
   if (minval(e_HC) <= e0 .AND. maxval(e_HC) >= e0)    then 

      call Cornelius(e0, e_HC, dSigma, Nsurf, Vmid, dt, dx, dy, dz, Nambi, Ndisc) ! find pieces of hypersurface
 
      do i=1,Nsurf

       if (iz==0 .and. iy==0) then; write(8,*)it*dt+Vmid(0,i), ix*dx+Vmid(1,i); endif   
       if (ix==0 .and. iy==0) then; write(9,*)it*dt+Vmid(0,i), iz*dz+Vmid(3,i); endif
        
       !Calculate flows through hypersurface
       do mu=0,3
        Tmn_hlp(mu,0) = Cube4Intpl(Vmid(0:3,i), Tmn(mu,it:it+1,ix:ix+1,iy:iy+1,iz:iz+1))
        jB_hlp(mu)    = Cube4Intpl(Vmid(0:3,i),  jB(mu,it:it+1,ix:ix+1,iy:iy+1,iz:iz+1))
        jS_hlp(mu)    = Cube4Intpl(Vmid(0:3,i),  jS(mu,it:it+1,ix:ix+1,iy:iy+1,iz:iz+1))
       end do    
       E_fl(it+1) = E_fl(it+1) + EuclidProduct(Tmn_hlp(0:3,0), dSigma(0:3,i))
       B_fl(it+1) = B_fl(it+1) + EuclidProduct(   jB_hlp(0:3), dSigma(0:3,i))
       S_fl(it+1) = S_fl(it+1) + EuclidProduct(   jS_hlp(0:3), dSigma(0:3,i))

       !Calculate e, nb, T muB, muS in the hypercell
       do a=0,1; do b=0,1; do c=0,1; do d=0,1
        nb_HC(a,b,c,d) = EuclidProduct (jB(0:3,it+a,ix+b,iy+c,iz+d), umu(0:3,it+a,ix+b,iy+c,iz+d))
        T_HC(a,b,c,d) = GetTemp(e_HC(a,b,c,d), nb_HC(a,b,c,d))
        mub_HC(a,b,c,d) = GetMuB(e_HC(a,b,c,d), nb_HC(a,b,c,d))
        muS_HC(a,b,c,d) = GetMuS(e_HC(a,b,c,d), nb_HC(a,b,c,d))
       end do; end do; end do; end do

       !Interpolate e, nb, T, muB, muS, vx, vy, vz, gamma to Vmid - to hypersurface
       txyz_intpl(0:3) = Vmid(0:3,i) / (/dt,dx,dy,dz/) !txyz_intpl - coordinates in hypercube, from 0 to 1
       e_intpl   =   Cube4Intpl(txyz_intpl, e_HC)
       
       nb_intpl  =   Cube4Intpl(txyz_intpl, nb_HC)
       T_intpl   =   Cube4Intpl(txyz_intpl, T_HC)*1.d-3  !MeV to GeV
       muB_intpl =   Cube4Intpl(txyz_intpl, mub_HC)*1.d-3 !MeV to GeV
       muS_intpl =   Cube4Intpl(txyz_intpl, mus_HC)*1.d-3 !MeV to GeV
       presx_intpl  =   Cube4Intpl(txyz_intpl, TL11(it:it+1, ix:ix+1, iy:iy+1, iz:iz+1))
       presy_intpl  =   Cube4Intpl(txyz_intpl, TL22(it:it+1, ix:ix+1, iy:iy+1, iz:iz+1))
       presz_intpl  =   Cube4Intpl(txyz_intpl, TL33(it:it+1, ix:ix+1, iy:iy+1, iz:iz+1))
       do mu = 1,3
        v_intpl(mu)  =  Cube4Intpl(txyz_intpl, umu(mu,it:it+1,ix:ix+1,iy:iy+1,iz:iz+1)/umu(0,it:it+1,ix:ix+1,iy:iy+1,iz:iz+1)) 
        if (abs(v_intpl(mu))> 1.d0) then
         print *,"Warning! In cell ",it,ix,iy,iz, " interpolated v is: ",v_intpl(1:3)
        endif
       end do
       u_intpl(1:3) = v_intpl(1:3); u_intpl(0) = 1.d0
       u_intpl = u_intpl/sqrt(LorentzProduct(u_intpl,u_intpl))

       !Calculate Cooper-Frye E, B, S assuming thermal equilibrium
       umu_sigmu = LorentzProduct(u_intpl(0:3), dSigma(0:3,i))
       E_CF(it+1) = E_CF(it+1) + E_hadgas(T_intpl, muB_intpl, muS_intpl)*u_intpl(0)*umu_sigmu+&
                                 pressure_hadgas(T_intpl, muB_intpl, muS_intpl)*(u_intpl(0)*umu_sigmu - dSigma(0,i))
       B_CF(it+1) = B_CF(it+1) + B_hadgas(T_intpl, muB_intpl, muS_intpl)*umu_sigmu
       S_CF(it+1) = S_CF(it+1) + S_hadgas(T_intpl, muB_intpl, muS_intpl)*umu_sigmu

       !Fill histograms with values on surface
       do mu=0,3
        call HistAdd(hyd_surf_sig(mu), dSigma(mu,i), 1.d0)      
        call HistAdd(hyd_surf_Vmid(mu), Vmid(mu,i),  1.d0)      
       end do 
       do mu=1,3
        call HistAdd(hyd_surf_v(mu), v_intpl(mu), 1.d0)
       end do
       call HistAdd(hyd_surf_gam, u_intpl(0), 1.d0)
       call HistAdd(hyd_surf_e,      e_intpl, 1.d0)
       call HistAdd(hyd_surf_nb,    nb_intpl, 1.d0)
       call HistAdd(hyd_surf_T,      T_intpl, 1.d0)
       call HistAdd(hyd_surf_mub,  muB_intpl, 1.d0)
       call HistAdd(hyd_surf_mus,  muS_intpl, 1.d0)
       call HistAdd(pres_CF, pressure_hadgas(T_intpl, muB_intpl, muS_intpl), 1.d0)
       call HistAdd(pres_x, presx_intpl, 1.d0)
       call HistAdd(pres_y, presy_intpl, 1.d0)
       call HistAdd(pres_z, presz_intpl, 1.d0)

       !Fill spectra histos
       if (calculate_hyd_spectra) then
       do j=1,Nphyd

         !dN+/dy
         do k = - hyd_spect_y(j)%Nbin, hyd_spect_y(j)%Nbin
          bin_val = BinToValue(hyd_spect_y(j),k)
          hyd_spect_y(j)%h(k) = hyd_spect_y(j)%h(k) + &
           dNdy(bin_val, phyd(j)%m, T_intpl, muB_intpl * phyd(j)%B + muS_intpl * phyd(j)%S,&
                u_intpl(0:3), dSigma(0:3,i), .TRUE., 2.d0*abs(phyd(j)%B) - 1.d0 , phyd(j)%g*1.d0) 
         end do

         !dN+/dy - d|N-|/dy
         do k = - hyd_neg_y(j)%Nbin, hyd_neg_y(j)%Nbin
          bin_val = BinToValue(hyd_neg_y(j),k)
          hyd_neg_y(j)%h(k) = hyd_neg_y(j)%h(k) + &
           dNdy(bin_val, phyd(j)%m, T_intpl, muB_intpl * phyd(j)%B + muS_intpl * phyd(j)%S,&
                u_intpl(0:3), dSigma(0:3,i), .FALSE., 2.d0*abs(phyd(j)%B) - 1.d0 , phyd(j)%g*1.d0) 
         end do

         !dN+/pt dpt
         do k = - hyd_spect_pt(j)%Nbin, hyd_spect_pt(j)%Nbin
          bin_val = BinToValue(hyd_spect_pt(j),k)
          hyd_spect_pt(j)%h(k) = hyd_spect_pt(j)%h(k) + &
           dNptdpt(bin_val, phyd(j)%m, T_intpl, muB_intpl * phyd(j)%B + muS_intpl * phyd(j)%S,&
                u_intpl(0:3), dSigma(0:3,i), .TRUE., 2.d0*abs(phyd(j)%B) - 1.d0 , phyd(j)%g*1.d0) 
         end do

         !dN+/pt dpt - d|N-|/pt dpt
         do k = - hyd_neg_pt(j)%Nbin, hyd_neg_pt(j)%Nbin
          bin_val = BinToValue(hyd_neg_pt(j),k)
          hyd_neg_pt(j)%h(k) = hyd_neg_pt(j)%h(k) + &
           dNptdpt(bin_val, phyd(j)%m, T_intpl, muB_intpl * phyd(j)%B + muS_intpl * phyd(j)%S,&
                u_intpl(0:3), dSigma(0:3,i), .FALSE., 2.d0*abs(phyd(j)%B) - 1.d0 , phyd(j)%g*1.d0) 
         end do
         

       end do
       endif

      end do

   endif


  end do; end do; end do
 end do

 close(8)
 close(9)

!======================================================================================================
!  Get particle crossings of hypersurface by explicit counting
!======================================================================================================
 print *
 print *,"BY PARTICLES"
 call pprt_init()

 open(unit=8,  file="output/"//trim(out_dir)//"/prt/surfTX.txt")
 open(unit=9,  file="output/"//trim(out_dir)//"/prt/surfTZ.txt")
 open(unit=10, file="output/"//trim(out_dir)//"/prt/neg_cr.txt")

 E_in = 0.d0; B_in = 0.d0; S_in = 0.d0

 fnum = 0
 nevt = 0

 do !Input files cycle 

   Npart_file = 0
   call make_in_name(input_prename,E_collision,b_collision,fnum,input_fname)
   inquire(file=input_fname, exist = op)
   if (.not. op) then; exit; endif
   
   print *,"Reading file: ",input_fname
   open(unit=11, form="unformatted", access="sequential", file = input_fname)
   do !event cycle
     read(11, iostat = io)Npart_evt
     if (io .ne. 0) then; exit; endif
     nevt = nevt + 1
     Npart_file = Npart_file + Npart_evt

     do i=1,Npart_evt
       read(11)part
       
       !Get E, B and S inside the hypersurface by particles for each time step
       do it = tstart_step, tend_step
        if (part%ri(0) > it*dt .or. part%rf(0) < it*dt) then; cycle; endif
        call GetTrajPoint(it*dt, r1, part%ri, part%rf)
        ar = GetEdens(r1)


        if (ar .ge. e0) then
          E_in(it) = E_in(it) + part%p(0)
          B_in(it) = B_in(it) + part%B
          S_in(it) = S_in(it) + part%S
        endif
       end do 

       !Get particle crossings of the hyper
       if (calculate_prt_spectra) then
       dt_ptraj = 0.1d0*dt
       t_ptraj = part%ri(0)

       do while (t_ptraj<=part%rf(0))
         
         call GetTrajPoint(t_ptraj,                           r1, part%ri, part%rf)
         call GetTrajPoint(min(t_ptraj+dt_ptraj, part%rf(0)), r2, part%ri, part%rf)
         ar = GetEdens(r1)
         br = GetEdens(r2)       
         if (ar > 0.d0) then

            if (ar > e0 .and. br < e0) then !Crossing from inside to outside
             if (abs(r1(1)/dx) < 0.3d0 .and. abs(r1(2)/dy) < 0.3d0) then; write(9,*)r1(0),r1(3); endif
             if (abs(r1(3)/dz) < 0.3d0 .and. abs(r1(2)/dy) < 0.3d0) then; write(8,*)r1(0),r1(1); endif
             call prt_addto_spectr_hist(part, .FALSE.)
            endif

            if (ar < e0 .and. br > e0) then !Crossing from outside to inside
             write(10,*)r1(0:3)
             call prt_addto_spectr_hist(part, .TRUE.)

            endif
         endif
         t_ptraj = t_ptraj + dt_ptraj
       end do
       endif

     end do !particles in event cycle
   end do ! end event cycle
   close(11)
   print *,Npart_file, " lines read from file ", input_fname
   !print *,"ri(Npart), rf(Npart): ", part(395)%ri(0:3), part(395)%rf(0:3)
   fnum = fnum + 1  

 end do !End input files cycle
 print *,"Total ", nevt, " events read from ", fnum, " files."

 E_in = E_in/nevt
 B_in = B_in/nevt
 S_in = S_in/nevt

 do i=1, Npprt
  prt_spect_y(i)%h  = prt_spect_y(i)%h/nevt
  prt_spect_pt(i)%h = prt_spect_pt(i)%h/nevt
  prt_neg_y(i)%h    = prt_neg_y(i)%h/nevt
  prt_neg_pt(i)%h   = prt_neg_pt(i)%h/nevt
 end do
 
 close(8)
 close(9)
 close(10)

!======================================================================================================
!  Finalize and output
!======================================================================================================

 open(unit=8, file = "output/"//trim(out_dir)//"/EBS_cons.txt")
  do it = tstart_step, tend_step
   write(8,*)it*dt, E_fl(it+1), E_in(it) - E_in(it+1), B_fl(it+1), B_in(it) - B_in(it+1), S_fl(it+1), S_in(it) - S_in(it+1)
  end do
 close(8)

 open(unit=8, file = "output/"//trim(out_dir)//"/CF_flows.txt")
  do it = tstart_step, tend_step
   write(8,*)it*dt, E_CF(it+1), B_CF(it+1), S_CF(it+1)
  end do
 close(8)


 call allhist_finalize_and_output(out_dir)


 deallocate(Tmn)
 deallocate(jB)
 deallocate(jS)
 deallocate(umu)
 deallocate(EdensL)
 deallocate(TL11)
 deallocate(TL22)
 deallocate(TL33)
 deallocate(mom_tot)
 deallocate(B_tot)
 deallocate(S_tot)
 deallocate(I3_tot)
 deallocate(E_fl)
 deallocate(E_in)
 deallocate(B_fl)
 deallocate(B_in)
 deallocate(S_fl)
 deallocate(S_in)

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
 call system('mkdir -p ./output/'//run_folder)
 call system('mkdir -p ./output/'//trim(run_folder)//"/hyd")
 call system('mkdir -p ./output/'//trim(run_folder)//"/prt")
 call system('mkdir -p ./output/'//trim(run_folder)//"/hyd/spectra")
 call system('mkdir -p ./output/'//trim(run_folder)//"/hyd/spectra/y")
 call system('mkdir -p ./output/'//trim(run_folder)//"/hyd/spectra/pt")
 call system('mkdir -p ./output/'//trim(run_folder)//"/hyd/neg_contr")
 call system('mkdir -p ./output/'//trim(run_folder)//"/hyd/neg_contr/y")
 call system('mkdir -p ./output/'//trim(run_folder)//"/hyd/neg_contr/pt")
 call system('mkdir -p ./output/'//trim(run_folder)//"/hyd/distr_on_hyper")
 call system('mkdir -p ./output/'//trim(run_folder)//"/prt/spectra")
 call system('mkdir -p ./output/'//trim(run_folder)//"/prt/spectra/y")
 call system('mkdir -p ./output/'//trim(run_folder)//"/prt/spectra/pt")
 call system('mkdir -p ./output/'//trim(run_folder)//"/prt/neg_contr")
 call system('mkdir -p ./output/'//trim(run_folder)//"/prt/neg_contr/y")
 call system('mkdir -p ./output/'//trim(run_folder)//"/prt/neg_contr/pt")
 
end subroutine


!===========================================================================================================
subroutine allhist_init()
implicit none
 !rapidity
 double precision, parameter :: yl=-3.d0, yu=3.d0
 integer, parameter :: Nbiny = 30
 !pt - transverse momentum
 double precision, parameter :: ptl=0.d0, ptu=3.d0
 integer, parameter :: Nbinpt = 30

 integer i

 do i=1, max_plot_part
  call HistNew(hyd_spect_y (i),   Nbiny,  yl,  yu)
  call HistNew(hyd_spect_pt(i),  Nbinpt, ptl, ptu)
  call HistNew(  hyd_neg_y (i),   Nbiny,  yl,  yu)
  call HistNew(  hyd_neg_pt(i),  Nbinpt, ptl, ptu)
  call HistNew(prt_spect_y (i),   Nbiny,  yl,  yu)
  call HistNew(prt_spect_pt(i),  Nbinpt, ptl, ptu)
  call HistNew(  prt_neg_y (i),   Nbiny,  yl,  yu)
  call HistNew(  prt_neg_pt(i),  Nbinpt, ptl, ptu)
 end do

 do i=0,3;  call HistNew(hyd_surf_sig(i),  20, -1.5d0, 1.5d0); end do
 call HistNew(hyd_surf_Vmid(0), 20, -0.2d0*dt, 1.2d0*dt);
 call HistNew(hyd_surf_Vmid(1), 20, -0.2d0*dx, 1.2d0*dx);  
 call HistNew(hyd_surf_Vmid(2), 20, -0.2d0*dy, 1.2d0*dy);
 call HistNew(hyd_surf_Vmid(3), 20, -0.2d0*dz, 1.2d0*dz);  
 do i=1,3;  call HistNew(hyd_surf_v(i),    20, -1.d0,  1.d0); end do
 call HistNew(hyd_surf_gam, 40, 0.5d0, 4.d0)
 call HistNew(hyd_surf_e,   50,  0.d0, 0.4d0)
 call HistNew(hyd_surf_nb,  50,  0.d0, 0.3d0)
 call HistNew(hyd_surf_T,   50,  0.d0, 0.25d0)
 call HistNew(hyd_surf_mub, 50 , 0.d0, 0.8d0)
 call HistNew(hyd_surf_mus, 50 , 0.d0, 0.2d0)

 call HistNew(pres_CF, 50 , 0.d0, 0.4d0)
 call HistNew(pres_x , 50 , 0.d0, 0.4d0)
 call HistNew(pres_y , 50 , 0.d0, 0.4d0)
 call HistNew(pres_z , 50 , 0.d0, 0.4d0)

end subroutine

!===========================================================================================================
subroutine prt_addto_spectr_hist(par_to_add, add_to_neg)
!Adds particle to the corresponding histograms
!For instance if par_to_add is K- then it must be added to corresponding y and pt histograms
!add_to_neg = .FALSE. (do not add to negative contribution histos) or .TRUE. (add to neg. contr.)
implicit none
  type(gen_part), intent(in) :: par_to_add
  logical, intent(in) :: add_to_neg
  integer num, i, j 
  logical i3match
  
  num = -1

  do i=1, Npprt
    if (par_to_add%id .ne. pprt(i)%ityp) then; cycle; endif

    i3match = .FALSE.
    do j=1, pprt(i)%Ni3
     if (pprt(i)%i3(j) == par_to_add%I3) then; i3match = .TRUE.; endif
    end do
    if (i3match) then; num = i; endif
  end do

  if (num > 0) then
  !Add particle to histograms with index num
   if (add_to_neg) then
    call HistAdd(prt_neg_y(num),   GetY(par_to_add%p(0:3)), 1.d0)
    call HistAdd(prt_neg_pt(num), GetPt(par_to_add%p(0:3)), 1.d0)
   else
    call HistAdd(prt_spect_y(num),   GetY(par_to_add%p(0:3)), 1.d0)
    call HistAdd(prt_spect_pt(num), GetPt(par_to_add%p(0:3)), 1.d0)
   endif
  endif

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
 integer i,j
 character(len=*)dname

 !Histograms of values on the the surface
 !NO: 1) From delta N(delta y) histograms make dN/dy, where N - number of hypersurface pieces
 ! 2) Write histograms to files

 print *,"Writing histograms of values on the hyper surface to directory: output/",trim(dname)
 !do i=0,3;  hyd_surf_sig(i)%h = hyd_surf_sig(i)%h/HistBinSize(hyd_surf_sig(i)); end do
 !do i=0,3;  hyd_surf_Vmid(i)%h = hyd_surf_Vmid(i)%h/HistBinSize(hyd_surf_Vmid(i)); end do
 !do i=1,3;  hyd_surf_v(i)%h = hyd_surf_v(i)%h/HistBinSize(hyd_surf_v(i)); end do
 !hyd_surf_gam%h = hyd_surf_gam%h/HistBinSize(hyd_surf_gam)
 !hyd_surf_e%h = hyd_surf_e%h/HistBinSize(hyd_surf_e)
 !hyd_surf_nb%h = hyd_surf_nb%h/HistBinSize(hyd_surf_nb)
 !hyd_surf_T%h = hyd_surf_T%h/HistBinSize(hyd_surf_T)
 !hyd_surf_mub%h = hyd_surf_mub%h/HistBinSize(hyd_surf_mub)
 !hyd_surf_mus%h = hyd_surf_mus%h/HistBinSize(hyd_surf_mus)

 call HistPrint(hyd_surf_sig(0),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/sig0.txt")
 call HistPrint(hyd_surf_sig(1),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/sig1.txt")
 call HistPrint(hyd_surf_sig(2),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/sig2.txt")
 call HistPrint(hyd_surf_sig(3),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/sig3.txt")

 call HistPrint(hyd_surf_Vmid(0),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/Vmid0.txt")
 call HistPrint(hyd_surf_Vmid(1),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/Vmid1.txt")
 call HistPrint(hyd_surf_Vmid(2),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/Vmid2.txt")
 call HistPrint(hyd_surf_Vmid(3),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/Vmid3.txt")

 call HistPrint(hyd_surf_v(1),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/v1.txt")
 call HistPrint(hyd_surf_v(2),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/v2.txt")
 call HistPrint(hyd_surf_v(3),"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/v3.txt")
 call HistPrint(hyd_surf_gam ,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/gam.txt")

 call HistPrint(hyd_surf_e,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/e.txt")
 call HistPrint(hyd_surf_nb,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/nb.txt")
 call HistPrint(hyd_surf_T,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/T.txt")
 call HistPrint(hyd_surf_muB,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/muB.txt")
 call HistPrint(hyd_surf_muS,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/muS.txt")

 call HistPrint(pres_CF,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/p_CF.txt")
 call HistPrint(pres_x ,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/p_x.txt")
 call HistPrint(pres_y ,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/p_y.txt")
 call HistPrint(pres_z ,"output/"//trim(adjustl(dname))//"/hyd/distr_on_hyper/p_z.txt")

 !Histograms for spectra by particles

 if (calculate_prt_spectra) then
 
 do i=1, Npprt
  prt_neg_y(i)%h = prt_neg_y(i)%h / prt_spect_y(i)%h        !make (dN-/dy)/(dN+/dy) ratio
  prt_neg_pt(i)%h = prt_neg_pt(i)%h / prt_spect_pt(i)%h     !make (dN-/dpt)/(dN+/dpt) ratio
  prt_spect_y(i)%h = prt_spect_y(i)%h / HistBinSize(prt_spect_y(i))  !make dN+/dy from dN+/dy * delta y
  prt_spect_pt(i)%h = prt_spect_pt(i)%h / HistBinSize(prt_spect_pt(i))  !make dN+/dpt from dN+/dpt * delta pt
  !make dN+/(pt*dpt)
  do j=-prt_spect_pt(i)%Nbin, prt_spect_pt(i)%Nbin
   prt_spect_pt(i)%h(j) = prt_spect_pt(i)%h(j) / BinToValue(prt_spect_pt(i),j)
  end do
  
  call HistPrint(prt_spect_y(i),"output/"//trim(adjustl(dname))//"/prt/spectra/y/"//trim(pprt(i)%pname)//".txt")
  call HistPrint(prt_spect_pt(i),"output/"//trim(adjustl(dname))//"/prt/spectra/pt/"//trim(pprt(i)%pname)//".txt")
  call HistPrint(prt_neg_y(i),"output/"//trim(adjustl(dname))//"/prt/neg_contr/y/"//trim(pprt(i)%pname)//".txt")
  call HistPrint(prt_neg_pt(i),"output/"//trim(adjustl(dname))//"/prt/neg_contr/pt/"//trim(pprt(i)%pname)//".txt")
 end do

 endif


 ! Histograms for spectra bu Cooper-Frye


 if (calculate_prt_spectra) then
 do i=1, Nphyd
  !Was (neg: dN+/dy - dN-/dy,  spect: dN+/dy)
  !Make (neg: (dN-/dy)/(dN+/dy), spect: dN+/dy)
  hyd_neg_y(i)%h = (hyd_spect_y(i)%h - hyd_neg_y(i)%h)/ hyd_spect_y(i)%h
  hyd_neg_pt(i)%h = (hyd_spect_pt(i)%h - hyd_neg_pt(i)%h)/ hyd_spect_pt(i)%h  

  call HistPrint(hyd_spect_y(i),"output/"//trim(adjustl(dname))//"/hyd/spectra/y/"//trim(phyd(i)%pname)//".txt")
  call HistPrint(hyd_spect_pt(i),"output/"//trim(adjustl(dname))//"/hyd/spectra/pt/"//trim(phyd(i)%pname)//".txt")
  call HistPrint(hyd_neg_y(i),"output/"//trim(adjustl(dname))//"/hyd/neg_contr/y/"//trim(phyd(i)%pname)//".txt")
  call HistPrint(hyd_neg_pt(i),"output/"//trim(adjustl(dname))//"/hyd/neg_contr/pt/"//trim(phyd(i)%pname)//".txt") 

 end do
 endif


end subroutine



!===========================================================================================================
subroutine phyd_init()
implicit none
  integer i

  Nphyd = 6
  phyd(1)%pname = "pi";  phyd(1)%m = 0.138d0; phyd(1)%g =  3; phyd(1)%B = 0; phyd(1)%S =  0
  phyd(2)%pname = "N";   phyd(2)%m = 0.938d0; phyd(2)%g =  4; phyd(2)%B = 1; phyd(2)%S =  0
  phyd(3)%pname = "De";  phyd(3)%m = 1.232d0; phyd(3)%g = 16; phyd(3)%B = 1; phyd(3)%S =  0
  phyd(4)%pname = "La";  phyd(4)%m = 1.115d0; phyd(4)%g =  1; phyd(4)%B = 1; phyd(4)%S = -1
  phyd(5)%pname = "Kmi"; phyd(5)%m = 0.495d0; phyd(5)%g =  1; phyd(5)%B = 0; phyd(5)%S = -1
  phyd(6)%pname = "Kpl"; phyd(6)%m = 0.495d0; phyd(6)%g =  1; phyd(6)%B = 0; phyd(6)%S =  1

  print *,"=============================================================="
  print *,"Initializing particles for hydro-style calculations:          "
  print *,"  For them Cooper-Frye y and pt spectra will be calculated.   "
  print *,"=============================================================="

  write(*,'(A4,2X,A5,2X,3A3)')"name"," mass","  g", "  B", "  S"
  do i=1, Nphyd
   write(*,'(A4,2X,F5.3,2X,3I3)')phyd(i)%pname, phyd(i)%m, phyd(i)%g, phyd(i)%B, phyd(i)%S
  end do

end subroutine

!===========================================================================================================
subroutine pprt_init()
implicit none
  integer i

  Npprt = 6
  pprt(1)%pname = "pi";  pprt(1)%ityp = 101; pprt(1)%Ni3 = 3; pprt(1)%i3(1:3) = (/-2,0,2/)
  pprt(2)%pname = "N";   pprt(2)%ityp =   1; pprt(2)%Ni3 = 2; pprt(2)%i3(1:2) = (/-1,1/)
  pprt(3)%pname = "De";  pprt(3)%ityp =  17; pprt(3)%Ni3 = 4; pprt(3)%i3(1:4) = (/-3,-1,1,3/)
  pprt(4)%pname = "La";  pprt(4)%ityp =  27; pprt(4)%Ni3 = 1; pprt(4)%i3(1:1) = (/0/)
  pprt(5)%pname = "Kmi"; pprt(5)%ityp =-106; pprt(5)%Ni3 = 1; pprt(5)%i3(1:1) = (/-1/)
  pprt(6)%pname = "Kpl"; pprt(6)%ityp = 106; pprt(6)%Ni3 = 1; pprt(6)%i3(1:1) = (/1/)

  print *,"====================================================================="
  print *,"Initializing particles for by-particles calculations:                "
  print *,"  For them y and pt spectra at crossing hypersurface will be counted."
  print *,"====================================================================="

  write(*,'(A4,2X,A5,2X,A)')"name"," ityp"," 2*i3 list"
  do i=1, Npprt
   write(*,'(A4,2X,I5,2X,6I3)')pprt(i)%pname, pprt(i)%ityp, pprt(i)%i3( 1:pprt(i)%Ni3 )
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

!===============================================================================================
double precision function GetEdens(r0) result (res)
!Interpolates energy density to space-time point r0
!Uses global array EdensL - energy density on a grid
implicit none
double precision r0(0:3), v0(0:3)
integer it_r0, ix_r0, iy_r0, iz_r0

   it_r0=floor(r0(0)/dt)
   ix_r0=floor(r0(1)/dx)
   iy_r0=floor(r0(2)/dy)
   iz_r0=floor(r0(3)/dz)


   if (ix_r0 < -nx .or. ix_r0 >= nx .or. &
       iy_r0 < -ny .or. iy_r0 >= ny .or. &
       iz_r0 < -nz .or. iz_r0 >= nz .or. & 
       it_r0 >= tend_step .or. it_r0 < tstart_step) then
     res = -1.d0
     return
   endif

    v0(0)=r0(0)/dt-it_r0*1.d0
    v0(1)=r0(1)/dx-ix_r0*1.d0
    v0(2)=r0(2)/dy-iy_r0*1.d0
    v0(3)=r0(3)/dz-iz_r0*1.d0

    res = Cube4Intpl(v0, EdensL(it_r0:it_r0+1, ix_r0:ix_r0+1, iy_r0:iy_r0+1, iz_r0:iz_r0+1))

end function GetEdens

end program analyze_part_lines
