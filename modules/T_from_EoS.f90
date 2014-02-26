module T_from_EoS
!Purpose: 
!  1) Read equation of state from UrQMD files and get T(e,n), mu_b(e,n), mu_s(e,n)
!  2) Given T, mu_b and mu_s calculate energy density e and baryon density n
!
! T(e,n), mu_b(e,n), mu_s(e,n) are copied from 1fluid.f code and adapted to f90 style
!
!Author: Dmytro Oliinychenko

      !The tables
      double precision, parameter :: pi = 4d0*datan(1d0)
      double precision, parameter :: e_gr = 146.51751415742d-3   !Gev/fm^3 - ground state energy
      double precision, parameter :: n_gr = 0.15891d0            !fm^-3 - ground state baryon density

      double precision, private :: ttab(0:2000,0:400), mutab(0:2000,0:400), mustab(0:2000,0:400)
      double precision, private :: ttab2(0:2000,0:400), mutab2(0:2000,0:400), mustab2(0:2000,0:400)
      double precision, private :: ttab3(0:2000,0:400), mutab3(0:2000,0:400), mustab3(0:2000,0:400)  
      logical, private :: tables_read = .FALSE.

     !particle table
     integer, parameter, private :: Max_part=1000
     integer, private ::  Npart !total number of particles dealt with

     !particle charactristics
     double precision, private :: g(Max_part),mass(Max_part),bar(Max_part),&
                            str(Max_part), I3(Max_part)
     !particle names
     character*16, private :: Part_name(Max_part)

contains
!------------------------------------------------------------------------
!     
!     subroutine to test this module
!     
      subroutine T_from_EoS_test()
      implicit none
      double precision T,mu,mus,e,nb

      
    
      e=0.15d0
      nb = 0.06d0
      write(*,'(A,f10.4,A)')"Energy density is ",e," GeV/fm^3"  
      write(*,'(A,f10.4,A)')"Baryon density is ",nb," fm^-3"


      write(*,'(A)')"Getting T and mu for these energy and baryon density:"
      T = GetTemp(e,nb)
      mu = GetMuB(e,nb)
      mus = GetMuS(e,nb)
      write(*,'(A,f10.4,A)')"Temperature is ",T, " MeV"
      write(*,'(A,f10.4,A)')"Baryon chemical potential is ",mu, " MeV"
      write(*,'(A,f10.4,A)')"Strange chemical potential is ",mus, " MeV"

      write(*,'(/,A)')"Now from T and mu obtain energy density and baryon density."
      write(*,'(A,f10.4,A)')"Energy density from my hadgas:",E_hadgas(T*1d-3, mu*1d-3, mus*1d-3), " GeV/fm^3"
      write(*,'(A,f10.4,A)')"Baryon density from my hadgas:",B_hadgas(T*1d-3, mu*1d-3,mus*1d-3), " fm^-3"
      write(*,'(A,f10.4,A)')"Net strangeness from my hadgas:",S_hadgas(T*1d-3, mu*1d-3, mus*1d-3), " fm^-3"
            
      end subroutine

!------------------------------------------------------------------------
!     
!     subroutine which reads in EoS matrices.
!     
!     readeos1 reads in pure hadronic EoS
!
      subroutine read_eos_tables()
      implicit none
      integer in, ie, i, seq_num
      double precision  t,mu,e,p,n,s,mus,lam
  
 
      open(unit=53, file="/home/oliiny/Work/Tmn_project/hyd_vs_prt_v2/eosfiles/hadgas_eos.dat")
      open(unit=54, file="/home/oliiny/Work/Tmn_project/hyd_vs_prt_v2/eosfiles/hg_eos_small.dat")     
      open(unit=55, file="/home/oliiny/Work/Tmn_project/hyd_vs_prt_v2/eosfiles/hg_eos_mini.dat") 

      do in = 0,400
         do ie = 0,2000
            read(53,'(2(1x,f8.3),6(1x,e15.7))') t,mu,e,p,n,s,mus,lam
            ttab(ie,in) = t
            mutab(ie,in) = mu
            mustab(ie,in) = mus
         end do
      end do
      close(53)

      do in = 0,200
         do ie = 0,200
            read(54,'(2(1x,f8.3),6(1x,e15.7))') t,mu,e,p,n,s,mus,lam
            ttab2(ie,in) = t
            mutab2(ie,in) = mu
            mustab2(ie,in) = mus
          end do
      end do
      close(54)      

      do in = 0,200
         do ie = 0,200
            read(55,'(2(1x,f8.3),6(1x,e15.7))') t,mu,e,p,n,s,mus,lam
            ttab3(ie,in) = t
            mutab3(ie,in) = mu
            mustab3(ie,in) = mus
         end do
      end do
      close(55)


      open(unit=8, file = "particles.txt")

      read (8,*) Npart !reading number of particles from file
      do i = 1, Npart ! reading all the particle table
        read(8,'(I3,2X,A16,F7.0,F11.5,F6.1,2F6.0)')seq_num,&
                      Part_name(i),g(i),mass(i),I3(i),bar(i),str(i)
        !print *,Part_name(i),g(i),mass(i),I3(i),bar(i),str(i)
      end do
    
      close(8)
      print *,"Particle table consumed, total particles:", Npart;


      tables_read = .TRUE.

 
      return
      end subroutine

!=====================GetTemp=================================================

double precision function GetTemp(e_in,n_in) result (temp)
      implicit none     
!------------------------------------------------------------------------------
!     This function determines the temperature of the
!     Hadron Gas EoS.
   
      double precision, intent(in) :: e_in, n_in
      double precision e,n
      double precision de,dn,p1,p2,p3,p4,p13,p24,p12,p34

         if (.not. tables_read) then
           print *,"Error: First read in the table then ask for value from it."
           call exit()
         endif
         
         e = e_in / e_gr    ! from GeV/fm^3 to table values
         n = n_in / n_gr     ! from GeV/fm^3 to table values

         if (e.le.400.0d0 .and. n.lt.40.d0 .and. n.gt.0.d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = ttab3(idint(e/de),idint(n/dn))
               p2 = ttab3(idint(e/de)+1,idint(n/dn))
               p3 = ttab3(idint(e/de),idint(n/dn)+1)
               p4 = ttab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.&
                   ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = ttab2(idint(e/de),idint(n/dn))
               p2 = ttab2(idint(e/de)+1,idint(n/dn))
               p3 = ttab2(idint(e/de),idint(n/dn)+1)
               p4 = ttab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then
               de = 0.2d0
               dn = 0.1d0

               p1 = ttab(idint(e/de),idint(n/dn))
               p2 = ttab(idint(e/de)+1,idint(n/dn))
               p3 = ttab(idint(e/de),idint(n/dn)+1)
               p4 = ttab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if 
            
         else 
            temp = 0d0         
         end if

      return
      end function GetTemp       
              
!
!==================GetMuB==============================================
!

double precision function GetMuB(e_in,n_in) result (chem)
!
!  This function-subprogram determines the baryon chemical potential of the
!  underlying EoS.
!
      implicit none
      double precision, intent(in) :: e_in, n_in
      double precision e,n
      double precision de,dn,p1,p2,p3,p4,p13,p24,p12,p34


         if (.not. tables_read) then
           print *,"Error: First read in the table then ask for value from it."
           call exit()
         endif

         e = e_in /  e_gr   ! from GeV/fm^3 to table values
         n = n_in / n_gr     ! from GeV/fm^3 to table values


         if (e.lt.400.0d0 .and. n.lt.40.d0 .and. n.gt.0.d0) then


            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mutab3(idint(e/de),idint(n/dn))
               p2 = mutab3(idint(e/de)+1,idint(n/dn))
               p3 = mutab3(idint(e/de),idint(n/dn)+1)
               p4 = mutab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            endif
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and. &
                ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mutab2(idint(e/de),idint(n/dn))
               p2 = mutab2(idint(e/de)+1,idint(n/dn))
               p3 = mutab2(idint(e/de),idint(n/dn)+1)
               p4 = mutab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            endif
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then       
               de = 0.2d0
               dn = 0.1d0
               p1 = mutab(idint(e/de),idint(n/dn))
               p2 = mutab(idint(e/de)+1,idint(n/dn))
               p3 = mutab(idint(e/de),idint(n/dn)+1)
               p4 = mutab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
               
            endif
         else 
            chem = 0d0
         endif

      !Turn quark potential to baryon one
      chem = 3d0 *chem 
         
      return
      end function GetMuB 




!
!==================GetMuS==============================================
!
double precision function GetMuS(e_in,n_in) result (chem)
   implicit none
   double precision, intent(in) :: e_in, n_in
   chem = GetMuB(e_in,n_in)/3d0 - GetMuS_tab(e_in,n_in) 
end function GetMuS

double precision function GetMuS_tab(e_in,n_in) result (chem)
!
!  This function-subprogram determines the baryon chemical potential of the
!  underlying EoS.
!
      implicit none
      double precision, intent(in) :: e_in, n_in
      double precision e,n
      double precision de,dn,p1,p2,p3,p4,p13,p24,p12,p34


         if (.not. tables_read) then
           print *,"Error: First read in the table then ask for value from it."
           call exit()
         endif

         e = e_in / e_gr    ! from GeV/fm^3 to table values
         n = n_in / n_gr     ! from GeV/fm^3 to table values


         if (e.lt.400.0d0 .and. n.lt.40.d0 .and. n.gt.0.d0) then


            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mustab3(idint(e/de),idint(n/dn))
               p2 = mustab3(idint(e/de)+1,idint(n/dn))
               p3 = mustab3(idint(e/de),idint(n/dn)+1)
               p4 = mustab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            endif
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and. &
                ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mustab2(idint(e/de),idint(n/dn))
               p2 = mustab2(idint(e/de)+1,idint(n/dn))
               p3 = mustab2(idint(e/de),idint(n/dn)+1)
               p4 = mustab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            endif
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then       
               de = 0.2d0
               dn = 0.1d0
               p1 = mustab(idint(e/de),idint(n/dn))
               p2 = mustab(idint(e/de)+1,idint(n/dn))
               p3 = mustab(idint(e/de),idint(n/dn)+1)
               p4 = mustab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
               
            endif
         else 
            chem = 0d0
         endif

        
      return
      end function GetMuS_tab



!=======================================================================
double precision function E_hadgas(T, mu, mus) result (res)
implicit none
 ! T - temperature [GeV]
 ! u - 4-velocity with lower index
 ! dSigma - normal 4-vector with lower index [fm^3]
 ! mu - baryon chemical potential
 ! mus - strange chemical potential
 ! Calculates energy flow through dSigma from all particles in a table using Cooper-Frye formula
 double precision, intent(in) :: T, mu, mus
 double precision bose_fermi_boltz
 integer i

 if (.not. tables_read) then
   print *,"Error: First read in the table then ask for value from it."
   call exit()
 endif
 
 res =0d0
 do i=1, Npart
  if (abs(abs(bar(i)) - 1d0) < 1d-9) then
   bose_fermi_boltz = 1d0
  else if (abs(bar(i)) < 1d-9) then
   bose_fermi_boltz = -1d0
  else
   print *,"Particle table trouble: bar = ", bar(i)
   call exit()
  endif
  
  res = res + E_part(mass(i), T, mu*bar(i) + mus*str(i), bose_fermi_boltz, g(i))
 end do


end function E_hadgas

!=======================================================================
double precision function B_hadgas(T, mu, mus) result (res)
implicit none
 ! T - temperature [GeV]
 ! u - 4-velocity with lower index
 ! dSigma - normal 4-vector with lower index [fm^3]
 ! mu - baryon chemical potential
 ! mus - strange chemical potential
 ! Calculates B flow through dSigma from all particles in a table using Cooper-Frye formula
 double precision, intent(in) :: T, mu, mus
 double precision bose_fermi_boltz
 integer i

 if (.not. tables_read) then
   print *,"Error: First read in the table then ask for value from it."
   call exit()
 endif

 
 res =0d0
 do i=1, Npart
  if (abs(abs(bar(i)) - 1d0) < 1d-9) then
   bose_fermi_boltz = 1d0
  else if (abs(bar(i)) < 1d-9) then
   bose_fermi_boltz = -1d0
  else
   print *,"Particle table trouble: bar = ", bar(i)
   call exit()
  endif
  
  res = res + bar(i)*N_part(mass(i), T, mu*bar(i) + mus*str(i), bose_fermi_boltz, g(i))
 end do


end function B_hadgas


!=======================================================================
double precision function S_hadgas(T, mu, mus) result (res)
implicit none
 ! T - temperature [GeV]
 ! u - 4-velocity with lower index
 ! dSigma - normal 4-vector with lower index [fm^3]
 ! mu - baryon chemical potential
 ! mus - strange chemical potential
 ! Calculates S flow from all particles in a table using Cooper-Frye formula
 double precision, intent(in) :: T, mu, mus
 double precision bose_fermi_boltz
 integer i

 if (.not. tables_read) then
   print *,"Error: First read in the table then ask for value from it."
   call exit()
 endif

 
 res =0d0
 do i=1, Npart
  if (abs(abs(bar(i)) - 1d0) < 1d-9) then
   bose_fermi_boltz = 1d0
  else if (abs(bar(i)) < 1d-9) then
   bose_fermi_boltz = -1d0
  else
   print *,"Particle table trouble: bar = ", bar(i)
   call exit()
  endif
  
  res = res + str(i)*N_part(mass(i), T, mu*bar(i) + mus*str(i), bose_fermi_boltz, g(i))
 end do


end function S_hadgas


!=======================================================================
double precision function pressure_hadgas(T, mu, mus) result (res)
implicit none
 ! T - temperature [GeV]
 ! u - 4-velocity with lower index
 ! dSigma - normal 4-vector with lower index [fm^3]
 ! mu - baryon chemical potential
 ! mus - strange chemical potential
 ! Calculates S flow from all particles in a table using Cooper-Frye formula
 double precision, intent(in) :: T, mu, mus
 double precision bose_fermi_boltz
 integer i

 if (.not. tables_read) then
   print *,"Error: First read in the table then ask for value from it."
   call exit()
 endif

 
 res =0d0
 do i=1, Npart
  if (abs(abs(bar(i)) - 1d0) < 1d-9) then
   bose_fermi_boltz = 1d0
  else if (abs(bar(i)) < 1d-9) then
   bose_fermi_boltz = -1d0
  else
   print *,"Particle table trouble: bar = ", bar(i)
   call exit()
  endif
  
  res = res + N_part(mass(i), T, mu*bar(i) + mus*str(i), bose_fermi_boltz, g(i))
  res = res*T
 end do


end function pressure_hadgas


!==========================================================================
double precision function E_part(m, T, mu, bose_fermi_boltz, degen)
 implicit none
 ! m - mass [GeV]
 ! T - temperature [GeV]
 ! u - 4-velocity with lower index
 ! dSigma - normal 4-vector with lower index [fm^3]
 ! Returns energy flow calculated by Cooper-Frye formula
 ! bose_fermi_boltz = 1 (Fermi distribution), 0 (Boltzmann), -1 (Bose)
 ! degen - degeneracy = (2S+1)(2I+1)
 double precision, intent(in) :: m, T, mu, bose_fermi_boltz, degen
 double precision hbarc, factor, tot_sum, BK(0:10)

 integer j, NSeries

 tot_sum = 0d0
 NSeries = 3
 
 if (abs(bose_fermi_boltz) < 1d-9) then

    call BESSK(3,m/T,BK)
    tot_sum = exp(mu/T)*(3d0*BK(3)+BK(1))

 else if (abs(bose_fermi_boltz+1d0) < 1d-9) then

    do j=1,NSeries
     call BESSK(3,j*m/T,BK)
     tot_sum = tot_sum + exp(j*mu/T)*(3d0*BK(3)+BK(1))/j
    end do

 else if (abs(bose_fermi_boltz-1d0) < 1d-9) then

    do j=1,NSeries
     call BESSK(3,j*m/T,BK)
     ! (-1)^(j-1)
     if (mod(j,2)==1) then
      tot_sum = tot_sum + exp(j*mu/T)*(3d0*BK(3)+BK(1))/j
     else
      tot_sum = tot_sum - exp(j*mu/T)*(3d0*BK(3)+BK(1))/j
     endif
    end do
 else

   print *,"Wrong option: bose_fermi_boltzmann. Exiting."
   call exit();

 endif

 hbarc = 0.19732d0 !hbar*c [GeV*fm]
 factor = pi* degen /(2.d0*pi*hbarc)**3 * T * m**3
 E_part = tot_sum*factor
 return

end function E_part

double precision function N_part(m, T, mu, bose_fermi_boltz, degen)
 implicit none
 ! m - mass [GeV]
 ! T - temperature [GeV]
 ! u - 4-velocity with lower index
 ! dSigma - normal 4-vector with lower index [fm^3]
 ! Returns energy flow calculated by Cooper-Frye formula
 ! bose_fermi_boltz = 1 (Fermi distribution), 0 (Boltzmann), -1 (Bose)
 ! degen - degeneracy = (2S+1)(2I+1)
 double precision, intent(in) :: m, T, mu, bose_fermi_boltz, degen
 double precision hbarc, factor, tot_sum, BK(0:10)

 integer j, NSeries

 tot_sum = 0d0
 NSeries = 3
 
 if (abs(bose_fermi_boltz) < 1d-9) then

    call BESSK(2,m/T,BK)
    tot_sum = exp(mu/T)*BK(2)

 else if (abs(bose_fermi_boltz+1d0) < 1d-9) then

    do j=1,NSeries
     call BESSK(2,j*m/T,BK)
     tot_sum = tot_sum + exp(j*mu/T)*BK(2)/j
    end do

 else if (abs(bose_fermi_boltz-1d0) < 1d-9) then

    do j=1,NSeries
     call BESSK(2,j*m/T,BK)
     ! (-1)^(j-1)
     if (mod(j,2)==1) then
      tot_sum = tot_sum + exp(j*mu/T)*BK(2)/j
     else
      tot_sum = tot_sum - exp(j*mu/T)*BK(2)/j
     endif
    end do

 else

   print *,"Wrong option: bose_fermi_boltzmann. Exiting."
   call exit();

 endif

 hbarc = 0.19732d0 !hbar*c [GeV*fm]
 factor = 4d0*pi* degen /(2.d0*pi*hbarc)**3 * T * m**2

 N_part = tot_sum*factor
 return

end function N_part

!=========================================================================
! ------------------------------------------------------------------------
!*             BESSEL   FUNCTION     BLOCK                 *
!*                                                         *
!*              From Numath Library By Tuan Dang Trong     *
!*              in Fortran 77 [BIBLI 18]                   *
!*                                                         *
!*                      F90 Version By J-P Moreau, Paris.  *
!*                             (www.jpmoreau.fr)           *
! ------------------------------------------------------------------------
!=========================================================================

      SUBROUTINE BESSK(N,X,BK)
      !Returns array BK containing K0(X), K1(X), ..., KN(X)
      implicit none
      integer, parameter :: maxNBES=10
      double precision, intent(in) :: X
      double precision, intent(out) :: BK(0:maxNBES)
      double precision TOX
      integer, intent (in) :: N
      integer J

      BK=0d0
      IF (N .GE. 0) THEN; BK(0) = BESSK0(X); ENDIF
      IF (N .GE. 1) THEN; BK(1) = BESSK1(X); ENDIF
      IF (X==0D0) THEN; BK = 1.D30; RETURN; ENDIF
      TOX = 2.D0/X

      DO J=1,min(N-1,maxNBES)
       BK(J+1) = BK(J-1)+DFLOAT(J)*TOX*BK(J)
      END DO
      return

      END SUBROUTINE
! ----------------------------------------------------------------------
      FUNCTION BESSK0(X)
      implicit none
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DU 3EME ESPECE D'ORDRE 0
!     POUR TOUT X REEL NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ZERO FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
      double precision X,BESSK0,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,&
                                  0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,&
                                -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      
      IF(X.EQ.0.D0) THEN
       BESSK0=1.D30
       RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
       Y=X*X/4.D0
       AX=-LOG(X/2.D0)*BESSI0(X)
       BESSK0=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
       Y=(2.D0/X)
       AX=EXP(-X)/SQRT(X)
       BESSK0=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
! ----------------------------------------------------------------------
      FUNCTION BESSK1(X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DE 3EME ESPECE D'ORDRE 1
!     POUR TOUT X REEL POSITF NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ONE FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
      implicit none
      double precision X,BESSK1,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,0.15443144D0,-0.67278579D0,-0.18156897D0,&
                                -0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,0.1504268D-1,&
                                -0.780353D-2,0.325614D-2,-0.68245D-3/
      IF(X.EQ.0.D0) THEN
       BESSK1=1.D32
       RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
       Y=X*X/4.D0
       AX=LOG(X/2.D0)*BESSI1(X)
       BESSK1=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
       Y=(2.D0/X)
       AX=EXP(-X)/SQRT(X)
       BESSK1=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
!
!     Bessel Function of the 1st kind of order zero.
!
      FUNCTION BESSI0(X)
      implicit none
      double precision X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,&
      0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X).LT.3.75D0) THEN
       Y=(X/3.75D0)**2
       BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
       AX=ABS(X)
       Y=3.75D0/AX
       BX=EXP(AX)/SQRT(AX)
       AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
       BESSI0=AX*BX
      ENDIF
      RETURN
      END
!
!     Bessel Function of the 1st kind of order one.
!
      FUNCTION BESSI1(X)
      implicit none
      double precision X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
      -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,        &
      -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(X).LT.3.75D0) THEN
       Y=(X/3.75D0)**2
       BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
       AX=ABS(X)
       Y=3.75D0/AX
       BX=EXP(AX)/SQRT(AX)
       AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
       BESSI1=AX*BX
      ENDIF
      RETURN
      END
  

end module T_from_EoS
