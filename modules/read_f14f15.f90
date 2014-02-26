module read_f14f15
! Purpose:
!      1) read 1 event from file15 of UrQMD output
!      2) read initial and final stage from file14 of UrQMD output
!      3) Create structure that contains straight line pieces of particle trajectories in event
! Author:
!      Dmytro Oliinychenko (oliiny@fias.uni-frankfurt.de)
!============================================================================================
type gen_part
 sequence

 integer :: id                          !Some id, this may be mc number or UrQMD ityp
 double precision :: ri(0:3), rf(0:3)  !space coordinates
 integer :: filled                      ! 0 - no ri, no rf; 1 - ri, no rf; 2 - both ri and rf
 double precision :: p(0:3), m         !4-momentum, mass
 integer :: B, S, I3                    !Baryon charge, strangeness, doubled isospin 3rd projection

end type gen_part

contains

subroutine read_f14_f15_event(uread_f14, uread_f15, Np_evt, part_evt)
!-----------------------------------------------------------------------------------
!
! Output: read_f14_f15_event returns structure that
!         contains all the straight line pieces of particle trajectories in event
!         in format r_initial, r_final, p, quantum numbers, enclosed in gen_part type.
!
!-----------------------------------------------------------------------------------
 implicit none
 integer :: Np_evt
 type (gen_part), dimension(:) :: part_evt
 integer uread_f14, uread_f15
 

   !Nullify everything
   call part_init(Np_evt, part_evt)

   !Read initial particles from f14
   call read_f14init(uread_f14,Np_evt, part_evt)

   !Read f15 collisions and match them to f14 initial particles and to each other
   call read_f15event(uread_f15, Np_evt, part_evt)

   !Read f14 final particles and match them to all previous
   call read_f14fin(uread_f14,Np_evt, part_evt)

   ! Test if all particles are matched and if matched properly
   call test_part(Np_evt, part_evt)

end subroutine

!===========================================================================================
subroutine read_f15_event_header(uread, E, b_imp)
implicit none

 integer, intent(in) :: uread
 double precision  b_imp, sqrtS, sig, E
 integer min_one, ev_num, mproj, mtarg, io

 read(uread,'(i8,i8,i4,i7,f8.3,3e12.4)', iostat=io)min_one, ev_num, mproj, mtarg, b_imp, sqrtS, sig, E
 if (min_one .ne. -1) then; print *,"Error! This is not an event header!"; call exit(); endif
end subroutine

!===========================================================================================
subroutine read_f15_collision_header(uread, nin, nout, proc_id, col_num, event_end)
implicit none

 integer, intent(in) :: uread
 integer nin, nout, proc_id, col_num, io
 logical event_end

 event_end = .FALSE.
 read(uread,'(i8,i8,i4,i7)', iostat=io)nin, nout, proc_id, col_num
 if ((io .ne. 0) .or. (nin==-1)) then; event_end = .TRUE.; endif
end subroutine

!===========================================================================================
subroutine part_try_add(part,id,r,p,m,B,S,I3, added)
!Subroutine tries to add particle with 
! coordinate r, !momentum p, mass m
! baryon charge B, strangeness S, doubled isospin projection I3 and id-number id
! as a final point of particle part
 implicit none
 type (gen_part) part
 double precision r(0:3), p(0:3), m
 integer B, S, I3, id
 logical added, coincide
 integer mu

 coincide = .TRUE.
 added = .FALSE.

 if (part%filled == 2) then
  !This segment is already filled.
  return

 else if (part%filled == 1) then
    !Initial point of the segment filled, need to fill final one.
    !But first check if particle matches.
    do mu=0,3
     if (dabs(part%p(mu) - p(mu)) > 1d-8) then; coincide=.FALSE.; endif  
    end do

    if ( part%id .ne. id)   then;  coincide = .FALSE.; endif 

    !if (part%B  .ne.  B)   then;  coincide = .FALSE.; endif
    !if (part%S  .ne.  S)   then;  coincide = .FALSE.; endif
    !if (part%I3 .ne. I3)   then;  coincide = .FALSE.; endif
    if (abs(part%m  -  m) > 1d-8)   then;  coincide = .FALSE.; endif

    if (coincide) then
     part%rf(0:3) = r(0:3)
     added = .TRUE.
     part%filled = 2
     if (part%B .ne. B) then
      print *,"Particle coincides, but B, part%B = ",B,part%B
     endif
     if (part%S .ne. S) then
      print *,"Particle coincides, but S, part%S = ",S,part%S
     endif
    endif

 else if (part%filled == 0) then
    !Nothing is filled, add initial point, no consistency checks needed
    part%id = id
    part%B  =  B
    part%S  =  S  
    part%I3 = I3 
    part%ri(0:3) = r(0:3)
    part%p(0:3) = p(0:3)
    part%m = m

    part%filled = 1
    added = .TRUE.
 else
    !Wrong value of part%filled. Report error.
    print *,"Error. Wrong value of part%filled: ", part%filled 
 endif

end subroutine

!===========================================================================================
subroutine read_f14init(uread14, Np, part)
!uread - number of unit to read
! Np - number of particles
! part - array of particles
implicit none

 integer, intent(in) :: uread14
 integer, intent(inout) :: Np
 type (gen_part), intent(inout), dimension(:) :: part
 integer ityp, I3, el_ch, S, B, io
 double precision r(0:3), p(0:3), m
 logical added

 integer i

 read(uread14,*,iostat = io);
 if (io .ne. 0) then; Np = 0; return; endif
 do i=1,13; read(uread14,*); end do
 read(uread14,*)Np
 read(uread14,*)

 do i=1, Np
  read(uread14,'(9e16.8,i11,2i3)') r(0:3), p(0:3), m, ityp, I3, el_ch
  B = BfromItyp(ityp)
  S = SfromItyp(ityp)
  call part_try_add(part(i),ityp,r,p,m,B,S,I3, added)
  if (.not. added) then; print *,"Error, initial particle ", i," from f14 not added."; endif
 end do

end subroutine

!===========================================================================================
subroutine read_f14fin(uread14, Np, part)
!uread - number of unit to read
! Np - number of particles
! part - array of particles
implicit none

 integer, intent(in) :: uread14
 integer, intent(inout) :: Np
 type (gen_part), intent(inout), dimension(:) :: part
 integer ityp, I3, el_ch, B, S, Nfin, n_added, io
 double precision r(0:3), p(0:3), m
 logical added

 integer i,j

 read(uread14,*,iostat = io)Nfin
 if (io.ne.0) then; return; endif
 read(uread14,*)

 do i=1, Nfin

      read(uread14,'(9e16.8,i11,2i3)', iostat = io) r(0:3), p(0:3), m, ityp, I3, el_ch
      if (io .ne. 0) then; Np=0; return; endif

      B = BfromItyp(ityp)
      S = SfromItyp(ityp)

      n_added = 0
      do j=1,Np
       call part_try_add(part(j),ityp,r,p,m,B,S,I3, added) 
       if (added) then; n_added = n_added + 1; endif
      end do

      if (n_added > 1) then
        print *,"Error: suitable particle found more than once. Collision number, final particle: ", i
      else if (n_added == 0) then
        Np = Np + 1
        call part_try_add(part(Np),ityp,r,p,m,B,S,I3, added)          
      endif

 end do

end subroutine


!===========================================================================================
subroutine part_init(Np, part)
 implicit none
 integer Np,i
 type (gen_part), dimension(:) :: part

 Np=0

 part(:)%filled = 0
 do i=0,3
  part(:)%ri(i) = 0d0
  part(:)%rf(i) = 0d0  
  part(:)%p(i)  = 0d0 
 end do
 part(:)%B  = 0
 part(:)%id = 0
 part(:)%S  = 0
 part(:)%I3 = 0

end subroutine

!===========================================================================================
subroutine read_f15event(uread, Np, part)
!uread - number of unit to read
! Np - number of particles
! part - array of particles
implicit none

 integer, intent(in) :: uread
 integer, intent(inout) :: Np
 type (gen_part), intent(inout), dimension(:) :: part
 integer B_in, B_out, S_in, S_out
 double precision E_in, E_out
 integer nin,nout,proc_id, col_num, i,j, n_added
 logical event_end, added

 integer part_index, ityp, I3, el_ch, parent_index, Ncoll, S, B
 double precision r(0:3), p(0:3), m

 event_end = .FALSE.

 
 do 

   call read_f15_collision_header(uread, nin, nout, proc_id, col_num, event_end)
   if (event_end) then; exit; endif
   B_in = 0; B_out = 0
   S_in = 0; S_out = 0
   E_in = 0d0; E_out = 0d0

   !Ignore Pauli-blocked collisions
   if (nin == 2 .AND. nout == 0) then
    do i=1, nin; read(uread,*); end do
    cycle
    !call read_f15_collision_header(uread, nin, nout, proc_id, col_num, event_end)
   endif

   !Ignore Pauli-blocked decays
   if (nin == 1 .AND. nout == 0) then
    do i=1, nin; read(uread,*); end do
    cycle 
    !call read_f15_collision_header(uread, nin, nout, proc_id, col_num, event_end)
   endif


   do i=1, nin+nout
    read(uread,'(i5,9e16.8,i11,2i3,i9,i5,i3)') part_index, r(0:3), p(0:3), m, ityp, I3, el_ch, parent_index, Ncoll, S
    B = BfromItyp(ityp)
    if (S.ne. SfromItyp(ityp)) then
     print *,"Reading f15, S ne S(ityp)"
    endif

    if (i .le. nin) then

      B_in = B_in + B
      S_in = S_in + S
      E_in = E_in + p(0)

      n_added = 0
      do j=1,Np
       call part_try_add(part(j),ityp,r,p,m,B,S,I3, added) 
       if (added) then; n_added = n_added + 1; endif
      end do

      if (n_added > 1) then
        print *,"Error: suitable particle found more than once. Collision number, particle: ", col_num, i
      else if (n_added == 0) then
        print *,"f15 ingoing particle not added. Collision number, particle: ", col_num, i
      endif
    
    else

      B_out = B_out + B
      S_out = S_out + S
      E_out = E_out + p(0)

      Np = Np + 1
      call part_try_add(part(Np),ityp,r,p,m,B,S,I3, added)

    endif

  

   end do
   
   if ((B_in .ne. B_out) .or. (S_in .ne. S_out) .or. (abs(E_out-E_in) > 1d-5)) then
    print *,"Error in collision ", col_num
    print *,"B in, out: ", B_in, B_out
    print *,"S in, out: ", S_in, S_out
    print *,"E in, out: ", E_in, E_out
    !call exit()
   endif

 end do

! print *,"Number of particles: ", Np
! print *,"Number of collisions", col_num

end subroutine 


!=========================================================================================================
subroutine test_part(Np, part)
 implicit none
 type (gen_part), intent(in), dimension(:) :: part
 integer, intent(in) :: Np
 integer N_filled(0:2), i

 N_filled = 0

 do i=1, Np


  if (part(i)%filled == 0) then
    N_filled(0) = N_filled(0) + 1
  else if (part(i)%filled == 1) then
    N_filled(1) = N_filled(1) + 1
  else if (part(i)%filled == 2) then
    N_filled(2) = N_filled(2) + 1
  else
   print *,"part(i)filled, i - ",part(i)%filled,i
  endif

 end do

 if (N_filled(0) .ne. 0) then
  print *, N_filled(0), "empty spaces, Np = ", Np, " is too large."
 endif

 if (N_filled(1) .ne. 0) then
  print *, N_filled(1), "unpaired particle points"
 endif

 if (N_filled(2) .ne. Np) then
  print *,"N_filled(2), Np: ", N_filled(2), Np
 endif

 !For initial 394 nucleons transverse momenta are non-zero, but x and y do not change.
 !This is UrQMD feature.
 do i=395,Np
  call test_part_consistency(part(i))
 end do
 
end subroutine


!=========================================================================================================
subroutine test_part_consistency(part)
 implicit none
 type (gen_part), intent(in) :: part
 double precision test(0:3), p(0:3),m
 integer mu
 logical good_part
 
 good_part = .TRUE.
 
 p(0:3) = part%p(0:3)
 m = part%m

 if (abs(part%B) > 1 .or. abs(part%S) > 3) then; good_part = .FALSE.; endif
 if (part%filled == 2) then
   do mu=0,3
    test(mu) = (part%rf(mu) - part%ri(mu))/(p(mu)/p(0))
   end do
   do mu = 1,3
    if (abs(test(0) - test(mu)) > 1d-2 .AND. &
        !to avoid getting error when p(mu) is very small
        abs(p(mu)) > 1.d-4) then
         good_part = .FALSE.
    endif
   end do
 
   if (part%rf(0) < part%ri(0)) then; good_part = .FALSE.; endif
 endif
 
 if (abs(p(0)*p(0) - p(1)*p(1) - p(2)*p(2) - p(3)*p(3) - m*m) > 1.d-3) then; good_part = .FALSE.; endif

 if (.not. good_part) then
  print *,"Warning. Particle failed tests."
  print *,"ri: ", part%ri(0:3)
  print *,"rf: ", part%rf(0:3)
  print *,"p: ", part%p(0:3)
  print *,"test: ", test(0:3)
  print *,"B: ", part%B
  print *,"S: ", part%S
  print *,"sqrt(pmu*pmu), m: ",sqrt(p(0)*p(0) - p(1)*p(1) - p(2)*p(2) - p(3)*p(3)),m
  print *,"filled = ", part%filled
  !call exit()
 endif 

end subroutine

!=========================================================================================================
integer function BfromItyp(ityp) result (Bch)
 implicit none
 integer, intent(in) :: ityp

    if (abs(ityp) .ge. 100) then
      Bch = 0
    else if (ityp > 0) then
      Bch = 1
    else
      Bch = -1
    endif 
    return

end function BfromItyp

!=========================================================================================================
integer function SfromItyp(ityp) result (Sch)
 implicit none
 integer, intent(in) :: ityp
 integer ia

 integer :: strres(1:55), strmes(100:137)
 data strres/ 26*0,13*1,9*1,6*2,3/
 data strmes/ 6*0,-1,0,-1,0,-1,0,0,-1,0,0,0,-1,0,0,0,-1,0,0,0,-1,0,0,0,-1,0,0,0,0,0,0,0,0/

  !return zero strangeness in case of PYTHIA PDG particles 
  if (abs(ityp).gt.1000)then
      Sch=0
      return
  endif

  !standart UrQMD particles
  if(ityp.eq.0) then;  Sch=0; return;  endif
  ia=iabs(ityp)
  if(ia .ge. 100)then;  Sch = strmes(ia);  else;   Sch = strres(ia);  endif
  if (ityp < 0) then; Sch = - Sch; endif
  return
  
end function SfromItyp



end module read_f14f15
