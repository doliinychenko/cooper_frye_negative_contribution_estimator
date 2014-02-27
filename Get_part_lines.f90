program GetPartLines
!Purpose: 
!      1) Read output by transport code (urqmd file14+file15) event by event
!      2) Create a structure of straight line pieces of particle trajectories
!      3) Write the structure to binary file (set of binary files with max_events_per_output)
!Uses: 
!      read_f14_f15 module
!Author:
!      Dmytro Oliinychenko (oliiny@fias.uni-frankfurt.de)

use read_f14f15, only: gen_part, read_f15_event_header, read_f14_f15_event

 implicit none

 integer, parameter :: max_lines_per_event   = 25000
 integer, parameter :: max_events_per_output = 500
 integer, parameter :: maxnevt = 4500
 integer, parameter :: fnum_start = 0
 character(LEN=*), parameter :: urqmdf14_input = "/tmp/Tmn_proj_data/UrQMD_f14f15/E10b0ev4500.f14"
 character(LEN=*), parameter :: urqmdf15_input = "/tmp/Tmn_proj_data/UrQMD_f14f15/E10b0ev4500.f15"
 character(LEN=*), parameter :: output_prename = "/tmp/Tmn_proj_data/TrajLines/"
 character(LEN=100) :: output_name
 logical op

 integer fout_iter, ev, Np_evt, i
 !type (gen_part) :: hlp_part
 type (gen_part), dimension(max_lines_per_event) :: part_evt
 double precision Elab, b_imp


 print *,"Transforming UrQMD files to binary structure of straight lines of particle trajectories."
 print *,"Reading from: UrQMD file 14 - ",urqmdf14_input
 print *,"Reading from: UrQMD file 15 - ",urqmdf15_input

 open(unit=8, file=urqmdf15_input)
 open(unit=9, file=urqmdf14_input)

 call read_f15_event_header(8, Elab, b_imp)
 write(*,'(A,f8.2,2X,A,f8.2)')"Elab[AGeV]: ",Elab,"b = ", b_imp

 fout_iter = fnum_start
 ev=0
 Np_evt = -1

 do while (Np_evt .ne. 0 .and. ev < maxnevt ) !event cycle 

   !read event from urqmd
   call read_f14_f15_event(9, 8, Np_evt, part_evt) 
   write(*,'(A,I5,A,I5,A)',advance="no")"event: ", ev, " Np: ", Np_evt, char(13);   
   
   !write to next file if limit exceeded for previous
   if (mod(ev,max_events_per_output) == 0 ) then
    inquire(unit = 10, opened = op)
    if (op) then; close(10); endif
    call make_out_name(output_prename, Elab, b_imp, fout_iter, output_name)    
    open(unit=10, file=output_name, form='UNFORMATTED')
    print *,"Event: ",ev, "   Output to file: ",output_name
    fout_iter = fout_iter + 1    
   endif

   !Dump lines array to file
   write(10)Np_evt
   do i=1, Np_evt
    write(10)part_evt(i)
   end do

   ev = ev + 1
 
 end do !end event cycle

 print *,"Transfortmation finshed. ", ev, " events transformed."

 close(8)
 close(9) 
 close(10)

contains

!===========================================================================================================
subroutine make_out_name(prename, E, b, cnt, out_name)
implicit none
 character(LEN=*), intent(in) :: prename 
 double precision, intent(in) :: E, b     !Collision energy[GeV], centrality[fm]
 integer, intent(in) :: cnt                !file marker
 character(LEN=*), intent(out) :: out_name 
 character(LEN=10) :: Estr, bstr, nstr

 write(Estr,'(F10.1)')E
 write(bstr,'(F10.1)')b
 write(nstr,'(I10)')cnt

 out_name = prename//"E"//trim(adjustl(Estr))//"_b"//trim(adjustl(bstr))//"_n"//trim(adjustl(nstr))//".dat"

end subroutine

end program GetPartLines
