module Histogram_module
implicit none

integer, parameter :: max_hist_size=3500

type histo
sequence
 double precision :: l,u
 integer :: Nbin
 double precision, dimension(-max_hist_size-1 : max_hist_size+1) :: h
end type histo

contains
!==========================================================================
subroutine HistNew(hist,Nbin,l,u)
!Creates empty histogram hist with bounds l,u and number of bins Nbin
implicit none
   type (histo), intent(inout) :: hist
   double precision, intent(in) :: l,u
   integer, intent(in) :: Nbin

   hist%l = l; hist%u = u; hist%Nbin = Nbin 
   hist%h = 0.d0;
end subroutine
!==========================================================================
integer function ValueToBin(hist,val) result (bin)
!Finds proper bin for val in histogram
implicit none
 type (histo), intent(in) :: hist
 double precision, intent(in) :: val

 if (hist%u.le.hist%l) then; print *,"Error: Upper bound - ",hist%u," should be larger than lower bound."; call exit(); endif
 bin=nint( (val - 0.5d0*(hist%l+hist%u))/HistBinSize(hist) )
 if (bin>hist%Nbin) then; bin=hist%Nbin+1; endif
 if (bin<-hist%Nbin) then; bin=-hist%Nbin-1; endif
 
end function ValueToBin
!==========================================================================
double precision function BinToValue(hist,bin) result (val)
!Finds value val of bin in histogram hist
implicit none
 type (histo), intent(in) :: hist
 integer, intent(in) :: bin

 if (hist%u.le.hist%l) then; print *,"Error: Upper bound - ",hist%u," should be larger than lower bound."; call exit(); endif
 val = bin*HistBinSize(hist)+ 0.5d0*(hist%l+hist%u)

end function BinToValue

!==========================================================================
double precision function HistBinSize(h) result (bsize)
!returns bin size of the histogram h
 implicit none
 type (histo), intent(in) :: h
 bsize = (h%u-h%l) *0.5d0 / h%Nbin
end function HistBinSize

!==========================================================================
subroutine HistAdd(hist,val,w)
!Adds w to hist bin corresponding to val
implicit none
  type (histo), intent(inout) :: hist
  double precision, intent(in) :: val,w
  integer i
 
  i=ValueToBin(hist,val)
  hist%h(i) = hist%h(i) + w
end subroutine
!==========================================================================
subroutine HistPrint(hist,fname)
!Prints histogram hist to file fname
implicit none
 type (histo), intent(in) :: hist
 character(len=*), intent(in) :: fname
 integer i

 open(unit=75, file=trim(adjustl(fname)))
  write(75,'(A,2f12.3)')"#histogram, range :",hist%l-0.5d0*HistBinSize(hist), &
                                               hist%u+0.5d0*HistBinSize(hist)
  write(75,'(A,i5)')    "#Number of bins: ", 2*hist%Nbin
  write(75,'(A,f18.9)') "#Below histogram range:", hist%h(-hist%Nbin-1)
  do i=-hist%Nbin,hist%Nbin
   write(75,'(2f18.9)')BinToValue(hist,i), hist%h(i)
  end do
  write(75,'(A,f18.9)') "#Above histogram range:", hist%h(hist%Nbin+1)
  write(75,'(A,f18.9)') "#Integral is (only in range of histogram): ", HistInt(hist)
 close(75)
end subroutine
!==========================================================================
subroutine HistSumup(h1,h2)
!Performs h1=h1+h2 for histos
implicit none
type (histo), intent(in) :: h2
type (histo), intent(inout) :: h1
 if (h1%l == h2%l .and. h1%u == h2%u .and. h1%Nbin==h2%Nbin) then
  h1%h=h1%h + h2%h
 else
  print *,"Error: attempt to add histos with different binning or bounds."
 endif
end subroutine
!==========================================================================


double precision function HistInt(h) result (Hint)
!Calculates integral of histogram
implicit none
type (histo), intent(in) :: h
integer i
 
 Hint=0d0
 do i=-h%Nbin, h%Nbin
  Hint = Hint + h%h(i)*HistBinSize(h)
 end do

end function HistInt

end module Histogram_module
