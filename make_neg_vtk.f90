program vtk_neg
! Reads points t,x,y,z
!makes vtk file for each [t] with unstructured point x,y,z
real, dimension(5000000) :: t,x,y,z
integer it(5000000), ni(200)
integer i,j, nlines,io
integer tstart_step, tend_step
real dt, dx, dy, dz
character*100 file_in,s_hlp
real tmax


dt = 0.3d0
dz = 0.3d0
dx = 1.d0
dy = 1.d0
tstart_step = 10
tend_step = 90
ni=0

file_in="/home/oliiny/Work/Tmn_project/hyd_vs_prt_git/output/E160.0_b0.0prod_e0_0.3/prt/neg_cr.txt"
print *,"Reading file ",trim(adjustl(file_in))," ..."
open(unit=8,file=trim(adjustl(file_in)))
 i=1
 tmax=0E0
 do 
  read(8,*,iostat=io) t(i),x(i),y(i),z(i)
  if (io.ne.0) then; exit; endif
  if (t(i)>tmax) then; tmax=t(i); endif
  it(i) = nint(t(i)/dt)
  ni(it(i)) = ni(it(i)) + 1
  i=i+1
  !print *,i
 end do
 nlines=i-1
close(8)
print *,nlines," lines read. Now converting to vtk."
print *,"ni: ", (ni(i),i=1,60)

do i=tstart_step, ceiling(tmax/dt)
 write(s_hlp,*)i;
 open(unit=8+i, file="/home/oliiny/Work/Tmn_project/hyd_vs_prt_git/output/E160.0_b0.0prod_e0_0.3/prt/vtk/neg.vtk."&
                 //trim(adjustl(s_hlp)) )
  write(8+i,'(A)')"# vtk DataFile Version 2.0"
  write(8+i,'(A)')"part_in"
  write(8+i,'(A)')"ASCII"
  write(8+i,'(A)')"DATASET UNSTRUCTURED_GRID"
  write(8+i,'(A,I8,A)')"POINTS ", ni(i), " float"
end do

do j=1, nlines;   write(8+it(j),'(3f10.2)') x(j)/dx,y(j)/dy,z(j)/dz;  end do

do i=tstart_step, ceiling(tmax/dt)
  write(8+i,'(A,2I8)')"CELLS ", ni(i), 2*ni(i)
  do j=1,ni(i); write(8+i,'(2I8)')1,j-1;  end do
  write(8+i,'(A,I8)')"CELL_TYPES", ni(i)
  do j=1,ni(i); write(8+i,'(I8)')1;  end do
end do


do i=tstart_step, ceiling(tmax/dt)
 close(8+i)
end do


end program vtk_neg
