! by urbanjost https://fortran-lang.discourse.group/t/speed-of-array-intrinsics/3251/12
module min_max_mod
implicit none
integer, parameter :: wp = kind(1.0d0)
contains

pure subroutine mymax_temp(max_val,max_ind,x)
real(kind=wp), intent(in)  :: x(:) 
real(kind=wp), intent(out) :: max_val 
integer, intent(out)       :: max_ind 

integer                    :: ix 
real(kind=wp)              :: temp 

max_val = -huge(0.d0)
max_ind = -1

do ix = 1,size(x)
   temp = x(ix)
   if (temp>max_val) then
       max_val = temp
       max_ind = ix
   endif
enddo

end subroutine mymax_temp

pure subroutine mymax_loc(max_val,max_ind,x)
real(kind=wp), intent(in)  :: x(:) 
real(kind=wp), intent(out) :: max_val 
integer, intent(out)       :: max_ind 

integer                    :: ix  

max_val = -huge(0.d0)
max_ind = 1

do ix = 1,size(x)
   if (x(ix)>x(max_ind)) then
       max_ind=ix         
   end if
end do

max_val=x(max_ind)

end subroutine mymax_loc

pure subroutine mymax_intrinsic(max_val,max_ind,x)
real(kind=wp), intent(in)  :: x(:) 
real(kind=wp), intent(out) :: max_val 
integer, intent(out)       :: max_ind 

max_ind = maxloc(x,dim=1)
max_val = x(max_ind)

end subroutine mymax_intrinsic

pure subroutine mymax_notemp(max_val,max_ind,x)
real(kind=wp), intent(in)  :: x(:) 
real(kind=wp), intent(out) :: max_val 
integer, intent(out)       :: max_ind 

integer                    :: nx, ix 

nx = size(x)

max_val = -huge(0.d0)
max_ind = -1

do ix = 1,nx
   if (x(ix)>max_val) then
       max_val = x(ix)
       max_ind = ix
   endif
enddo

end subroutine mymax_notemp

end module min_max_mod
!
program xmin_max
use min_max_mod
implicit none
integer                      :: n 
real(kind=wp), allocatable   :: x(:) 
real(kind=wp)                :: t(2) 
real(kind=wp)                :: max_val
integer                      :: max_ind 
character(len=:),allocatable :: title
character(len=*),parameter   :: fmt_cr = "(a20,a20,1x,f16.12,i9,f8.4)" 

n = 10**8
allocate(x(n))
call random_number(x)

!do i=1,4

   call setup('temp')
   call mymax_temp(max_val,max_ind,x)
   call printme()
   
   call setup('notemp')
   call mymax_notemp(max_val,max_ind,x)
   call printme()
   
   call setup('intrinsic')
   call mymax_intrinsic(max_val,max_ind,x)
   call printme()

   call setup('loc')
   call mymax_loc(max_val,max_ind,x)
   call printme()

!enddo

contains
subroutine setup(header)
character(len=*),intent(in) :: header
   title=header
   call cpu_time(t(1))
end subroutine setup

subroutine printme()
   call cpu_time(t(2))
   print fmt_cr,title,"max, loc:",max_val, max_ind,t(2)-t(1)
end subroutine printme

end program xmin_max
