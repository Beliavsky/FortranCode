module string_funcs_mod
implicit none
private
public :: count_char
contains
elemental function count_char(s, c) result(ncount)
! count the number of times character c appears in string s
    character(len=*), intent(in) :: s
    character, intent(in) :: c
    integer :: ncount, i, len_trim_s
    len_trim_s = len_trim(s)
    ncount = 0
    do i = 1, len_trim_s
       if (s(i:i) == c) ncount = ncount + 1
    end do
    if (c == ' ') ncount = ncount + len(s) - len_trim_s
end function count_char
end module string_funcs_mod