program xstring_funcs
use string_funcs_mod
implicit none
print "(*(1x,i0))", count_char("hello world ", ["l", "a", " "])
end program xstring_funcs
