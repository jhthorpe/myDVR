program conv
  real(kind=8) :: rad,deg,pi,ax
  integer :: i,j
  pi = 3.14159265358979D0
  open(file='full_grid.txt',unit=100)
  open(file='V.in',unit=101)
  open(file='foo',unit=102)
  do i=0,99
    read(100,*) j,rad
    read(101,*) j,deg,ax
    write(102,*) j,rad,ax
  end do
  close(unit=100)
  close(unit=101)
  close(unit=102)
end program conv
