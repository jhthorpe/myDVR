program fixer
  implicit none
  integer :: i,j
  real(kind=8) :: x,energy
  open(file='V_temp.txt',unit=100)
  open(file='full_grid.txt',unit=101)
  open(file='V.in',unit=102)
  do i=0,198
    read(100,*) j,energy
    read(101,*) j,x
    write(102,*) j,x,energy
  end do
  close(unit=100)
  close(unit=101)
  close(unit=102)

end program fixer
