program gen
  implicit none
  integer :: i,j,N
  real(kind=8) :: om,phi,q
  
  write(*,*) "Enter omega"
  read(*,*) om
  phi = om
  write(*,*) "Enter number of points"
  read(*,*) N
  open(file='full_grid.txt',unit=100)
  open(file='V.in',unit=101)
  do i=0,N-1
    read(100,*) j,q
    write(101,*) j,q,0.5D0*om*q**2.0D0
  end do
  close(unit=100)
  close(unit=101)

end program gen
