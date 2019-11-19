program foo
  implicit none
  integer :: i,j,N
  real(kind=8) :: a,b,v,x

  write(*,*) "enter a"
  read(*,*) a
  write(*,*) "enter b"
  read(*,*) b
  write(*,*) "Enter N"
  read(*,*) N
  open(file='full_grid.txt',unit=100)
  open(file='V.in',unit=101)
  do i=0,N-1
    read(100,*) j,x
    v = a*x**2.0D0 + b*x**4.0D0
    write(101,*) j,x,v 
  end do
  close(unit=100)
  close(unit=101)

end program foo
