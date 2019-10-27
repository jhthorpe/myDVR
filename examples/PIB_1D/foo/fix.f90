program fix
  integer :: N,i,dummy
  real(kind=8) :: off,x,V

  write(*,*) "N="
  read(*,*) N
  write(*,*) "off ="
  read(*,*) off
  open(file='wave_2.dat',unit=100)
  open(file='foo.out',unit=101)
  do i=0,N-1
    read(100,*) dummy,x,V
    write(101,*) x+off,V
  end do 
  close(unit=100)
  close(unit=101)
end program 
