program fix
  implicit none
  integer :: i,dummy,N
  real(kind=8) :: x,y,c,yold
  character(len=1024) :: fname
  write(*,*) "enter N"
  read(*,*) N
  write(*,*) "enter fname"
  read(*,*) fname
  open(file=trim(fname),unit=100)
  open(file='foo',unit=101)
  read(100,*) dummy,x,y,c
  write(101,*) dummy,x,y,c
  yold = y
  do i=1,N-1 
    read(100,*) dummy,x,y,c
    if ( y .ne. yold) write(101,*)
    write(101,*) dummy,x,y,c
    yold = y
  end do
  close(unit=100)
  close(unit=101)
end program fix
