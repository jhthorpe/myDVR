program pot
  implicit none
  integer :: N,i,num
  real(kind=8) :: x,y,ql,qr,kl,kr,ky
  real(kind=8) :: V,T,D,g,Vr,a,b
  
  !ql and qr are the locations of the left and right oscillators in x dimension
  !T is the trace of the 2x2 matrix
  !D is the determinate of the 2x2 matrix
  !g is the coupling between ql and qr
  !Vr is the potential offset of the right oscillator
 
  write(*,*) "enter ql"
  read(*,*) ql
  write(*,*) "enter qr"
  read(*,*) qr
  write(*,*) "enter kl"
  read(*,*) kl 
  write(*,*) "enter kr"
  read(*,*) kr 
  write(*,*) "enter Vr"
  read(*,*) Vr
  write(*,*) "enter g"
  read(*,*) g
  write(*,*) "enter ky"
  read(*,*) ky
  write(*,*) "enter number of lines"
  read(*,*) N
  
  open(file='full_grid.txt',unit=100,status='old')
  open(file='V.in',unit=101,status='replace')
  do i=0,N-1
    read(100,*) num,x,y
    a = 0.5D0*kl*(x - ql)**2.0D0
    b = 0.5D0*kr*(x - qr)**2.0D0 + Vr
    T = a + b
    D = a*b - g*g
    V = 0.5D0*T - sqrt(0.25d0*T*T - D) + 0.5D0*ky*y*y
    write(101,*) num,x,y,V
  end do
  close(unit=100)
  close(unit=101) 
  

end program pot
