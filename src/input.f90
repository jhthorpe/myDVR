module input
  implicit none

contains

subroutine input_get(ndim,box,delx,pot,npoints,Vc)
  implicit none
  real(kind=8), dimension(:), allocatable, intent(inout) :: box,delx
  integer, dimension(:), allocatable, intent(inout) :: npoints
  real(kind=8), intent(inout) :: Vc
  integer, intent(inout) :: ndim,pot
 
  integer :: i
  
  write(*,*) "Enter potential type"
  write(*,*) "Options:"
  write(*,*) "1 -> Harmonic Oscillator"
  write(*,*) "2 -> Particle in a Box" 
  write(*,*) "3 -> Particle in a Stadium"
  write(*,*) "4 -> General potential (-4 to generate)"
  read(*,*) pot
  write(*,*)
  
  write(*,*) "Enter number of dimensions"
  read(*,*) ndim
  write(*,*)

  allocate(box(0:ndim-1))
  allocate(delx(0:ndim-1))
  allocate(npoints(0:ndim-1))

  if (pot .eq. 1 .or. pot .eq. 4) then
    write(*,*) "Enter V cutoff..."
    read(*,*) Vc
    write(*,*)
  else
    Vc = 0.1D0
  end if

  write(*,*) "Enter bounding box..."
  do i=0,ndim-1
    write(*,*) "dimension ",i+1
    read(*,*) box(i) 
  end do
  write(*,*)

  write(*,*) "Enter number of gridpoints..." 
  do i=0,ndim-1
    write(*,*) "dimension ",i+1
    read(*,*) npoints(i)
    npoints(i) = 2*npoints(i)
  end do
  write(*,*)

  write(*,*) "Using the following gridsize"
  do i=0,ndim-1
    delx(i) = 2*box(i)/npoints(i)
    write(*,*) "dimension",i+1,delx(i)
  end do
  npoints = npoints+1
 
  write(*,*) "Number of gridpoints per dimension"
  do i=0,ndim-1
    write(*,*) npoints(i)
  end do 
  write(*,*)

end subroutine
  

end module input
