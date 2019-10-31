module input
  implicit none

contains


!---------------------------------------------------------------------
! input_get
!	- reads input for dvr program, and generates Δx and bounding
!	  conditions
!---------------------------------------------------------------------
! nimd		: int, number of dimensions
! lb		: 1D real*8, lower bound for each dim
! ub		: 1D real*8, upper bound for each dim
! delx		: 1D real*8, delta x for each dim
! pot		: int, potential type for the system
! coord		: 1D int, coordinate types for each dim
! npoints	: 1D int, number of points for each dim
! Vc		: real*8, potential energy cutoff for trimming

subroutine input_get(ndim,lb,ub,delx,pot,coord,npoints,Vc)
  implicit none
  real(kind=8), dimension(:), allocatable, intent(inout) :: lb,ub,delx
  integer, dimension(:), allocatable, intent(inout) :: npoints,coord
  real(kind=8), intent(inout) :: Vc
  integer, intent(inout) :: ndim,pot
 
  real(kind=8) :: pi
  integer :: i,N
  pi = 3.14159265358979D0
  
  !potential type -- probably remove this later once testing is done
  write(*,*) "Enter potential type"
  write(*,*) "Options:"
  write(*,*) "1 -> Harmonic Oscillator"
  write(*,*) "2 -> Particle in a Box" 
  write(*,*) "3 -> Particle in a Stadium"
  write(*,*) "4 -> General potential (-4 to generate)"
  write(*,'(1x,A1,1x)',advance='no') ">"
  read(*,*) pot
  if ((pot .lt. 1 .and. pot .ne. -4) .or. pot .gt. 4) then
    write(*,*) "That was a bad potential, try again."
    stop
  end if
  write(*,*)
  
  write(*,*) "Enter number of dimensions"
  write(*,'(1x,A1,1x)',advance='no') ">"
  read(*,*) ndim
  write(*,*)

  allocate(lb(0:ndim-1))
  allocate(ub(0:ndim-1))
  allocate(delx(0:ndim-1))
  allocate(npoints(0:ndim-1))
  allocate(coord(0:ndim-1))

  !potential energy cutoff
  if (pot .eq. 1 .or. pot .eq. 4) then
    write(*,*) "Enter V cutoff..."
    write(*,'(1x,A1,1x)',advance='no') ">"
    read(*,*) Vc
    write(*,*)
  else
    Vc = 0.1D0
  end if

  !coordinate information
  write(*,*) "Enter information for each coordinate"
  if (pot .ne. 4) then
  !if (pot .eq. 4 .or. pot .eq. -4) then
    write(*,*) "Coordinate Types:"
    write(*,*) "1 -> -inf, +inf   (cartesian)"
    write(*,*) "2 ->    0, +inf   (radial)"
    write(*,*) "3 ->    0,    π   (polar)"
    write(*,*) "4 ->    0,   2π   (azimuthal, periodic)"
    write(*,*) "5 ->    a,    b   (box)"
    write(*,*) 

    do i=0,ndim-1
      write(*,'(1x,A9,1x,I3)') "dimension",i+1
      write(*,*) "------------------------"
      write(*,*) "Enter coordinate type"
      write(*,'(1x,A1,1x)',advance='no') ">"
      read(*,*) coord(i)

      !cartesian
      if (coord(i) .eq. 1) then
        write(*,*) "type : cartesian"
        write(*,*) "Enter gridpoint spacing"
        write(*,'(1x,A1,1x)',advance='no') ">"
        read(*,*) delx(i)
        write(*,*) "Enter maximum N"
        write(*,'(1x,A1,1x)',advance='no') ">"
        read(*,*) N 
        lb(i) = -1.0D0*delx(i)*(N)
        ub(i) = delx(i)*(N)
        npoints(i) = 2*N - 1

      !radial
      else if (coord(i) .eq. 2) then
        write(*,*) "type : radial"
        write(*,*) "Enter gridpoint spacing"
        write(*,'(1x,A1,1x)',advance='no') ">"
        read(*,*) delx(i)
        write(*,*) "Enter maximum N" 
        write(*,'(1x,A1,1x)',advance='no') ">"
        read(*,*) N 
        lb(i) = 0.0d0
        ub(i) = delx(i)*N
        npoints(i) = N - 1

      !polar
      else if (coord(i) .eq. 3) then
        write(*,*) "type : polar"
        write(*,*) "Enter N" 
        write(*,'(1x,A1,1x)',advance='no') ">"
        read(*,*) N 
        lb(i) = 0.0D0
        ub(i) = pi 
        delx(i) = pi/N
        npoints(i) = N - 1

      !azimuthal
      else if (coord(i) .eq. 4) then
        write(*,*) "type : azimuthal"
        write(*,*) "Enter N" 
        write(*,'(1x,A1,1x)',advance='no') ">"
        read(*,*) N 
        lb(i) = 0.0D0
        ub(i) = 2.0D0*pi
        delx(i) = 2.0D0*pi/(2*N+1) 
        npoints(i) = 2*N+1

      !box
      else if (coord(i) .eq. 5) then
        write(*,*) "type : box"
        write(*,*) "Enter lower bound" 
        write(*,'(1x,A1,1x)',advance='no') ">"
        read(*,*) lb(i)
        write(*,*) "Enter upper bound"
        write(*,'(1x,A1,1x)',advance='no') ">"
        read(*,*) ub(i)
        write(*,*) "Enter N" 
        write(*,'(1x,A1,1x)',advance='no') ">"
        read(*,*) N
        delx(i) = (ub(i) - lb(i))/N
        npoints(i) = N - 1

      !bad
      else 
        write(*,*) "Sorry, that coordinate type is not supported"
        stop 
      end if 
     
      write(*,*)
    end do
  end if

  !print output
  write(*,*) "Grid information for each coordinate"
  write(*,'(1x,A9,4x,A4,6x,A11,5x,A11,9x,A3)') "dimension","type","lower bound","upper bound","Δx"    
  do i=0,ndim-1
    write(*,'(1x,6x,I3,6x,I2,5x,ES12.5,4x,ES12.5,4x,ES12.5)') i+1,coord(i),lb(i),ub(i),delx(i)
  end do
  write(*,*) 
  call input_save(ndim,npoints,coord,lb,ub,delx)
  write(*,*) "Grid data saved to grid.dat"
  write(*,*) 

end subroutine

!---------------------------------------------------------------------
! input_save
!	- saves the input given by user in grid.dat so it can be
!	  used later
!---------------------------------------------------------------------
! ndim		: int, number of dimensions
! npoints	: 1D int, npoints in each dimension
! coord		: 1D int, coordinate type of each dimension
! lb		: 1D real*8, lower bounds of each dimension
! ub		: 1D real*8, upper bounds of each dimension  
! delx		: 1D real*8, Δx of each dimension

subroutine input_save(ndim,npoints,coord,lb,ub,delx)
  implicit none
  real(kind=8), dimension(0:), intent(in) :: lb,ub,delx
  integer, dimension(0:), intent(in) :: npoints,coord
  integer, intent(in) :: ndim
  integer :: i
  open(file='grid.dat',unit=100,status='replace')
  write(100,*) ndim
  do i=0,ndim-1
    write(100,*) i,coord(i),npoints(i),lb(i),ub(i),delx(i)
  end do
  close(unit=100)
end subroutine input_save
  
!---------------------------------------------------------------------

end module input
