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
  integer :: i
  pi = 3.14159265358979D0
  
  !potential type -- probably remove this later once testing is done
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

  allocate(lb(0:ndim-1))
  allocate(ub(0:ndim-1))
  allocate(delx(0:ndim-1))
  allocate(npoints(0:ndim-1))
  allocate(coord(0:ndim-1))

  !potential energy cutoff
  if (pot .eq. 1 .or. pot .eq. 4) then
    write(*,*) "Enter V cutoff..."
    read(*,*) Vc
    write(*,*)
  else
    Vc = 0.1D0
  end if

  write(*,*) "Enter information for each coordinate"
  !coordinate information
  if (pot .eq. 4 .or. pot .eq. -4) then
    write(*,*) "Coordinate Types:"
    write(*,*) "1 -> -inf, +inf   (cartesian)"
    write(*,*) "2 ->    0, +inf   (radial)"
    write(*,*) "3 ->    0,    π   (polar)"
    write(*,*) "4 ->    0,   2π   (azimuthal, periodic)"
    write(*,*) "5 ->    a,    b   (box)"
    write(*,*) 

    do i=0,ndim-1
      write(*,'(1x,A9,1x,I3)') "dimension",i+1
      write(*,*) "Enter coordinate type"
      read(*,*) coord(i)

      !cartesian
      if (coord(i) .eq. 1) then
        write(*,*) "type : cartesian"
        write(*,*) "Enter gridpoint spacing"
        read(*,*) delx(i)
        write(*,*) "Enter maximum points"
        read(*,*) npoints(i)
        lb(i) = -1.0D0*delx(i)*(npoints(i)/2)
        ub(i) = delx(i)*(npoints(i)/2)

      !radial
      else if (coord(i) .eq. 2) then
        write(*,*) "type : radial"
        write(*,*) "Enter gridpoint spacing"
        read(*,*) delx(i)
        write(*,*) "Enter maximum points" 
        read(*,*) npoints(i)
        lb(i) = 0.0d0
        ub(i) = delx(i)*npoints(i)

      !polar
      else if (coord(i) .eq. 3) then
        write(*,*) "type : polar"
        write(*,*) "Enter number of grid points"
        read(*,*) npoints(i)
        lb(i) = 0.0D0
        ub(i) = pi 
        delx(i) = pi/npoints(i)

      !azimuthal
      else if (coord(i) .eq. 4) then
        write(*,*) "type : azimuthal"
        write(*,*) "Enter number of grid points"
        read(*,*) npoints(i)
        lb(i) = 0.0D0
        ub(i) = 2.0D0*pi
        delx(i) = 2.0D0*pi/npoints(i) 

      !box
      else if (coord(i) .eq. 5) then
        write(*,*) "type : box"
        write(*,*) "Enter lower bound" 
        read(*,*) lb(i)
        write(*,*) "Enter upper bound"
        read(*,*) ub(i)
        write(*,*) "Enter number of gridpoints"
        read(*,*) npoints(i) 
        delx(i) = (ub(i) - lb(i))/npoints(i)

      !bad
      else 
        write(*,*) "Sorry, that coordinate type is not supported"
        stop 
      end if 
     
      write(*,*)
    end do
  end if

  !print output
  write(*,*) "Grid information for each system"
  write(*,'(1x,A9,4x,A4,6x,A11,5x,A11,9x,A3)') "dimension","type","lower bound","upper bound","Δx"    
  do i=0,ndim-1
    write(*,'(1x,6x,I3,6x,I2,5x,ES12.5,4x,ES12.5,4x,ES12.5)') i+1,coord(i),lb(i),ub(i),delx(i)
  end do
  write(*,*) 

end subroutine
!---------------------------------------------------------------------
  

end module input
