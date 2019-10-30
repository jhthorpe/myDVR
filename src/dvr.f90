! H		: 2D real*8, hamiltonian
! npoints	: 1D int, number of gridnpoints per dim
! Vc		: real*8, potential cutoff
! delx		: 1D real*8, gridpoint size
! ub		: 1D real*8, upper bounds
! lb		: 1D real*8, lower bounds
! N		: int, number of untrimmed points
! Np		: int, number of trimmed points
! id_vec	: 1D int, list of original gridpoint numbers

program dvr
  use input
  use V
  use T
  use points
  use linal
  use wave

  implicit none
  real(kind=8), dimension(:,:), allocatable :: H,Psi  
  real(kind=8), dimension(:), allocatable :: lb,ub,delx,eval
  integer, dimension(:), allocatable :: npoints,id_vec,coord
  real(kind=8) :: Vc
  integer :: ndim,pot,N,Np,i,neig,lwork,error

  call input_get(ndim,lb,ub,delx,pot,coord,npoints,Vc)

  N = 1
  do i=0,ndim-1
    N = (npoints(i))*N
  end do

  call points_full(ndim,N,npoints,delx,lb,ub,coord) 
  stop
  if (pot .eq. -4) stop
  call V_calc(ndim,npoints,delx,N,pot,coord,Vc,Np,id_vec,H)
  call points_trim(ndim,Np,npoints,delx,coord,id_vec,H) 
  call T_calc(ndim,npoints,delx,coord,Np,id_vec,H)

  allocate(eval(0:Np-1))
  allocate(Psi(0:Np-1,0:Np-1))
  write(*,*) "Enter number of eigenvalues"
  read(*,*) neig
  if (neig .gt. Np) neig = Np
  call linal_dsyevx_lwork(Np,neig,lwork,error)
  call linal_dsyevx(Np,neig,lwork,H,eval,Psi,error)
 
  write(*,*) "Eigenvalues"
  write(*,*) "        QN      λ"
  do i=0,neig-1
    write(*,*) i,eval(i)
  end do
  write(*,*)

  call wave_print(ndim,Np,npoints,delx,coord,id_vec,neig,eval,Psi)

  Vc = 0.0D0
  do i=0,Np-1
    Vc = Vc + Psi(i,0)**2.0D0 
  end do
  write(*,*) "Psi^2 of ground state is", Vc

end program dvr
