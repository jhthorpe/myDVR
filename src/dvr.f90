! H		: 2D real*8, hamiltonian
! npoints	: 1D int, number of gridnpoints per dim
! Vc		: real*8, potential cutoff
! delx		: 1D real*8, gridpoint size
! ub		: 1D real*8, upper bounds
! lb		: 1D real*8, lower bounds
! m		: 1D real*8, mass 
! N		: int, number of untrimmed points
! Np		: int, number of trimmed points
! id_vec	: 1D int, list of original gridpoint numbers

program dvr
  use input
  use V
  use T
  use points
  use linal
  use fprint
  use crd

  implicit none
  real(kind=8), dimension(:,:), allocatable :: H,Psi  
  real(kind=8), dimension(:), allocatable :: lb,ub,delx,eval,m
  integer, dimension(:), allocatable :: npoints,id_vec,coord
  real(kind=8) :: Vc
  integer :: ndim,pot,N,Np,i,neig,lwork,sys,error

  !Get input and coordinate system
  call input_get(ndim,lb,ub,delx,m,pot,coord,npoints,Vc)
  call crd_get(ndim,coord,sys,error)
  if (error .ne. 0) then
    call dvr_clean(lb,ub,delx,eval,npoints,id_vec,coord,H,Psi)
    stop 1
  end if

  !Get size of untrimmed Hamiltonian
  N = 1
  do i=0,ndim-1
    N = (npoints(i))*N
  end do
  write(*,*) "Full Hamiltonian is dimension",N
  write(*,*)

  !generate points and hamiltonian
  call points_full(ndim,N,npoints,delx,lb,ub,coord) 
  if (pot .eq. -4) stop
  call V_calc(ndim,npoints,delx,lb,ub,N,pot,coord,Vc,Np,id_vec,H)
  call points_trim(ndim,Np,npoints,delx,lb,coord,id_vec,H) 
  call T_calc(ndim,npoints,sys,delx,lb,ub,m,coord,Np,id_vec,H,error)
  if (error .ne. 0) then
    call dvr_clean(lb,ub,delx,eval,npoints,id_vec,coord,H,Psi)
    stop 1
  end if

  !get desired number of eigenvalues
  allocate(eval(0:Np-1))
  allocate(Psi(0:Np-1,0:Np-1))
  write(*,*) "Enter number of eigenvalues"
  read(*,*) neig
  if (neig .gt. Np) neig = Np
  call linal_dsyevx_lwork(Np,neig,lwork,error)
  call linal_dsyevx(Np,neig,lwork,H,eval,Psi,error)
 
  !print out results
  write(*,*) "Eigenvalues"
  write(*,*) "        QN      λ"
  do i=0,neig-1
    write(*,*) i,eval(i)
  end do
  write(*,*)
  call fprint_espc(ndim,lb,ub,neig,eval)
  call fprint_wave(ndim,Np,npoints,delx,lb,coord,id_vec,neig,eval,Psi)
  Vc = 0.0D0
  do i=0,Np-1
    Vc = Vc + Psi(i,0)**2.0D0 
  end do
  write(*,*) "Psi^2 of ground state is", Vc

contains

!---------------------------------------------------------------------
! dvr_clean
!	- cleans up memory 
!---------------------------------------------------------------------
subroutine dvr_clean(lb,ub,delx,eval,npoints,id_vec,coord,H,Psi)
  implicit none
  real(kind=8), dimension(:,:), allocatable, intent(inout) :: H,Psi  
  real(kind=8), dimension(:), allocatable, intent(inout) :: lb,ub,&
                                                            delx,eval
  integer, dimension(:), allocatable, intent(inout) :: npoints,&
                                                       id_vec,coord
  if (allocated(H)) deallocate(H)
  if (allocated(Psi)) deallocate(Psi)
  if (allocated(lb)) deallocate(lb)
  if (allocated(ub)) deallocate(ub)
  if (allocated(delx)) deallocate(delx)
  if (allocated(eval)) deallocate(eval)
  if (allocated(npoints)) deallocate(npoints)
  if (allocated(id_vec)) deallocate(id_vec)
  if (allocated(coord)) deallocate(coord)

end subroutine dvr_clean
!---------------------------------------------------------------------

end program dvr
