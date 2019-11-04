module linal
  use sort
  implicit none

contains

subroutine linal_dgeev(N,H,R_eval,I_eval,Psi)
  implicit none
  real(kind=8), dimension(0:,0:), intent(inout) :: H,Psi
  real(kind=8), dimension(0:), intent(inout) :: R_eval,I_eval
  integer, intent(in) :: N

  real(kind=8), dimension(:), allocatable :: WORK
  real(kind=8), dimension(0:1,0:1) :: VL
  character(len=1) :: JOBVL,JOBVR
  integer :: LDA,LDVL,LDVR,LWORK,INFO

  LWORK = -1
  JOBVL = 'N'
  JOBVR = 'V'
  LDA = N
  LDVL = 1
  LDVR = N 

  write(*,*) "Diagonalizing the Hamiltonian"
  allocate(WORK(0:1))
  call DGEEV(JOBVL,JOBVR,N,H(0:N-1,0:N-1),LDA,R_eval,I_eval,VL,LDVL,Psi,LDVR,WORK,LWORK,INFO)
  LWORK = INT(WORK(0))+1
  deallocate(WORK)
  allocate(WORK(0:LWORK-1))
  call DGEEV(JOBVL,JOBVR,N,H(0:N-1,0:N-1),LDA,R_eval,I_eval,VL,LDVL,Psi,LDVR,WORK,LWORK,INFO)
  write(*,*) "Hamiltonian diagonalized with info",INFO

  call sort_eval(N,R_eval,I_eval,Psi)
end subroutine linal_dgeev

!------------------------------------------------------------
! linal_dsyevx_lwork
!       - calculates the needed length of working vector,
!         lwork
!------------------------------------------------------------
! N             : int, size of matrix
! enum          : int, number of eigenvalues to find
! lwork         : int, length of lwork
! error         : int, exit code

SUBROUTINE linal_dsyevx_lwork(N,enum,lwork,error)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: lwork,error
  INTEGER, INTENT(IN) :: N,enum
  REAL(KIND=8), DIMENSION(0:1,0:1) :: A
  REAL(KIND=8), DIMENSION(0:1) :: W,WORK,Z
  INTEGER, DIMENSION(0:1) :: IWORK,IFAIL
  CHARACTER(LEN=1) :: JOBZ,RANG,UPLO
  REAL(KIND=8) :: VL,VU,ABSTOL
  INTEGER :: LDA,IL,IU,LDZ,INFO,M
  error = 0

  JOBZ = 'V'
  RANG = 'I'
  UPLO = 'U'
  LDA = N
  VL = 0.0
  VU = 0.0
  IL = 1
  IU = enum
  ABSTOL = 1.0D-15
  LDZ = N
  lwork = -1

  CALL DSYEVX(JOBZ,RANG,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,&
              M,W,Z,LDZ,WORK,lwork,IWORK,IFAIL,INFO)
  IF (INFO .NE. 0) THEN
    WRITE(*,*) "linal_dsyevx_lwork  : ERROR"
    WRITE(*,*) "Something went wrong while determining the "
    WRITE(*,*) "  length of working vectors"
    error = 1
    RETURN
  END IF

  lwork = CEILING(WORK(0))

END SUBROUTINE linal_dsyevx_lwork

!------------------------------------------------------------
! linal_dsyevx
!       - diagonalizes matrix up to a certain eigenvalue
!       - assumes lower triangular matrix
!------------------------------------------------------------
! N             : int, size of matrix
! IU            : int, upper integer
! LWORK         : int, length of work array
! A             : 2D real*8, symmetrix matrix to diagonalize
! W             : 1D real*8, eigenvalues
! Z             : 2D real*8, eigenvectors
! error         : int, exit code

SUBROUTINE linal_dsyevx(N,IU,LWORK,A,W,Z,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: A,Z
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: W
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,IU,LWORK

  REAL(KIND=8), DIMENSION(0:LWORK-1) :: WORK
  INTEGER, DIMENSION(0:5*N-1) :: IWORK
  INTEGER, DIMENSION(0:N-1) :: IFAIL
  CHARACTER(LEN=1) :: JOBZ,RANG,UPLO
  REAL(KIND=8) :: VL,VU,ABSTOL
  INTEGER :: LDA,IL,M,LDZ,INFO
  INTEGER :: i,j

  error = 0
  JOBZ = 'V'
  RANG = 'I'
  UPLO = 'U'
  LDA = N
  VL = 0.0
  VU = 0.0
  IL = 1
  ABSTOL = 1.0D-15
  LDZ = N

  CALL DSYEVX(JOBZ,RANG,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,&
              M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
  IF (INFO .LT. 0) THEN
    WRITE(*,*) "linal_dsyevx  : ERROR"
    WRITE(*,*) "This input had a bad value ",INFO
    error = 1
  END IF
  IF (INFO .GT. 0) THEN
    WRITE(*,*) "linal_dsyevx  : ERROR"
    WRITE(*,*) "The following eigenvalues failed to converge"
    WRITE(*,*) IFAIL
    error = 1
  END IF

  !Standardize the eigenvectors
  !DO j=0,IU-1
  !  i = MAXLOC(ABS(Z(0:N-1,j)),1)-1
  !  IF (Z(i,j) .LT. 0.0D0) Z(0:N-1,j) = -1.0D0*Z(0:N-1,j)
  !END DO

  open(file='eigs.dat',unit=99,status='replace')
  OPEN(FILE='evec.dat',UNIT=100,STATUS='REPLACE')
  DO i=0,IU-1
    write(99,*) i,W(i)
    WRITE(100,*) "Eigenvalue #",i,W(i)
    DO j=0,N-1
      WRITE(100,*) j,Z(j,i)
    END DO
    WRITE(100,*) "------------------------------------------------------"
  END DO
  CLOSE(UNIT=100)
  close(unit=99)

END SUBROUTINE linal_dsyevx

end module linal
