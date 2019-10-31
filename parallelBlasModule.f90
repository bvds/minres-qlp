!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     File minresqlpBlasModule.f90
!
!     This file contains the following BLAS subroutines
!        mpi_dot, mpi_dnrm2
!     required by subroutine MINRESQLP.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module parallelBlasModule

  use mpi
  use minresqlpDataModule,    only : dp, ip, zero, one
  implicit none

  public   :: mpi_dot, mpi_dnrm2
#ifdef USE_MKL
  include 'mkl_blas.fi'
#else
  !
  !  from "cblas.h"
  !  double cblas_ddot(f77_int N, const double *X, f77_int incX,
  !                const double *Y, f77_int incY);
  interface
     double precision function ddot(N, X, incX, Y, incY)
       use iso_c_binding
       integer :: N, incX, incY
       double precision :: X(N), Y(N)
     end function ddot
  end interface
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! mpi_dot forms the dot product of two vectors.
!
! This is intended to be an MPI extension of the Fortran
! function dot_product

    real(dp) function mpi_dot(n, dx, dy)

      implicit none
      integer(ip), intent(in) :: n
      real(dp),    intent(in) :: dx(*), dy(*)

      integer     :: ierr
      real(dp)    :: dtemp

      dtemp = ddot(n, dx, 1, dy, 1)

      call mpi_allreduce(dtemp, mpi_dot, 1, MPI_DOUBLE, MPI_SUM, &
           MPI_COMM_WORLD, ierr )
      return
end function mpi_dot

!*****************************************************************************
!
!! mpi_dnrm2 returns the euclidean norm of a vector using mpi.
!

real(dp) function mpi_dnrm2 ( n, x, incx )

  implicit none
  integer(ip), intent(in) :: n, incx
  real(dp),    intent(in) :: x(*)
  intrinsic :: sqrt

  real(dp)                :: dtemp, gnorm2
  integer                 :: ierr

  dtemp = ddot(n, x, incx, x, incx)
  call mpi_allreduce(dtemp, gnorm2, 1, MPI_DOUBLE, MPI_SUM, &
       MPI_COMM_WORLD, ierr)
  mpi_dnrm2 = sqrt(gnorm2)
  return
end function mpi_dnrm2

end module parallelBlasModule
