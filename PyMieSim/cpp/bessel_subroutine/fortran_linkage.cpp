#pragma once

// namespace sp_bessel {
/*! @name Fortran Linkage
 * We link the Fortran subroutines to our C++ code.
 * Since Fortran treats all variables by reference, we must
 * pass our C++ arguments as pointers.
 *
 * Most subroutines have flags to signal errors such as underflow
 * and overflow. Exponential scaling is also available.
 * The KODE flag is a parameter that indicates the scaling option.
 * KODE = 1 means no scaling while KODE activates the exponential scaling.
 * Integer N controls the size of the array that is returned in the Fortran code.
 * We use N=1 throughout as to return scalars.  The NZ integer counts the number of
 * elements set to zero due to underflow. We ignore it. The IERR integer
 * is for error signaling.
 */
///@{
extern "C"
{
  /*! Bessel function of the first kind. */
  extern void zbesj_wrap(double, double, double, int, int, double*, double*, int*, int*);

  /*! Bessel function of the second kind. */
  extern void zbesy_wrap(double, double, double, int, int, double*, double*, int*, double*, double*, int*);

  /*! Hankel function of both kinds. Kind determined by integer argument. */
  extern void zbesh_wrap(double, double, double, int, int, int, double*, double*, int*, int*);
}
///@}

// } // namespace sp_bessel
