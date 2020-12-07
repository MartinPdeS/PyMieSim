#include <vector>
#include <complex>
#include <boost/math/special_functions.hpp>
#include <cmath>


#if __has_include("python3.8/Python.h")
#include "python3.8/Python.h"
#elif _has_include("python3.6/Python.h")
#include "python3.6/Python.h"
#endif

typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;

static void
LowFrequencyMie_ab(const double m,
                        const double x,
                        iVec *an,
                        iVec *bn);



static void
HighFrequencyMie_ab(const double m,
                         const double x,
                         const long unsigned int nmax,
                         const std::vector<double>* n,
                         iVec *an,
                         iVec *bn);



static void
MiePiTau(const double mu,
              const long unsigned int nmax,
              iVec *pin,
              iVec *taun);




static void
C_GetS1S2(const double m,
          const double x,
          const double*  phi,
          const long unsigned int lenght,
          complex128* v);














// -
