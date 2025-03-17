
// Factorial function using recursion
unsigned long long factorial(int n);

// Function to calculate Bessel function of the first kind J_n(x) using series expansion
double bessel_J(int n, double x);

inline int msta1(double x, int mp);

inline int msta2(double x, int n, int mp);

template <typename T> void jynbh(int n, int nmin, T x, int *nm, T *bj, T *by);


inline void jyndd(int n, double x, double *bjn, double *djn, double *fjn, double *byn, double *dyn, double *fyn);


inline void bessel_zeros(int n, int nt, double *rj0, double *rj1, double *ry0, double *ry1);