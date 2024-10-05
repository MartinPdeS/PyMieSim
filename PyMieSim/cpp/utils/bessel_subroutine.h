
// Factorial function using recursion
unsigned long long factorial(int n) {
    return (n == 0 || n == 1) ? 1 : n * factorial(n - 1);
}

// Function to calculate Bessel function of the first kind J_n(x) using series expansion
double bessel_J(int n, double x) {
    const int num_terms = 20;  // Number of terms in the series expansion
    double sum = 0.0;
    double x_half_pow_n = std::pow(x / 2.0, n);  // (x/2)^n part of the series

    for (int m = 0; m < num_terms; ++m) {
        double term = std::pow(-1.0, m) * std::pow(x / 2.0, 2 * m) / (factorial(m) * factorial(m + n));
        sum += term;
    }

    return x_half_pow_n * sum;
}


inline int msta1(double x, int mp) {

    // ===================================================
    // Purpose: Determine the starting point for backward
    //          recurrence such that the magnitude of
    //          Jn(x) at that point is about 10^(-MP)
    // Input :  x     --- Argument of Jn(x)
    //          MP    --- Value of magnitude
    // Output:  MSTA1 --- Starting point
    // ===================================================

    int it, nn, n0, n1;
    double a0, f, f0, f1;

    a0 = fabs(x);
    n0 = (int)(1.1*a0) + 1;
    f0 = 0.5*log10(6.28*n0) - n0*log10(1.36*a0/n0)- mp;
    n1 = n0 + 5;
    f1 = 0.5*log10(6.28*n1) - n1*log10(1.36*a0/n1) - mp;
    for (it = 1; it <= 20; it++) {
        nn = n1 - (n1 - n0) / (1.0 - f0/f1);
        f = 0.5*log10(6.28*nn) - nn*log10(1.36*a0/nn) - mp;
        if (abs(nn-n1) < 1) { break; }
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn;
}

inline int msta2(double x, int n, int mp) {

    // ===================================================
    // Purpose: Determine the starting point for backward
    //          recurrence such that all Jn(x) has MP
    //          significant digits
    // Input :  x  --- Argument of Jn(x)
    //          n  --- Order of Jn(x)
    //          MP --- Significant digit
    // Output:  MSTA2 --- Starting point
    // ===================================================

    int it, n0, n1, nn;
    double a0, hmp, ejn, obj, f, f0, f1;

    a0 = fabs(x);
    hmp = 0.5*mp;
    ejn = 0.5*log10(6.28*n) - n*log10(1.36*a0/n);
    if (ejn <= hmp ) {
        obj = mp;
        n0 = (int)(1.1*a0) + 1;
    } else {
        obj = hmp + ejn;
        n0 = n;
    }
    f0 = 0.5*log10(6.28*n0) - n0*log10(1.36*a0/n0) - obj;
    n1 = n0 + 5;
    f1 = 0.5*log10(6.28*n1) - n1*log10(1.36*a0/n1) - obj;
    for (it = 1; it <= 20; it++) {
        nn = n1 - (n1 - n0) / (1.0 - f0/f1);
        f = 0.5*log10(6.28*nn) - nn*log10(1.36*a0/nn) - obj;
        if (abs(nn-n1) < 1) { break; }
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn;
}

template <typename T>
void jynbh(int n, int nmin, T x, int *nm, T *bj, T *by) {

    // =====================================================
    // Purpose: Compute Bessel functions Jn(x), Yn(x)
    // Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
    //          n --- Highest order of Jn(x) and Yn(x) computed  ( n ≥ 0 )
    //          nmin -- Lowest order computed  ( nmin ≥ 0 )
    // Output:  BJ(n-NMIN) --- Jn(x)   ; if indexing starts at 0
    //          BY(n-NMIN) --- Yn(x)   ; if indexing starts at 0
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 to calculate the starting
    //          point for backward recurrence
    // =====================================================

    int k, m, ky;
    T pi = 3.141592653589793;
    T r2p = 0.63661977236758;
    T bs, s0, su, sv, f2, f1, f;
    T bj0, bj1, ec, by0, by1, bjk, byk;
    T p0, q0, cu, t1, p1, q1, t2;

    T a[4] = { -0.0703125, 0.112152099609375, -0.5725014209747314, 0.6074042001273483e+01 };
    T b[4] = { 0.0732421875, -0.2271080017089844, 0.1727727502584457e+01, -0.2438052969955606e+02 };
    T a1[4] = { 0.1171875, -0.144195556640625, 0.6765925884246826, -0.6883914268109947e+01 };
    T b1[4] = { -0.1025390625, 0.2775764465332031, -0.1993531733751297e+01, 0.2724882731126854e+02 };

    *nm = n;
    if (x < 1.0e-100) {
        for (k = nmin; k <= n; k++) {
            bj[k - nmin] = 0.0;
            by[k - nmin] = -1.0e+300;
        }

        if (nmin == 0) { bj[0] = 1.0; }
        return;
    }

    if ((x <= 300.0) || (n > (int)(0.9 * x))) {
        // Backward recurrence for Jn
        if (n == 0) {
            *nm = 1;
        }
        m = msta1(x, 200);
        if (m < *nm) {
            *nm = m;
        } else {
            m = msta2(x, *nm, 15);
        }
        bs = 0.0;
        su = 0.0;
        sv = 0.0;
        f2 = 0.0;
        f1 = 1.0e-100;
        f = 0.0;

        for (k = m; k >= 0; k--) {
            f = 2.0*(k + 1.0)/x*f1 - f2;
            if ((k <= *nm) && (k >= nmin)) {
                bj[k - nmin] = f;
            }
            if (k == 2 * (int)(k / 2) && k != 0) {
                bs += 2.0 * f;
                su += pow(-1, k / 2) * f / k;
            } else if (k > 1) {
                sv += pow(-1, k / 2) * k / (k * k - 1.0) * f;
            }
            f2 = f1;
            f1 = f;
        }
        s0 = bs + f;

        for (k = nmin; k <= *nm; k++) {
            bj[k - nmin] /= s0;
        }
        // Estimates for Yn at start of recurrence
        bj0 = f1 / s0;
        bj1 = f2 / s0;
        ec = log(x / 2.0) + 0.5772156649015329;
        by0 = r2p * (ec * bj0 - 4.0*su/s0);
        by1 = r2p * ((ec - 1.0)*bj1 - bj0/x - 4.0*sv/s0);

        if (0 >= nmin) { by[0 - nmin] = by0; }
        if (1 >= nmin) { by[1 - nmin] = by1; }
        ky = 2;
    } else {
        // Hankel expansion
        t1 = x - 0.25*pi;
        p0 = 1.0;
        q0 = -0.125/x;

        for (k = 1; k <= 4; k++) {
            p0 += a[k - 1] * pow(x,-2*k);
            q0 += b[k - 1] * pow(x, -2*k - 1);
        }

        cu = sqrt(r2p / x);
        bj0 = cu * (p0*cos(t1) - q0*sin(t1));
        by0 = cu * (p0*sin(t1) + q0*cos(t1));

        if (0 >= nmin) {
            bj[0 - nmin] = bj0;
            by[0 - nmin] = by0;
        }

        t2 = x - 0.75*pi;
        p1 = 1.0;
        q1 = 0.375/x;

        for (k = 1; k <= 4; k++) {
            p1 += a1[k - 1] * pow(x, -2*k);
            q1 += b1[k - 1] * pow(x, -2*k - 1);
        }

        bj1 = cu * (p1*cos(t2) - q1*sin(t2));
        by1 = cu * (p1*sin(t2) + q1*cos(t2));

        if (1 >= nmin) {
            bj[1 - nmin] = bj1;
            by[1 - nmin] = by1;
        }

        for (k = 2; k <= *nm; k++) {
            bjk = 2.0*(k - 1.0)/x*bj1 - bj0;
            if (k >= nmin) { bj[k - nmin] = bjk; }
            bj0 = bj1;
            bj1 = bjk;
        }
        ky = 2;
    }
    // Forward recurrence for Yn
    for (k = ky; k <= *nm; k++) {
        byk = 2.0 * (k - 1.0) * by1 / x - by0;

        if (k >= nmin)
            by[k - nmin] = byk;

        by0 = by1;
        by1 = byk;
    }
}



inline void jyndd(int n, double x, double *bjn, double *djn, double *fjn, double *byn, double *dyn, double *fyn) {

    // ===========================================================
    // purpose: compute bessel functions jn(x) and yn(x), and
    //          their first and second derivatives
    // input:   x   ---  argument of jn(x) and yn(x) ( x > 0 )
    //          n   ---  order of jn(x) and yn(x)
    // output:  bjn ---  jn(x)
    //          djn ---  jn'(x)
    //          fjn ---  jn"(x)
    //          byn ---  yn(x)
    //          dyn ---  yn'(x)
    //          fyn ---  yn"(x)
    // routines called:
    //          jynbh to compute jn and yn
    // ===========================================================

    int nm = 0;
    double bj[2], by[2];

    jynbh(n+1, n, x, &nm, bj, by);
    // compute derivatives by differentiation formulas
    *bjn = bj[0];
    *byn = by[0];
    *djn = -bj[1] + n*bj[0]/x;
    *dyn = -by[1] + n*by[0]/x;
    *fjn = (n*n/(x*x) - 1.0)*(*bjn) - (*djn)/x;
    *fyn = (n*n/(x*x) - 1.0)*(*byn) - (*dyn)/x;
    return;
}



inline void bessel_zeros(int n, int nt, double *rj0, double *rj1, double *ry0, double *ry1) {

    // ======================================================
    // Purpose: Compute the zeros of Bessel functions Jn(x),
    //          Yn(x), and their derivatives
    // Input :  n  --- Order of Bessel functions  (n >= 0)
    //          NT --- Number of zeros (roots)
    // Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
    //          RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
    //          RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
    //          RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
    // Routine called: JYNDD for computing Jn(x), Yn(x), and
    //                 their first and second derivatives
    // ======================================================

    /*
     * SciPy Note:
     * See GH-18859 for additional changes done by SciPy for
     * better initial condition selection in Newton iteration
     */

    int L;
    double b, h, x, x0, bjn, djn, fjn, byn, dyn, fyn;
    const double pi = 3.141592653589793;
    // -- Newton method for j_{N,L}
    // initial guess for j_{N,1}
    if (n == 0) {
        x = 2.4;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = n + 1.85576*pow(n, 0.33333) + 1.03315/ pow(n, 0.33333);
    }
    // iterate
    L = 0;
L10:
    x0 = x;
    jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= bjn/djn;
    if (fabs(x - x0) > 1e-11) { goto L10; }

    L += 1;
    rj0[L - 1] = x;
    // initial guess for j_{N,L+1}
    if (L == 1) {
        if (n == 0) {
            x = 5.52;
        } else {
            // Expansion from https://dlmf.nist.gov/10.21#E32 and
            // coefficients from Olver 1951
            x= n + 3.24460 * pow(n, 0.33333) + 3.15824 / pow(n, 0.33333);
        }
    } else {
        // growth of roots is approximately linear (https://dlmf.nist.gov/10.21#E19)
        x = rj0[L - 1] + (rj0[L - 1] - rj0[L - 2]);
    }
    if (L <= (n + 10)) {
        jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
        h = atan(fabs(djn) / sqrt(fabs(fjn * bjn)));
        b = -djn / (bjn * atan(h));
        x -= (h - pi/2) / b;
    }

    if (L < nt) { goto L10; }

    // -- Newton method for j_{N,L+1}'
    if (n == 0) {
        x = 3.8317;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = n + 0.80861 * pow(n, 0.33333) + 0.07249 / pow(n, 0.33333);
    }
    // iterate
    L=0;
L15:
    x0 = x;
    jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= djn/fjn;
    if (fabs(x-x0) > 1e-11) goto L15;
    L += 1;
    rj1[L - 1] = x;
    if (L < nt) {
        // https://dlmf.nist.gov/10.21#E20
        x = rj1[L - 1] + (rj0[L] - rj0[L - 1]);
        goto L15;
    }

    // -- Newton method for y_{N,L}
    // initial guess for y_{N,1}
    if (n == 0) {
           x = 0.89357697;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = n + 0.93158 * pow(n, 0.33333) + 0.26035 / pow(n, 0.33333);
    }
    // iterate
    L=0;
L20:
    x0 = x;
    jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= byn/dyn;
    if (fabs(x - x0) > 1.0e-11) goto L20;
    L += 1;
    ry0[L - 1] = x;
    // initial guess for y_{N,L+1}
    if (L == 1) {
        if (n == 0) {
            x = 3.957678419314858;
        } else {
            // Expansion from https://dlmf.nist.gov/10.21#E33 and
            // coefficients from Olver 1951
            x = n + 2.59626 * pow(n, 0.33333) + 2.022183 / pow(n, 0.33333);
        }
    } else {
        // growth of roots is approximately linear (https://dlmf.nist.gov/10.21#E19)
        x = ry0[L - 1] + (ry0[L - 1] - ry0[L - 2]);
    }
    if (L <= n+10) {
        jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
        h = atan(fabs(dyn) / sqrt(fabs(fyn * byn)));
        b = -dyn / (byn * tan(h));
        x -= (h - pi/2) / b;
    }

    if (L < nt) goto L20;

    // -- Newton method for y_{N,L+1}'
    if (n == 0) {
        x = 2.67257;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = n + 1.8211 * pow(n, 0.33333) + 0.94001 / pow(n, 0.33333);
    }
    // iterate
    L=0;
L25:
    x0 = x;
    jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= dyn/fyn;
    if (fabs(x-x0) > 1.0e-11) goto L25;
    L += 1;
    ry1[L - 1] = x;
    if (L < nt) {
        // https://dlmf.nist.gov/10.21#E20
        x=ry1[L - 1] + (ry0[L] - ry0[L - 1]);
        goto L25;
    }
    return;
}
