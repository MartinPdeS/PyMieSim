#pragma once

#include <cmath>

/**
 * @brief Factorial function using recursion.
 * @param n Non-negative integer.
 * @return Factorial of n.
 */
unsigned long long factorial(int n);


/**
 * @brief Computes the Bessel function of the first kind Jn(x) using series expansion.
 * @param n Order of the Bessel function.
 * @param x Argument of the Bessel function.
 * @return Value of Jn(x).
 */
double bessel_J(int n, double x);

/**
 * @brief Determines the starting point for backward recurrence such that the magnitude of Jn(x) at that point is about 10^(-MP).
 * @param x Argument of Jn(x).
 * @param mp Value of magnitude.
 * @return Starting point for backward recurrence.
 */
int msta1(double x, int mp);

/**
 * @brief Determines the starting point for backward recurrence such that all Jn(x) has MP significant digits.
 * @param x Argument of Jn(x).
 * @param n Order of Jn(x).
 * @param mp Significant digits.
 * @return Starting point for backward recurrence.
 */
int msta2(double x, int n, int mp);

/**
 * @brief Computes Bessel functions jn(x) and yn(x) using backward recurrence.
 * @param n Order of the Bessel functions.
 * @param nmin Minimum order to compute.
 * @param x Argument of the Bessel functions (x > 0).
 * @param nm Pointer to store the number of computed orders.
 * @param bj Pointer to store jn(x).
 * @param by Pointer to store yn(x).
 */
template <typename T> void jynbh(int n, int nmin, T x, int *nm, T *bj, T *by);

/**
 * @brief Computes Bessel functions jn(x) and yn(x), and their first and second derivatives.
 * @param n Order of the Bessel functions.
 * @param x Argument of the Bessel functions (x > 0).
 * @param bjn Pointer to store jn(x).
 * @param djn Pointer to store jn'(x).
 * @param fjn Pointer to store jn"(x).
 * @param byn Pointer to store yn(x).
 * @param dyn Pointer to store yn'(x).
 * @param fyn Pointer to store yn"(x).
 */
void jyndd(int n, double x, double *bjn, double *djn, double *fjn, double *byn, double *dyn, double *fyn);

/**
 * @brief Computes the zeros of Bessel functions Jn(x), Yn(x), and their derivatives.
 * @param n Order of the Bessel functions (n >= 0).
 * @param nt Number of zeros (roots).
 * @param rj0 Pointer to store the zeros of Jn(x).
 * @param rj1 Pointer to store the zeros of Jn'(x).
 * @param ry0 Pointer to store the zeros of Yn(x).
 * @param ry1 Pointer to store the zeros of Yn'(x).
 */
void bessel_zeros(int n, int nt, double *rj0, double *rj1, double *ry0, double *ry1);