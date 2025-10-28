#define FE_1(M, X1)                             M(X1)
#define FE_2(M, X1, X2)                         M(X1) M(X2)
#define FE_3(M, X1, X2, X3)                     M(X1) M(X2) M(X3)
#define FE_4(M, X1, X2, X3, X4)                 M(X1) M(X2) M(X3) M(X4)
#define FE_5(M, X1, X2, X3, X4, X5)             M(X1) M(X2) M(X3) M(X4) M(X5)
#define FE_6(M, X1, X2, X3, X4, X5, X6)         M(X1) M(X2) M(X3) M(X4) M(X5) M(X6)
#define FE_7(M, X1, X2, X3, X4, X5, X6, X7)     M(X1) M(X2) M(X3) M(X4) M(X5) M(X6) M(X7)
#define FE_8(M, X1, X2, X3, X4, X5, X6, X7, X8) M(X1) M(X2) M(X3) M(X4) M(X5) M(X6) M(X7) M(X8)

#define GET_FE(_1,_2,_3,_4,_5,_6,_7,_8, NAME, ...)     NAME

#define FOR_EACH(FUNCTION, ...)         GET_FE(__VA_ARGS__, FE_8, FE_7, FE_6, FE_5, FE_4, FE_3, FE_2, FE_1)(FUNCTION, __VA_ARGS__)
