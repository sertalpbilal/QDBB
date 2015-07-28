NAME
* Written by MOSEK version 7.1.0.29
OBJSENSE
    MIN
ROWS
 N  OBJ
 G  C1      
 L  C2      
 L  C3      
 L  C4      
 L  C5      
 L  C6      
 L  C7      
 L  C8      
 L  C9      
 L  C10     
 L  C11     
 E  C12     
 L  C13     
 L  C14     
 G  C15     
 G  C16     
COLUMNS
    X1        OBJ       1e+0
    X1        C13       -1e+0
    X2        OBJ       0
    X2        C1        4.7e-2
    X2        C2        1e+0
    X2        C12       1e+0
    X3        OBJ       0
    X3        C1        9e-3
    X3        C3        1e+0
    X3        C12       1e+0
    X4        OBJ       0
    X4        C1        -1.4e-2
    X4        C4        1e+0
    X4        C12       1e+0
    X5        OBJ       0
    X5        C1        4.8e-2
    X5        C5        1e+0
    X5        C12       1e+0
    X6        OBJ       0
    X6        C1        -1e-3
    X6        C6        1e+0
    X6        C12       1e+0
    X7        OBJ       0
    X7        C1        2e-3
    X7        C7        1e+0
    X7        C12       1e+0
    X8        OBJ       0
    X8        C1        7.3e-2
    X8        C8        1e+0
    X8        C12       1e+0
    X9        OBJ       0
    X9        C1        2.1e-2
    X9        C9        1e+0
    X9        C12       1e+0
    X10       OBJ       0
    X10       C1        1.3e-2
    X10       C10       1e+0
    X10       C12       1e+0
    X11       OBJ       0
    X11       C1        8.6e-2
    X11       C11       1e+0
    X11       C12       1e+0
    0         'MARKER'                 'INTORG'
    X12       OBJ       0
    X12       C2        -1e+0
    X13       OBJ       0
    X13       C3        -1e+0
    X14       OBJ       0
    X14       C4        -1e+0
    X15       OBJ       0
    X15       C5        -1e+0
    X16       OBJ       0
    X16       C6        -1e+0
    X17       OBJ       0
    X17       C7        -1e+0
    X18       OBJ       0
    X18       C8        -1e+0
    X19       OBJ       0
    X19       C9        -1e+0
    X19       C16       -2.67949192e-1
    X20       OBJ       0
    X20       C10       -1e+0
    X21       OBJ       0
    X21       C11       -1e+0
    X21       C15       -2.67949192e-1
    1         'MARKER'                 'INTEND'
    X22       OBJ       0
    X22       C15       1e+0
    X23       OBJ       0
    X23       C16       1e+0
RHS
    rhs       C1        6e-2
    rhs       C12       1e+0
    rhs       C14       4e+0
    rhs       C15       -2e+0
    rhs       C16       -2e+0
RANGES
    ran       C1        0
    ran       C15       0
    ran       C16       0
BOUNDS
 UP bound     X2        1e+0
 UP bound     X3        1e+0
 UP bound     X4        1e+0
 UP bound     X5        1e+0
 UP bound     X6        1e+0
 UP bound     X7        1e+0
 UP bound     X8        1e+0
 UP bound     X9        1e+0
 UP bound     X10       1e+0
 UP bound     X11       1e+0
 BV bound     X12     
 BV bound     X13     
 BV bound     X14     
 BV bound     X15     
 BV bound     X16     
 BV bound     X17     
 BV bound     X18     
 BV bound     X19     
 BV bound     X20     
 BV bound     X21     
QSECTION      C13     
    X2        X2        2.24302624e+5
    X3        X2        3.7976526e+3
    X6        X2        3.11950035e+3
    X4        X2        1.9612498e+3
    X11       X2        1.71622954e+5
    X10       X2        1.01331016e+3
    X5        X2        1.71649075e+5
    X7        X2        2.75058553e+3
    X8        X2        7.85188452e+4
    X9        X2        2.8298205e+3
    X3        X3        8.20125e+2
    X4        X3        2.1346578e+3
    X6        X3        3.772575e+2
    X5        X3        3.3932034e+3
    X11       X3        4.941729e+3
    X10       X3        1.28672145e+3
    X7        X3        1.82953485e+3
    X8        X3        1.924803e+3
    X9        X3        8.2134e+2
    X6        X4        8.8944075e+2
    X4        X4        4.28717762e+4
    X5        X4        4.32940226e+3
    X7        X4        3.68773994e+4
    X8        X4        4.63885444e+3
    X9        X4        3.34034415e+3
    X10       X4        2.74664574e+4
    X11       X4        -1.78646554e+3
    X6        X5        7.1856072e+3
    X5        X5        4.85782531e+5
    X7        X5        1.07943788e+4
    X8        X5        1.31167198e+5
    X9        X5        4.5809478e+3
    X10       X5        1.49123527e+4
    X11       X5        3.00676756e+5
    X7        X6        8.3160675e+2
    X10       X6        1.16417655e+3
    X11       X6        5.18881545e+3
    X6        X6        8.20125e+2
    X8        X6        2.9513646e+3
    X9        X6        5.3044875e+2
    X7        X7        3.74777442e+4
    X8        X7        3.46977821e+3
    X10       X7        2.56805092e+4
    X9        X7        2.42911305e+3
    X11       X7        6.68121761e+3
    X9        X8        9.370543e+3
    X8        X8        2.00775171e+5
    X10       X8        1.72565005e+4
    X11       X8        1.43042678e+5
    X9        X9        3.570125e+3
    X10       X9        5.75280225e+3
    X11       X9        4.63973445e+3
    X10       X10       4.57773282e+4
    X11       X10       7.38404129e+3
    X11       X11       7.44419616e+5
QSECTION      C14     
    X12       X12       2e+0
    X13       X13       2e+0
    X14       X14       2e+0
    X15       X15       2e+0
    X16       X16       2e+0
    X17       X17       2e+0
    X18       X18       2e+0
    X19       X19       2e+0
    X20       X20       2e+0
    X21       X21       2e+0
CSECTION      K1        0              QUAD
    X22     
    X12     
    X13     
    X14     
    X15     
    X16     
    X17     
    X18     
    X19     
    X20     
ENDATA
