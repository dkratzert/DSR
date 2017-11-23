TITL 'OFIZAG'@ diox-OFIZAG.res
CELL  0.71073 6.847 7.0982 8.6412 106.676 105.448 90.679
ZERR 2 0.0004 0.0004 0.0003 0.005 0.005 0.004
SYMM  -x, -y, -z
SFAC C H N O S 
UNIT 10 24 4 2 2
REM ######################################################
REM This file exported by ShelXle is for information or 
REM visualiztion purposes only. You may run into trouble 
REM if you try to refine it against data.
REM DISABLE_REFINE
REM ######################################################
FVAR 1.00 1 
O1    4   0.190880   0.039150   0.622930   1.00000   0.05000
O2    4  -0.190880  -0.039150   0.377070   1.00000   0.05000
C1    1   0.009080  -0.004220   0.662560   1.00000   0.05000
C2    1  -0.172280   0.057560   0.551210   1.00000   0.05000
C3    1  -0.009080   0.004220   0.337440   1.00000   0.05000
C4    1   0.172280  -0.057560   0.448790   1.00000   0.05000



EQIV $1  -x, -y, 1-z
BIND C4 O1
BIND C5 O1
BIND_$1 C4 C5
BIND_$1 C4 O1
BIND_$1 C5 C4
BIND_$1 C5 O1
HKLF 0 !don't refine this
END
