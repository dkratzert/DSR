TITL 'HIXHIF05'@ acetone-HIXHIF05.res
CELL  0.71073 6.3925 5.3424 10.733 90 90 90
ZERR 8 0.0017 0.0014 0.003 0 0 0
SYMM  -x, -y, 1/2+z
SYMM  -x, 1/2+y, 1/2-z
SYMM  +x, 1/2-y, -z
SYMM  -x, -y, -z
SYMM  +x, +y, 1/2-z
SYMM  +x, 1/2-y, 1/2+z
SYMM  -x, 1/2+y, +z
SFAC C H O 
UNIT 12 24 4
REM ######################################################
REM This file exported by ShelXle is for information or 
REM visualiztion purposes only. You may run into trouble 
REM if you try to refine it against data.
REM DISABLE_REFINE
REM ######################################################
FVAR 1.00 1 
O1    3   0.296100   0.183000   0.250000   0.50000   0.05000
C1    1   0.259300   0.404400   0.250000   0.50000   0.05000
C2    1   0.232000   0.550600   0.131290   1.00000   0.05000
C3    1   0.232000   0.550600   0.368710   1.00000   0.05000

EQIV $1  +x, +y, 1/2-z

BIND C1 O1
BIND C2 C1
BIND_$1 C2 C1
HKLF 0 !don't refine this
END
