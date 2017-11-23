TITL 'HEXANE01'@ hexane-HEXANE01.res
CELL  0.71073 4.1309 4.6963 8.539 83.4 87.265 75.172
ZERR 2 0.0007 0.0010 0.002 0.02 0.015 0.015
SYMM  -x, -y, -z
SFAC C H 
UNIT 6 14

REM ######################################################
REM This file exported by ShelXle is for information or 
REM visualiztion purposes only. You may run into trouble 
REM if you try to refine it against data.
REM DISABLE_REFINE
REM ######################################################

FVAR 1.00 1 

C1    1   0.223300   0.731400   0.152630   1.00000   0.05000
C2    1   0.371900   0.492880   0.283660   1.00000   0.05000
C3    1   0.427500   0.618430   0.433610   1.00000   0.05000
C4    1   0.572500   0.381570   0.566390   1.00000   0.05000
C5    1   0.628100   0.507120   0.716340   1.00000   0.05000
C6    1   0.776700   0.268600   0.847370   1.00000   0.05000


EQIV $1  1-x, 1-y, 1-z


BIND H1 C1
BIND H2 C1
BIND H3 C1
BIND C2 C1
BIND H4 C2
BIND H5 C2
BIND C3 C2
BIND H6 C3
BIND H7 C3
BIND_$1 H1 C1
BIND_$1 H2 C1
BIND_$1 H3 C1
BIND_$1 C2 C1
BIND_$1 H4 C2
BIND_$1 H5 C2
BIND_$1 C3 C3
BIND_$1 C3 C2
BIND_$1 H6 C3
BIND_$1 H7 C3

HKLF 0 !don't refine this
END
