TITL 'PENTAN01'@ pentane-PENTAN01-export.res
CELL  0.71073 4.1357 9.025 14.816 90 90 90
ZERR 8 0.0008 0.003 0.005 0 0 0
SYMM  1/2-x, 1/2-y, 1/2+z
SYMM  -x, +y, 1/2-z
SYMM  1/2+x, 1/2-y, -z
SYMM  -x, -y, -z
SYMM  1/2+x, 1/2+y, 1/2-z
SYMM  +x, -y, 1/2+z
SYMM  1/2-x, 1/2+y, +z
SFAC C H 
UNIT 20 48

REM ######################################################
REM This file exported by ShelXle is for information or 
REM visualiztion purposes only. You may run into trouble 
REM if you try to refine it against data.
REM DISABLE_REFINE
REM ######################################################

FVAR 1.00 1 

C1    1   0.196000   0.368860   0.087800   1.00000   0.05000
C2    1   0.094100   0.461660   0.168320   1.00000   0.05000
C3    1   0.000000   0.369650   0.250000   0.50000   0.05000
C4    1  -0.094100   0.461660   0.331680   1.00000   0.05000
C5    1  -0.196000   0.368860   0.412200   1.00000   0.05000

EQIV $1  -x, +y, 1/2-z


BIND H1 C1
BIND H2 C1
BIND H3 C1
BIND C2 C1
BIND H4 C2
BIND H5 C2
BIND C3 C2
BIND H6 C3
BIND_$1 H1 C1
BIND_$1 H2 C1
BIND_$1 H3 C1
BIND_$1 C2 C3
BIND_$1 C2 C1
BIND_$1 H4 C2
BIND_$1 H5 C2
BIND_$1 H6 C3

HKLF 0 !don't refine this
END
