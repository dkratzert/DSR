TITL 'NAPHTA26'@ napht-NAPHTA26.res
CELL  0.71073 8.1125 5.9398 8.649 90 124.472 90
ZERR 4 0.0007 0.0006 0.0008 0 0.003 0
SYMM  1/2-x, 1/2+y, -z
SYMM  -x, -y, -z
SYMM  1/2+x, 1/2-y, +z
SFAC C H 
UNIT 20 16
REM ######################################################
REM This file exported by ShelXle is for information or 
REM visualiztion purposes only. You may run into trouble 
REM if you try to refine it against data.
REM DISABLE_REFINE
REM ######################################################
FVAR 1.00 1 
C1    1   0.076670   0.251940  -0.075630   1.00000   0.05000
C2    1   0.013470   0.190770  -0.254470   1.00000   0.05000
C3    1  -0.082520  -0.019250  -0.328940   1.00000   0.05000
C4    1  -0.112930  -0.163550  -0.222850   1.00000   0.05000
C4A   1  -0.048140  -0.105520  -0.037400   1.00000   0.05000
C5    1  -0.076670  -0.251940   0.075630   1.00000   0.05000
C6    1  -0.013470  -0.190770   0.254470   1.00000   0.05000
C7    1   0.082520   0.019250   0.328940   1.00000   0.05000
C8    1   0.112930   0.163550   0.222850   1.00000   0.05000
C8A   1   0.048140   0.105520   0.037400   1.00000   0.05000
EQIV $1  -x, -y, -z
BIND C2 C1
BIND C3 C2
BIND C4 C3
BIND C5 C1
BIND_$1 C2 C1
BIND_$1 C3 C3
BIND_$1 C3 C2
BIND_$1 C4 C5
BIND_$1 C4 C3
BIND_$1 C5 C4
BIND_$1 C5 C1
HKLF 0 !don't refine this
END
