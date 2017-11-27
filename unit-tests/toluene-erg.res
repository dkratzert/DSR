TITL toluene
REM This file was exported by DSR version 204
REM Name: Toluene, C7H8
REM Source: CCDC CESLUJ
CELL 0.71073  11.2460  14.1233  27.1842  90.0000 100.0790  90.0000
ZERR    1.00   0.000    0.000    0.000    0.000    0.000    0.000
LATT  -1
SFAC C
UNIT 1 
REM  RESIDUE: TOL
REM Sum formula: C7 
WGHT  0.1
FVAR  1.0
rem Restraints from DSR database:
SADI C2 C3  C3 C4  C4 C5  C5 C6  C6 C7  C7 C2
SADI 0.04 C2 C6  C2 C4  C7 C5  C3 C7  C4 C6  C3 C5
DFIX 1.51 C1 C2
SADI 0.04 C1 C7  C1 C3
FLAT C1 > C7
SIMU C1 > C7
RIGU C1 > C7
rem Restraints from atom connectivities:
DFIX 1.3922 C3   C2  
DFIX 1.3775 C3   C4  
DFIX 1.5058 C2   C1  
DFIX 1.3946 C2   C7  
DFIX 1.3802 C7   C6  
DFIX 1.3814 C6   C5  
DFIX 1.3897 C5   C4  
DANG 2.5246 C1   C3  
DANG 2.5243 C1   C7  
DANG 2.4183 C2   C4  
DANG 2.4124 C2   C6  
DANG 2.3878 C3   C7  
DANG 2.3909 C3   C5  
DANG 2.3961 C4   C6  
DANG 2.3967 C5   C7  
FLAT C2 C5 C6 C7
FLAT C1 C2 C6 C7
FLAT C3 C4 C5 C6
rem end of restraints

C1   1     0.34810   0.50619   0.44851   11.0   0.04
C2   1     0.37174   0.58816   0.41613   11.0   0.04
C3   1     0.27706   0.63878   0.38821   11.0   0.04
C4   1     0.29758   0.71355   0.35825   11.0   0.04
C5   1     0.41548   0.73951   0.35559   11.0   0.04
C6   1     0.51068   0.69033   0.38312   11.0   0.04
C7   1     0.48938   0.61536   0.41297   11.0   0.04

HKLF 0
END
