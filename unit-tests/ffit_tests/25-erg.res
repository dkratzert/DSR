TITL p-1_a.res in P-1
    p-1_a.res
    created by SHELXL-2018/3 at 14:30:49 on 26-Nov-2018
REM Old TITL p-1 in P-1
REM SHELXT solution in P-1: R1 0.134, Rweak 0.002, Alpha 0.043
REM <I/s>   0.000 for    0 systematic absences,  Orientation as input
REM Formula found by SHELXT:  C30 B3 F11 N16 Pt2
CELL  0.71073   9.2011  16.2534  16.2490   82.839   73.563   73.564
ZERR   22.000   0.0031   0.0055   0.0054    0.009    0.009    0.008
LATT  1
SFAC  C  H  B  F  N  PT  CL
UNIT  1  2  3  4  5  6  7
TEMP -33.150
L.S. 10
XNPD
RIGU
DELU
SIMU 0.04 0.08 1
BOND
LIST 4
ACTA
FMAP 2
PLAN 60
REM REM DSR REPLACE ACETONITRILE WITH N1 C1 C2 ON N00C C00R C010 OCC 11
REM !RESI MECN
REM REM DSR REPLACE MOFUG WITH CL1 C1 CL2 ON Q1 Q3 N01I OCC 61 RESI CCL2
REM Restraints for Fragment mofug, Dichloromethane, CH2Cl2, DCM SADI from:
REM . Please cite https://doi.org/10.1107/S1600576718004508
SADI_CCL2 CL1 C1 CL2 C1
SADI_CCL2 0.03 CL1 CL2
SIMU_CCL2 CL1 > C1
RIGU_CCL2 CL1 > C1
SAME_CCL2 CL1 > C1
REM REM DSR PUT BF4 WITH F2 B1 F3 ON F4_3A B1_3A Q40 PART 2 OCC -31 RESI
REM !BF4
REM REM DSR PUT BF4 WITH F1 B1 F3 ON Q10 Q54 F2_6A PART 2 OCC -51 RESI BF4
WGHT    0.100000
FVAR 0.19261 0.87434 0.53257 0.83155 0.47353 0.32938
REM REM DSR REPLACE CH2CL2 WITH CL1 C1 CL2 ON C01P C01M Q1 OCC 11 RESI CCL2

REM Restraints for Fragment ch2cl2, Dichloromethane, CH2Cl2, DCM from:
REM CCDC ASOMAX. Please cite https://doi.org/10.1107/S1600576718004508
DFIX_CCL2 1.771 0.01 CL1 C1 CL2 C1
DFIX_CCL2 2.916 0.03 CL1 CL2
RESI CCL2 14
CL1   7    0.180364    0.451428   0.648831   11.00000    0.03000
CL2   7    0.290421    0.266156   0.677484   11.00000    0.03000
C1    1    0.284665    0.352415   0.599805   11.00000    0.03000
RESI 0

RESI BF4 13
PART 2  -51
B1    3    0.482086    0.828746    0.646855   -51.00000    0.03770    0.09892 =
         0.05336    0.01913   -0.00965   -0.01308
F1    4    0.551745    0.789014    0.570774   -51.00000    0.14056    0.13951 =
         0.08504   -0.01012    0.00043   -0.04110
F2    4    0.557658    0.803815    0.706951   -51.00000    0.16651    0.24846 =
         0.13016    0.02173   -0.10330   -0.05096
F3    4    0.442127    0.913378    0.633478   -51.00000    0.10845    0.09592 =
         0.13749    0.01441    0.00768   -0.03958
F4    4    0.348708    0.801271    0.681566   -51.00000    0.08915    0.12882 =
         0.17996    0.02220   -0.00126   -0.05658
PART 0
RESI 0
RESI BF4 12
PART 2  -31
B1    3   -0.023136    0.647682    0.828070   -31.00000    0.03616    0.04481 =
         0.07570    0.03084   -0.01495   -0.02127
F1    4   -0.055220    0.627091    0.912095   -31.00000    0.09578    0.13698 =
         0.06919    0.04168   -0.03196   -0.01298
F2    4   -0.163371    0.684430    0.804191   -31.00000    0.06834    0.10821 =
         0.10029    0.04261   -0.05034   -0.03856
F3    4    0.066861    0.701215    0.806573   -31.00000    0.08937    0.14199 =
         0.17502    0.01770   -0.02714   -0.09614
F4    4    0.052054    0.572480    0.787751   -31.00000    0.15148    0.08152 =
         0.11577   -0.01063   -0.04124    0.00159
PART 0
RESI 0
RESI CCL2 11
CL1   7    0.775865    0.891301    1.139414    61.00000    0.14397    0.29513 =
         0.20560   -0.09092   -0.05851    0.07384
CL2   7    0.804862    0.789345    1.021201    61.00000    0.10959    0.12586 =
         0.15134    0.01616   -0.04564   -0.01256
C1    1    0.712253    0.849100    1.089783    61.00000    0.06146    0.16700 =
         0.12097    0.02861    0.01525   -0.03706
RESI 0
RESI MECN 10
RESI 0
REM REM DSR REPLACE ACETONITRILE WITH N1 C1 C2 ON N008 C00P C00V OCC 11
REM !RESI MECN
RESI MECN 9
RESI 0
REM REM DSR REPLACE ACETONITRILE WITH N1 C1 C2 ON N00B C00S C00U OCC 11
REM !RESI MECN
RESI MECN 8
RESI 0
REM REM DSR REPLACE ACETONITRILE WITH N1 C1 C2 ON N009 C00O C00W OCC 11
REM !RESI MECN
REM Restraints for Fragment acetonitrile, Acetonitrile, NMe, C2H3N from:
REM CCDC CISDOZ. Please cite https://doi.org/10.1107/S1600576718004508
DFIX_MECN 1.137 0.01 C1 N1
DFIX_MECN 1.456 0.01 C1 C2
DFIX_MECN 2.593 0.02 C2 N1
SIMU_MECN N1 > C2
RIGU_MECN N1 > C2
SAME_MECN N1 > C2
RESI MECN 7
RESI 0
REM REM DSR REPLACE BF4 WITH F3 B1 F2 ON F007 C01D N011 PART 1 OCC 51 RESI
REM !BF4
RESI BF4 6
PART 1  51
B1    3    0.453510    0.851638    0.671249    51.00000    0.04160    0.04494 =
         0.07503    0.02380   -0.02674   -0.02031
F1    4    0.484822    0.873820    0.588114    51.00000    0.14004    0.12561 =
         0.07843    0.04077   -0.02411   -0.08256
F2    4    0.412658    0.927052    0.711343    51.00000    0.14151    0.08269 =
         0.12699   -0.01652   -0.02215   -0.03477
F3    4    0.321579    0.818454    0.694106    51.00000    0.05077    0.12067 =
         0.10393    0.03094   -0.02397   -0.04684
F4    4    0.569692    0.795092    0.695562    51.00000    0.06116    0.10828 =
         0.19650    0.02722   -0.07287   -0.00361
PART 0
RESI 0
REM REM DSR REPLACE BF4 WITH F3 B1 F2 ON F00M B01O F018 PART 1 OCC 41 RESI
REM !BF4
RESI BF4 5
PART 1  41
B1    3   -0.246503    1.097761    0.904167    41.00000    0.08635    0.06830 =
         0.06524   -0.00048    0.00864    0.00422
F1    4   -0.180438    1.158938    0.890257    41.00000    0.13540    0.10799 =
         0.19799   -0.03516    0.04492   -0.03641
F2    4   -0.321661    1.107837    0.841987    41.00000    0.14815    0.21236 =
         0.10744   -0.03605   -0.04445    0.04689
F3    4   -0.150198    1.016192    0.894811    41.00000    0.10853    0.08441 =
         0.10745    0.00057   -0.00034    0.02951
F4    4   -0.349688    1.106105    0.984254    41.00000    0.11052    0.11163 =
         0.08260    0.00387    0.02730    0.00641
PART 0
RESI 0
REM REM DSR REPLACE BF4 WITH F3 B1 F2 ON N00J B01N C019 PART 1 OCC 21 RESI
REM !BF4
RESI BF4 4
PART 1  21
B1    3    0.252923    0.597456    0.404114    21.00000    0.14716    0.07513 =
         0.07928    0.00055   -0.05332   -0.05062
F1    4    0.259616    0.606043    0.484194    21.00000    0.22438    0.12524 =
         0.08567    0.00691   -0.07482   -0.07986
F2    4    0.370279    0.609754    0.340806    21.00000    0.23238    0.24394 =
         0.11168   -0.02451   -0.00863   -0.16549
F3    4    0.241352    0.515977    0.394652    21.00000    0.21009    0.07996 =
         0.12285   -0.00300   -0.06414   -0.06916
F4    4    0.129616    0.659442    0.391243    21.00000    0.23333    0.10376 =
         0.23387   -0.03654   -0.15630    0.00359
PART 0
RESI 0
REM REM DSR REPLACE BF4 WITH F3 B1 F2 ON C01Q B01J F00H PART 1 OCC 21 RESI
REM !BF4
REM Restraints for Fragment bf4, Tetrafluoroborate, [BF4]- from: TURBOMOLE
REM V6.0 B-P86/def-SV(P). Please cite
REM https://doi.org/10.1107/S1600576718004508
SIMU_BF4 B1 > F4
RIGU_BF4 B1 > F4
SADI_BF4 0.02 B1 F1 B1 F2 B1 F3 B1 F4
SADI_BF4 0.04 F1 F2 F2 F3 F3 F4 F4 F1 F2 F4 F1 F3
SAME_BF4 B1 > F4
RESI BF4 3
PART 1  31
B1    3   -0.042902    0.671295    0.852239    31.00000    0.07419    0.09467 =
         0.06613    0.01471   -0.03607   -0.04383
F1    4   -0.086757    0.711852    0.928340    31.00000    0.15903    0.14937 =
         0.08466   -0.01062   -0.03276   -0.02338
F2    4    0.061971    0.699777    0.790291    31.00000    0.11452    0.19317 =
         0.11491    0.02799    0.00645   -0.05245
F3    4   -0.003608    0.586428    0.864870    31.00000    0.17093    0.09908 =
         0.14935    0.02579   -0.09662   -0.02053
F4    4   -0.172143    0.694411    0.819337    31.00000    0.11016    0.17650 =
         0.16623    0.01968   -0.09326   -0.01246
PART 0
RESI 0
REM REM DSR REPLACE BIPYRIDINE WITH N1 N2 C6 ON F006 F003 N00E OCC 11 RESI
REM !BIP
RESI BIP 2
RESI 0
REM REM DSR REPLACE BIPYRIDINE WITH N1 N2 C6 ON F004 F005 N00D OCC 11 RESI
REM !BIP
REM Restraints for Fragment bipyridine, 2,2'-Bipyridine, bipy, C10H8N2
REM from: CCDC DOCYEB. Please cite
REM https://doi.org/10.1107/S1600576718004508
SADI_BIP N1 C1 N1 C5 N2 C6 N2 C10
SADI_BIP C2 C3 C3 C4 C4 C5 C6 C7 C7 C8 C8 C9 C9 C10
SADI_BIP 0.02 C5 C6
SADI_BIP 0.04 N1 C4 N1 C2 N2 C7 N2 C9
SADI_BIP 0.04 C1 C5 C1 C3 C2 C4 C3 C5 C6 C10 C6 C8 C7 C9
FLAT_BIP N1 > C10
SIMU_BIP N1 > C10
RIGU_BIP N1 > C10
SAME_BIP N1 > C10
RESI BIP 1
RESI 0
N01G  5    0.300500    1.024120    0.791998    11.00000    0.11330    0.15216 =
         0.13933    0.03445    0.00364   -0.03176
C01M  1    0.382976    0.291457    0.521379    11.00000    0.11345    0.10426 =
         0.11505    0.02453   -0.04422   -0.05113
HKLF 4




REM  p-1_a.res in P-1
REM wR2 = 0.1183, GooF = S = 1.013, Restrained GooF = 1.099 for all data
REM R1 = 0.0335 for 8839 Fo > 4sig(Fo) and 0.0379 for all 10249 data
REM 681 parameters refined using 2247 restraints

END 
 
WGHT      0.0841      0.7053 

REM Highest difference peak  2.490,  deepest hole -0.702,  1-sigma level  0.159
Q1    1   0.2512  0.3015  0.6978  11.00000  0.05    2.49
Q2    1   0.2517  0.3260  0.6692  11.00000  0.05    2.38
Q3    1   0.7474  0.7971  1.2008  11.00000  0.05    2.37
Q4    1   0.3149  0.3385  0.5849  11.00000  0.05    1.86
Q5    1   0.7658  0.8413  1.0879  11.00000  0.05    1.82
Q6    1   0.1628  0.4233  0.6207  11.00000  0.05    1.62
Q7    1   0.1776  0.4158  0.6645  11.00000  0.05    1.57
Q8    1   0.3285  0.3818  0.5741  11.00000  0.05    1.55
Q9    1   0.7244  0.8873  1.0720  11.00000  0.05    1.43
Q10   1   0.7383  0.9141  1.1658  11.00000  0.05    1.38
Q11   1  -0.2562  1.0818  0.9181  11.00000  0.05    1.08
Q12   1  -0.2523  1.1532  0.8468  11.00000  0.05    0.78
Q13   1   0.6060  0.8408  0.6599  11.00000  0.05    0.73
Q14   1  -0.2412  1.0238  0.8886  11.00000  0.05    0.72
Q15   1  -0.0751  0.8836  0.6144  11.00000  0.05    0.72
Q16   1   0.4342  0.6101  0.8808  11.00000  0.05    0.68
Q17   1  -0.0941  1.0960  0.9136  11.00000  0.05    0.68
Q18   1  -0.2474  1.1183  0.9705  11.00000  0.05    0.67
Q19   1   0.2495  0.6491  0.3444  11.00000  0.05    0.66
Q20   1  -0.4162  1.0897  0.9120  11.00000  0.05    0.65
Q21   1  -0.0049  0.9244  0.5747  11.00000  0.05    0.64
Q22   1   0.0855  0.5939  0.4154  11.00000  0.05    0.63
Q23   1   0.4938  0.5695  0.9212  11.00000  0.05    0.62
Q24   1   0.3947  0.6237  0.6140  11.00000  0.05    0.61
Q25   1   0.1073  0.6519  0.8394  11.00000  0.05    0.61
Q26   1  -0.2980  1.0925  1.0024  11.00000  0.05    0.59
Q27   1   0.2232  0.5955  0.4952  11.00000  0.05    0.59
Q28   1   0.0537  0.9096  0.5871  11.00000  0.05    0.58
Q29   1   0.3427  0.5291  0.3792  11.00000  0.05    0.57
Q30   1   0.1481  0.6233  0.4717  11.00000  0.05    0.56
Q31   1   0.4156  0.5807  0.4110  11.00000  0.05    0.53
Q32   1  -0.1036  0.6197  0.6167  11.00000  0.05    0.51
Q33   1  -0.0602  0.9124  0.6560  11.00000  0.05    0.50
Q34   1  -0.0367  0.5885  0.8927  11.00000  0.05    0.46
Q35   1   0.0120  0.8324  0.6566  11.00000  0.05    0.46
Q36   1   0.4240  0.6407  0.9287  11.00000  0.05    0.43
Q37   1   0.4504  0.8368  0.5911  11.00000  0.05    0.43
Q38   1  -0.3788  0.8775  0.8884  11.00000  0.05    0.42
Q39   1   0.0810  0.8481  0.5869  11.00000  0.05    0.42
Q40   1   0.2994  0.8021  0.6972  11.00000  0.05    0.41
Q41   1   0.5217  0.6621  0.8406  11.00000  0.05    0.41
Q42   1  -0.1317  0.6547  0.9076  11.00000  0.05    0.39
Q43   1   0.0552  0.8690  0.6543  11.00000  0.05    0.39
Q44   1   0.3224  0.3089  0.5046  11.00000  0.05    0.38
Q45   1   0.5780  0.7863  0.7171  11.00000  0.05    0.37
Q46   1   0.2455  0.9622  0.4642  11.00000  0.05    0.36
Q47   1   0.1489  0.5027  0.7027  11.00000  0.05    0.35
Q48   1   0.6788  0.5410  1.0391  11.00000  0.05    0.35
Q49   1  -0.3862  1.1266  0.9722  11.00000  0.05    0.35
Q50   1   0.3772  0.9269  0.4252  11.00000  0.05    0.34
Q51   1   0.0926  0.8642  0.6444  11.00000  0.05    0.34
Q52   1   0.6077  0.6527  1.1266  11.00000  0.05    0.34
Q53   1   0.1821  1.0470  0.5420  11.00000  0.05    0.32
Q54   1   0.5265  0.6127  0.8213  11.00000  0.05    0.32
Q55   1   0.0867  0.7141  0.7826  11.00000  0.05    0.32
Q56   1   0.4624  0.9067  0.5997  11.00000  0.05    0.32
Q57   1   0.5357  0.8988  0.6085  11.00000  0.05    0.32
Q58   1   0.4096  0.8541  0.3716  11.00000  0.05    0.32
Q59   1  -0.1511  0.6381  0.8235  11.00000  0.05    0.31
Q60   1  -0.2018  0.6982  0.8107  11.00000  0.05    0.31
