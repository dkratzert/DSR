
python -m cProfile -o dsr.profile -s time dsr.py -r p21c-tol.res -d

rem python -m pstats dsr.profile 