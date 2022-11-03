[![Unit tests](https://github.com/dkratzert/DSR/actions/workflows/pythonapp.yml/badge.svg)](https://github.com/dkratzert/DSR/actions/workflows/pythonapp.yml)
![Contributions](https://img.shields.io/badge/contributions-welcome-blue)

DSR
===

The program DSR consists of a text database with fragments of molecules and the DSR program. 
It acts as a preprocessor for SHELXL .res files. The user inserts a special command in the SHELXL .res file 
and the DSR program reads this information to put a molecule or fragment with the desired atoms on the position 
of the target atoms or q-peaks in the unit cell. Bond restraints can be either applied from the database to the molecule 
or automatically generated.

* [The homepage](https://dkratzert.de/dsr.html)

I apologize for the messy code. This was my first bigger project...

You have either a command line version:
```
C:\Users\daniel>dsr
----------------------------------------------------- D S R - v211 -------------------------------------
Disordered Structure Refinement (DSR)

Example DSR .res file command line:
REM DSR PUT/REPLACE "Fragment" WITH C1 C2 C3 ON Q1 Q2 Q3 PART 1 OCC -21 =
  RESI DFIX
---------------------------------------------------------------------------------------------------------
   PUT:     Just put the fragment source atoms here.
   REPLACE: Replace atoms of PART 0 in 1.3 A distance around target atoms.
---------------------------------------------------------------------------------------------------------

optional arguments:
  -h, --help            show this help message and exit
  -r "res file" ["res file" ...]
                        res file with DSR command
  -re "res file" ["res file" ...]
                        res file with DSR command (write restraints to external file)
  -e "fragment"         export fragment as .res/.png file
  -c "fragment"         export fragment to clipboard
  -t                    inverts the current fragment
  -i "tgz file"         import a fragment from GRADE (needs .tgz file)
  -l                    list names of all database entries
  -s "string"           search the database for a name
  -g                    keep group rigid (no restraints)
  -u                    Update DSR to the most current version
  -n                    do not refine after fragment transfer
```

Or a graphical user interface in [ShelXle](https://www.shelxle.org/shelx/eingabe.php):

![DSR main window](https://github.com/dkratzert/DSR/blob/master/pictures/dsr_shelxle.png?raw=true)

![DSR editor](https://github.com/dkratzert/DSR/blob/master/pictures/dsr_editor.png?raw=true)
