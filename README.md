Copyright 2020 Joshua Laniado and Todd O. Yeates.


NANOHEDRA SETUP
1) COMPILE 'orient_oligomer.f'
in the nanohedra directory execute the following commands:
cd orient
gfortran -o orient_oligomer orient_oligomer.f

2) INSTALL 'FreeSASA 2.0.3'
go to: https://freesasa.github.io for a quick installation guide

3) INSTALL the following Python modules that support Python 2.7:
'sklearn' version 0.20.1 (https://scikit-learn.org/0.20)
'biopython' version 1.72 (https://biopython.org/wiki/Download)


RUNNING NANOHEDRA
to access the user manual page: 
python nanohedra.py 


NOTES
nanohedra currently only supports Python 2.7
nanohedra currently runs on Linux/Unix and Mac

