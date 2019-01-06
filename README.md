# AlignAndCompareMolecules

Program to compare and align two molecules with the same number and type of atoms and recognize if the moleculas are Equal, Different or Enantiomer.

### Extra:
Program to rotate the coordinates of a molecule.

# Git instructions

You can get alignandcomparemolecules source code as follows.
In your bash terminal type:

~~~~~~~~~~
$cd /local/path/to/alignandcomparemolecules
$git clone https://github.com/LuisAlfredoNu/alignandcomparemolecules.git
~~~~~~~~~~

After this, git will transfer the source files to ```/local/path/to/alignandcomparemolecules```

# Building the package

## Auxiliary programs/dependencies
The following dependencies are needed:

~~~~~~~~~~
g++ std=c++11 make
~~~~~~~~~~

## Compilation

For building the binaries, type:

~~~~~~~~~~
$cd alignandcomparemolecules/src/
$make
~~~~~~~~~~

For installing, type:

~~~~~~~~~~
$cd alignandcomparemolecules/src/
$sudo make install
~~~~~~~~~~

# Testing the suite

How to test the program (should be):

~~~~~~~~~~
$make runQuickTest
~~~~~~~~~~

## Big test
If you want a huge test of the program with a report on PDF format, you need the following dependencies

~~~~~~~~~~
vmd tachyon_LINUXAMD64 pdflatex
~~~~~~~~~~

The library tachyon_LINUXAMD64 are in the installation folder of vmd, if the file exist in the path 

~~~~~~~~~~
/usr/local/lib/vmd/tachyon_LINUXAMD64
~~~~~~~~~~

Do nothing, if not update the file 

~~~~~~~~~~
alignandcomparemolecules/bigTest/Scripts/getPNGofXYZ
~~~~~~~~~~

And run the script 

~~~~~~~~~~
$alignandcomparemolecules/bigTest/scriptBig_Test.bsh
~~~~~~~~~~

The report are in the following folder 

~~~~~~~~~~
alignandcomparemolecules/bigTest/tex/report_of_comparations.pdf
~~~~~~~~~~
# Updating the program (git instructions)

If you only need to get the latest version of demdysigma2d, and you are not using personal modifications to the source, in the local main directory (using the above example should be ```/local/path/to/alignandcomparemolecules```):

   ~~~~~~~~~~
   $ git pull
   $ cd src
   $ make update
   $ sudo make install
   ~~~~~~~~~~


# Is there a manual?

#Developer instructions

Under construction

#Contributions and List of Contributors






