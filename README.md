# AlignAndCompareMolecules

AlignAndCompareMolecules is a program that compares and aligns two molecules, which have the same number and type of atoms, and recognizes if the molecules are Equal, Different or Enantiomers between each other.

### Extra programs:
A program to rotate the coordinates of a molecule.

# Git instructions

You can obtain the source code of alignandcomparemolecules as follows.
In your bash terminal type:

~~~~~~~~~~
$cd /local/path/to/alignandcomparemolecules
$git clone https://github.com/LuisAlfredoNu/alignandcomparemolecules.git
~~~~~~~~~~

After this, git will transfer the source files to ```/local/path/to/alignandcomparemolecules```

# Building and installing the package

## Auxiliary programs/dependencies
The following dependencies are needed:

~~~~~~~~~~
g++ std=c++11 make
~~~~~~~~~~


## Compilation

For building the binaries, type:

~~~~~~~~~~
$cd /local/path/to/alignandcomparemolecules/src/
$make
~~~~~~~~~~

## Installation
For installing the program, type:

~~~~~~~~~~
$cd /local/path/to/alignandcomparemolecules/src/
$sudo make install
~~~~~~~~~~

# Testing the suite

How to test the program (or how it should be):

~~~~~~~~~~
$make runQuickTest
~~~~~~~~~~

## Big test
If you want a huge test of the program, which produces a pfd report, you need the following dependencies to be installed in your system:

~~~~~~~~~~
vmd tachyon_LINUXAMD64 pdflatex
~~~~~~~~~~

The library tachyon_LINUXAMD64 is included the installation directory of vmd. You can check whether this file exist in the following path: 

~~~~~~~~~~
/usr/local/lib/vmd/tachyon_LINUXAMD64
~~~~~~~~~~

If the library exists and work OK, the following command should produce a png file:

~~~~~~~~~~
/local/path/to/alignandcomparemolecules/bigTest/Scripts/getPNGofXYZ
~~~~~~~~~~

If the above command works, then run the script 

~~~~~~~~~~
$cd alignandcomparemolecules/bigTest
$./scriptBig_Test.bsh
~~~~~~~~~~

The command will generate a report that will be placed in the following folder 

~~~~~~~~~~
/local/path/to/alignandcomparemolecules/bigTest/tex/report_of_comparisons.pdf
~~~~~~~~~~

# Updating the program (git instructions)

If you only need to get the latest version of alignandcomparemolecules, and you are not using personal modifications to the source, in the local main directory (using the above example should be ```/local/path/to/alignandcomparemolecules```) type:

~~~~~~~~~~
$git pull
$cd src
$make update
$sudo make install
~~~~~~~~~~


# User manual
To get rough instructions of how to use the program, type:

~~~~~~~~~~
$alignandcomparemolecules -h
~~~~~~~~~~

#Developer manual

Under construction ...

#Contributions and List of Contributors






