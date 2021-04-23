Custom version of the sncf Skyrme Hartree-Fock code.


How to use:
If you'd like to use a custom EDFs, some examples of which are contained in the folder Interactions, you can proceed in two ways.


1. copy the file containing the EDF (e.g. av4p_256.in) to interaction.in. Then run:
   ./sncf.x [A] [Z] [c0] [c1] [w0]
   All arguments are optional. c0, c1, wo are the gradient and spin-orbit coefficients, which by default are set to zero.

   E.g. cp Interactions/av4p_256.in interaction.in
        ./sncf.x 90 40 -150 20 140
     
   To read from sncf.in instead, it is enough to modify interaction.in and set the number of terms to a negative number.


2. use my_script.sh. 
   This script computes several different nuclei in one shot and outputs many files saved in the Results folder. Useful are the "tab_" files and the densities.
   The syntax is:
   ./my_script.sh [name interaction] [c0] [c1] [w0]

    E.g. ./my_script.sh nnlosat_23456.in -25 10 50
