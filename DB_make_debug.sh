#!/bin/bash -l
#Determines the Fortran compiler to use for DigiBog based on where the user is logged on
#Dylan Young, 15th May 2017.
#Version 150517
#Run this file by typing 'source filename.sh' at the bash prompt
# ==============================================================================

echo 'Setting Fortan compiler for debugging with gdb .....'

comp='gfortran'

#Export 'comp' so it can be used by the makefile
export comp
#Confirm which compiler has been selected
echo "Fortran compiler set to ${comp}"
#Remove old files created by the compiler and compile a new version
make -f Make_DB_debug clean && make -f Make_DB_debug
#Clear the compiler
unset comp
