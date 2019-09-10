# Welcome to arensi's project 1 for FYS3150 


In this project we solve the one-dimensional Poisson equation with Dirichlet bound-ary conditions by rewriting it as a set of linear equations. We use different algorithms found in this repository to find approximated values to a specific given matrix.

To compile most of the programs follow the following instructions:

Compile with:   g++ -02 -std=c++0x [C++ program] -o [Chosen name of program] -larmadillo

and then run with: 		./[Chosen name of program] [Output file name] [Wanted exponent of n, where n is the dimension of the matrix]

For example: 

Compile with:  g++ -02 -std=c++0x 1c.cpp -o test_c -larmadillo

and then run with: 		./test_c output_c 3


This formula will work with all of the programs in this repository

The output test files are provided to demonstrate the results of the programs. 
The index was removed for graphing simplicity, but are in order: X   Approximaton   Exact   Relative Error
