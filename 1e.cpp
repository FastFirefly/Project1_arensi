#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include "time.h"
using namespace std;
using namespace arma;



inline double f(double x){
	return 100.0*exp(-10.0*x);
}
inline double u(double x) {
	return 1.0-(1-exp(-10))*x-exp(-10*x);
}

int main(int argc, char *argv[]){
	ofstream outfile;
	string filename;
	int expo;
	if(argc <= 1){
    	printf("Error: Please include filename and n dimension of matrix");
    	exit(1);
    }
    else{
        filename = argv[1]; 
        expo = atoi(argv[2]);
    }

	clock_t start, finish;  //  declare start and final time
    start = clock();

	for (int i = 1; i <= expo; i++){
    	int n = (int) pow(10.0,i);

    	string fileout = filename;
    	string argument = to_string(i);
    	fileout.append(argument);
    	double h = (double) 1/(n+1);
    	double hh = h*h;
        n = n-1;
    	mat A = zeros<mat>(n,n);
        mat L, U;
    	vec b(n);  vec x(n);
    	A(0,0) = 2.0;  
    	A(0,1) = -1;  
    	x(0) = h;  
    	b(0) =  hh*f(x(0)); 
    	x(n-1) = x(0)+(n-1)*h; 
    	b(n-1) = hh*f(x(n-1));

    	for (int i = 1; i < n-1; i++) {
    		x(i) = x(i-1)+h; 
			b(i) = hh*f(x(i));
	        A(i,i-1)  = -1.0;
	        A(i,i)    = 2.0;
	        A(i,i+1)  = -1.0;
    	}

    	A(n-1,n-1) = 2.0; 
    	A(n-2,n-1) = -1.0; 
    	A(n-1,n-2) = -1.0;
        // Solve with armadillos solve function
    	vec ans = solve(A, b);
		finish = clock();
        printf("%f\n", ((finish - start)/(double) CLOCKS_PER_SEC ));
        // Write out results in a text file
    	outfile.open(fileout);
       	outfile << setiosflags(ios::showpoint | ios::uppercase);
	    for (int i = 0; i < n; i++) {
			double RelativeError = fabs((u(x(i))-ans(i))/u(x(i)));
			outfile << setw(15) << setprecision(8) << x(i);
			outfile << setw(15) << setprecision(8) << ans(i);
			outfile << setw(15) << setprecision(8) << u(x(i));
         	outfile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
	    }

	    outfile.close();
        // Print out time used for that exponent
        
    }

	return 0;
}