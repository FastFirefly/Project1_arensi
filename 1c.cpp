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

	

	for (int i = 1; i <= expo; i++){

    	int n = (int) pow(10.0,i);		// Setup n

    	string fileout = filename;		// Create filename
    	string argument = to_string(i);
    	fileout.append(argument);

    	double h = (double) 1/(n+1);	// Setup h and h*h
    	double hh = h*h;

    	n = n-1;				// Setup n

    	// Intialize variables and set values
    	vec d(n+1), ans(n+1), b(n+1), x(n+1);
    	x(0) = h; x(n) = 1.0;  
    	b(0) = hh*f(x(0)); b(n) = hh*f(x(n)); 
    	d(0) = 2.0;  d(n) = 2.0;

    	clock_t start, finish;  		//  declare start and final time for each exponent to test the time of the algorithm
    	start = clock();
    	for (int i = 1; i < n; i++){ 			// Setup up values along the diagonal and for the x and b value
			d(i) = (i+1.0)/( (double) i);  
	        x(i) = i*h;
			b(i) = hh*f(x(i));
		}

		// Forward substitution
		for (int i = 2; i < n; i++) {
			b(i) = b(i) + b(i-1)/d(i-1);
      	}

      	// Backward substitution
      	ans(n-1) = b(n-1)/d(n-1); 
    	for (int i = n-2; i > 0; i--){ 
        		ans(i) = (b(i)+ans(i+1))/d(i);
        }
        
		// Calculate the relative error and write the results in a text file
        outfile.open(fileout);
        outfile << setiosflags(ios::showpoint | ios::uppercase);
		for (int i = 1; i < n; i++) {
			double RelativeError = fabs((u(x(i))-ans(i))/u(x(i)));
			outfile << setw(15) << setprecision(8) << x(i);
			outfile << setw(15) << setprecision(8) << ans(i);
			outfile << setw(15) << setprecision(8) << u(x(i));
	       	outfile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
	    }
	    outfile.close();
	    
	    // Print out time used for that exponent
	    finish = clock();
	    printf("%f\n", ((finish - start)/(double) CLOCKS_PER_SEC ));
	}
	return 0;
}