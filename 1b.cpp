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
    	int n = (int) pow(10.0,i);

    	string fileout = filename;
    	string argument = to_string(i);
    	fileout.append(argument);
    	double h = (double) 1/(n+1);
    	double hh = h*h;

        n = n-1;				// Setup n
    	mat A = zeros<mat>(n,n);
    	vec b(n);  vec x(n);
    	A(0,1) = -1;  
    	x(0) = h;  
    	b(0) =  hh*f(x(0)); 
    	x(n-1) = x(0)+(n-1)*h; 
    	b(n-1) = hh*f(x(n-1));

    	clock_t start, finish;  //  declare start and final time for each exponent to test the time of the algorithm
    	start = clock();
    	for (int i = 1; i < n-1; i++) {
    		x(i) = x(i-1)+h; 
			b(i) = hh*f(x(i));
	        A(i,i-1)  = -1.0;
	        A(i,i+1)  = -1.0;
    	}
 
    	A(n-2,n-1) = -1.0; 
    	A(n-1,n-2) = -1.0;
		int maxiter = 10; double diff = 1.0; 
     	double e = 1.0e-7;  int iter = 0;

     	// Setup random vectors for the iterative Jacobi method
      	vec ansNew  = zeros<vec>(n);
      	vec ansOld  = randu<vec>(n);

      	// While loop to get the smallest acceptable difference between the exact value of u and the approximated value of ans
      	while (iter <= maxiter || diff > e){
			ansNew = (b - A*ansOld)*0.5; 
        	iter++; 
        	diff = fabs(sum(ansNew-ansOld)/n);
        	ansOld = ansNew;
      	}
      	vec ans = ansOld;

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
	    finish = clock();
	    printf("%f\n", ((finish - start)/(double) CLOCKS_PER_SEC ));
      	// delete [] x; delete [] d; delete [] b; delete [] ans;
	}

	return 0;
}