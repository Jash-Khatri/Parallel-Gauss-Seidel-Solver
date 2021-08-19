/**
Author: CS19S018
Code for: Computing the analytical solution for the poissons eqn
Status: 
complete: Yes,
compiles: Yes,
Runs: Yes,
Runs-and-gives-correct-result: Yes.
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define PI 3.14159265358		
#define delta 0.1			// Domain size. (can be changed..)

int main(){

	//set the value for N
	int N = ( ( 1.0 - (-1.0) ) / delta ) + 1;	
	// Declare and initialize the phi
	double phi[N][N];
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
		phi[i][j] = 0.0;	
		}
	}

	double start = -1.0;		// lower bound in range (-1,1)

	double x,y;
	x = start; 
	// Computing the exact solution
	for(int i=0;i<N;i++){
	y = start;
		for(int j=0;j<N;j++){
			//printf("(%lf,%lf)",x,y);
			phi[i][j] = ( (x*x) - 1 ) * ( (y*y) - 1 );
			y += delta;
		}
		//printf("\n");
	x += delta;
	}

	// Print the exact solution
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
		printf("%lf, ", phi[i][j]);
		}
		printf("\n");
	}

return 0;
}

