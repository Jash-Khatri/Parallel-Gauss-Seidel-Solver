/**
Author: CS19S018
Code for: Computing the numerical solution for the poissons eqn using Gauss-Seidel method + red-black approach with upto 1% error in results using MPI.
Status: 
complete: Yes,
compiles: Yes,
Runs: Yes,
Runs-and-gives-correct-result: Yes.
Note: ( (N-1)/p) should be an integer for the proper functioning of the below code.
where N is the problem size.
p is the number of processes with which the MPI code is executed.
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#define PI 3.14159265358
#define delta 0.1			// can be changed..

double func(double x,double y)		// q
{
 return ( 2*( 2-(x*x)-(y*y) ) );
}

int main(){

	//set the value for N
	int N = ( ( 1.0 - (-1.0) ) / delta ) + 1;

	int myid, size, tag=100;
  	MPI_Status status; 		/* data type that is defined in mpi.h... */
  	MPI_Init(NULL, NULL);

  	MPI_Comm_size(MPI_COMM_WORLD, &size); /* tells about the number of processes */
  	MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* will return the rank or id of the process that called it.  */

	// If the number of processes is 1 then run sequentially
	if(size == 1){

		// Measuring time starts
		double start_time = MPI_Wtime();

		// declare the phi
		double phi[N][N];
		// Initialization for phi
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
			phi[i][j] = 0.0;	
			}
		}
	
		double start = -1.0 ;	// lower bound in range (-1,1) 
	
		int count=0;	//for counting the number of iterations taken by particular method to converge..
		int flag=1; 	// To check if the the solution has converged or not
		
		while(flag){
		count++;
			double x,y;
			// compute the  phi for the red points
			
			for(int i=1;i<N;i++){
			
				for(int j=1;j<N;j++){
	
					if( (i+j)%2 == 0 ){
						//compute the value of phi[i][j] using the equation 12 given in PDF
						if(i < N-1 && j < N-1)
						phi[i][j] = (0.25*(phi[i+1][j]+phi[i][j+1]+phi[i][j-1]+phi[i-1][j]) + ( ((delta*delta)/4.0)*func( (start+(delta*i)), (start+(delta*j)) ) ));
						// for the last row the phi values will be zero as both x = 1 and y = 1						
						else	
						phi[i][j] = 0;
					}
					
				}
				//printf("\n");
			}
	
			// compute the phi for the black points	
		
			for(int i=1;i<N;i++){
	
				for(int j=1;j<N;j++){
	
					if( (i+j)%2 == 1 ){	
						//compute the value of phi[i][j] using the equation 12 given in PDF
						if(i < N-1 && j < N-1)	
						phi[i][j] = (0.25*(phi[i+1][j]+phi[i][j+1]+phi[i][j-1]+phi[i-1][j]) + ( ((delta*delta)/4.0)*func( (start+(delta*i)), (start+(delta*j)) ) ));
						// for the last row the phi values will be zero as both x = 1 and y = 1							
						else
						phi[i][j] = 0;
					}
	
				}
				//printf("\n");
			}
			
			flag = 0; 
			//check if amount of error at each coordinates location is less than 1%
			x = start + delta;
			for(int i=1;i<N-1;i++){
			y = start + delta;
				for(int j=1;j<N-1;j++){
				double exactval = ( (x*x) - 1 ) * ( (y*y) - 1 ) ;
					  if( ( (fabs(phi[i][j]- exactval)  ) / ( exactval ) ) * 100  > 1   )
						{
						flag = 1;
						break;
						}
				//printf("%lf %lf ", exactval , phi[i][j]  );
				y += delta;
				}
				if(flag == 1)
					break;
				//printf("\n");
			x += delta;
			}
		
		}

		// Stop Measuring time
		double end_time = MPI_Wtime();
			
			// Print the final solution for phi
			for(int i=0;i<N;i++){
				for(int j=0;j<N;j++){
				printf("%lf, ", phi[i][j]);
				}
				printf("\n");
			}
			
			printf("\n");
			// Print the number of itertions need to converge
			printf("Number of iterations needed to converge: %d", count);
			printf("\n");
		
		printf("\n\n");
		// print the total time taken for the computation
		printf("Total time taken to compute the derivative is %lf\n", end_time-start_time );

		MPI_Finalize();
		return 0;	
		}


	// Measuring time starts
	double start_time = MPI_Wtime();

	// declare the phi
	double phi[N/size][N];			//assuming that N/size is integer.. Each process will calc phi values for N/size rows..
	// Initialization for phi
	for(int i=0;i<(N/size);i++){
		for(int j=0;j<N;j++){
		phi[i][j] = 0.0;	
		}
	}

	double xstart = -1.0 + ( (N/size)*delta*myid ) + delta;		// computing the lower bound x value for each process
	double ystart = -1.0; 						// computing the lower bound y value for each process

	int count=0;	//for counting the number of iterations taken by particular method to converge..
	int flag=1; 	// To check if the the solution has converged or not

	double ghostpts1[N];				// Arrays for storing the ghost/Virtual points
	for(int i=0;i<N;i++){
		ghostpts1[i] = 0.0;
	}
	double ghostpts2[N];	
	for(int i=0;i<N;i++){
		ghostpts2[i] = 0.0;
	}	

	// Loop until the error in numerical solution is >1%. Number of steps for convergence. 	
	while(flag){
	count++;
		double x,y;

		// compute the  phi for the red points
		// each process computes the phi values for each rows assiged to it in each iteration..
		for(int i=0;i<(N/size);i++){
			// for the boundary rows we need to compute the values using ghost points based on the equation 13 and 14 given in PDF
			if(i == 0 || i == (N/size)-1){
				if(i == 0){
					if(myid == 0){
						for(int j=1;j<N-1;j++){
							// Based on the eqaution 14 given in the PDF
							if( (i+j)%2 == (myid%2)? 0 : 1 )
							phi[i][j] = (0.25*(phi[i+1][j]+phi[i][j+1]+phi[i][j-1]+0.0) + ( ((delta*delta)/4.0)*func( (xstart+(delta*i)), (ystart+(delta*j)) ) ));
						}
					}	
					else{
						for(int j=1;j<N-1;j++){
							// Based on the eqaution 14 given in the PDF
							if( (i+j)%2 == (myid%2)? 0 : 1 )
							phi[i][j] = (0.25*(phi[i+1][j]+phi[i][j+1]+phi[i][j-1]+ghostpts1[j]) + ( ((delta*delta)/4.0)*func( (xstart+(delta*i)), (ystart+(delta*j)) ) ));	
						}
					}			
				}
				if( i == ( (N/size)-1 ) ){
					if(myid == (size-1) ){
						for(int j=1;j<N-1;j++){
							phi[i][j] = 0.0;					
						}
					}
					else{
						for(int j=1;j<N-1;j++){
						// Based on the eqaution 13 given in the PDF						
						if( (i+j)%2 == (myid%2)? 0 : 1 )
						phi[i][j] = (0.25*(ghostpts2[j]+phi[i][j+1]+phi[i][j-1]+phi[i-1][j]) + ( ((delta*delta)/4.0)*func( (xstart+(delta*i)), (ystart+(delta*j)) ) ));		
						}
					}
				}
			}
			// If it is a internal row of the process then compute phi values for each row using eqaution 12
			else{
				for(int j=1;j<N;j++){

					if( (i+j)%2 == (myid%2)? 0 : 1 ){
						if(i < N-1 && j < N-1){
						// Based on equation 12						
						phi[i][j] = (0.25*(phi[i+1][j]+phi[i][j+1]+phi[i][j-1]+phi[i-1][j]) + ( ((delta*delta)/4.0)*func( (xstart+(delta*i)), (ystart+(delta*j)) ) ));
						}
						else
						phi[i][j] = 0;
					}
				
				}
			}
			
			//printf("\n");
		}

		// update the ghost point values using the other process data
		if(myid == 0){
			// Perform the communication with neighbouring process to get the ghost/virtual point values.
			if(myid%2 == 0){
				MPI_Send(phi[(N/size)-1], N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts2, N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
			}	
			else{
				MPI_Recv(ghostpts2, N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[(N/size)-1], N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);
			}	
		}
		if( myid == (size-1) ){
			// Perform the communication with neighbouring process to get the ghost/virtual point values.
			if(myid%2 == 0){
				MPI_Send(phi[0], N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts1, N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
			}	
			else{
				MPI_Recv(ghostpts1, N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[0], N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
			}	
		}
		if( myid != 0 && myid != (size-1) ){
			// Perform the communication with neighbouring process to get the ghost/virtual point values.
			if(myid%2 == 0){
				MPI_Send(phi[0], N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts1, N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[(N/size)-1], N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts2, N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
			}	
			else{
				MPI_Recv(ghostpts2, N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[(N/size)-1], N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts1, N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[0], N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
			}	
		}

		// compute the phi for the black points	
		// each process computes the phi values for each rows assiged to it in each iteration..
		for(int i=0;i<(N/size);i++){

			if(i == 0 || i == (N/size)-1){
				if(i == 0){
					if(myid == 0){
						for(int j=1;j<N-1;j++){
							// Based on the eqaution 14 given in the PDF
							if( (i+j)%2 == (myid%2)? 1 : 0 )
							phi[i][j] = (0.25*(phi[i+1][j]+phi[i][j+1]+phi[i][j-1]+0.0) + ( ((delta*delta)/4.0)*func( (xstart+(delta*i)), (ystart+(delta*j)) ) ));
						}
					}	
					else{
						for(int j=1;j<N-1;j++){
							// Based on the eqaution 14 given in the PDF
							if( (i+j)%2 == (myid%2)? 1 : 0 )
							phi[i][j] = (0.25*(phi[i+1][j]+phi[i][j+1]+phi[i][j-1]+ghostpts1[j]) + ( ((delta*delta)/4.0)*func( (xstart+(delta*i)), (ystart+(delta*j)) ) ));	
						}
					}			
				}
				if( i == ( (N/size)-1 ) ){
					if(myid == (size-1) ){
						for(int j=1;j<N-1;j++){
							phi[i][j] = 0.0;					
						}
					}
					else{
						for(int j=1;j<N-1;j++){
						// Based on the eqaution 13 given in the PDF
						if( (i+j)%2 == (myid%2)? 1 : 0 )
						phi[i][j] = (0.25*(ghostpts2[j]+phi[i][j+1]+phi[i][j-1]+phi[i-1][j]) + ( ((delta*delta)/4.0)*func( (xstart+(delta*i)), (ystart+(delta*j)) ) ));		
						}
					}
				}
			}
			else{
				for(int j=1;j<N;j++){

					if( (i+j)%2 == (myid%2)? 1 : 0 ){
						if(i < N-1 && j < N-1){
						// Based on equation 12	
						phi[i][j] = (0.25*(phi[i+1][j]+phi[i][j+1]+phi[i][j-1]+phi[i-1][j]) + ( ((delta*delta)/4.0)*func( (xstart+(delta*i)), (ystart+(delta*j)) ) ));
						}
						else
						phi[i][j] = 0;
					}
				}
			}
			//printf("\n");
		}

		// update the ghost point values using the other process data
		if(myid == 0){
			// Perform the communication with neighbouring process to get the ghost/virtual point values.
			if(myid%2 == 0){
				MPI_Send(phi[(N/size)-1], N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts2, N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
			}	
			else{
				MPI_Recv(ghostpts2, N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[(N/size)-1], N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);
			}	
		}
		if( myid == (size-1) ){
			// Perform the communication with neighbouring process to get the ghost/virtual point values.
			if(myid%2 == 0){
				MPI_Send(phi[0], N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts1, N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
			}	
			else{
				MPI_Recv(ghostpts1, N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[0], N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
			}	
		}
		if( myid != 0 && myid != (size-1) ){
			// Perform the communication with neighbouring process to get the ghost/virtual point values.
			if(myid%2 == 0){
				MPI_Send(phi[0], N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts1, N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[(N/size)-1], N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts2, N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
			}	
			else{
				MPI_Recv(ghostpts2, N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[(N/size)-1], N, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);
				MPI_Recv(ghostpts1, N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
				MPI_Send(phi[0], N, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
			}	
		}

		//if(count == 1000)
		//	flag=0;
		
		// calculate the residual error
		if(myid == 0){
			flag = 0; 
			//check if amount of error at each coordinates location is less than 1%
			x = -1.0 + ( (N/size)*delta*myid ) + delta;
			for(int i=0;i<(N/size);i++){
			y = -1.0 + delta;
				for(int j=1;j<N-1;j++){
				double exactval = ( (x*x) - 1 ) * ( (y*y) - 1 ) ;
					  if( ( (fabs(phi[i][j]- exactval)  ) / ( exactval ) ) * 100  > 1   )
						{
						flag = 1;
						break;
						}
				//printf("%lf %lf ", exactval , phi[i][j]  );
				y += delta;
				}
				if(flag == 1)
					break;
			//printf("\n");
			x += delta;
			}
		}
		// send the flag to other process so that they can know when to halt
		MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);
	
	}

	// Measure time ends
	double end_time = MPI_Wtime();
	double time = end_time-start_time;
	double final_time;

	// Collecting all the measured time on the all processes.
	MPI_Allreduce(&time,&final_time,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	// declare the finalphi array to store all the phi values computed on all the process to the process 0.
	double *finalphi;
	if(myid == 0){
		finalphi = (double *)malloc(N*N*sizeof(double));
		for(int i=0;i<N*N;i++){
			finalphi[i] = 0.0;
		}
	}	

	// copying the phi values computed on each processes on to the proc 0
	MPI_Gather(&phi[0][0], (N/size)*N, MPI_DOUBLE, &finalphi[N], (N/size)*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(myid == 0)
	{
		
		// Print the final solution for phi on the process 0..
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
			printf("%lf, ", finalphi[i*N+j]);
			}
			printf("\n");
		}
		
		printf("\n");
		// Print the number of itertions need to converge on the process 0..
		printf("Number of iterations needed to converge: %d", count);
		printf("\n");
		printf("\n\n");
		// print the total time taken for the computation
		printf("Total time taken to compute the derivative is %lf\n", final_time/size );
	}

	MPI_Finalize(); 

return 0;
}


