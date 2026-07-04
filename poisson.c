#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <fenv.h>


#define electron 1.6e-19
const double Kt = (electron*25.9e-3);
const double Nc =2.8e25 ;
const double Nv =1.04e25;
const double Eg = (1.12*electron);
const double Esisio2 = 3.1*electron;
const double epsi = 11.7*8.85e-12;
const double epsio2 = 3.9*8.85e-12;
const double intrinsic2 = 2.25e32;


void thomas( long N, double x[N+1],double a[N+1],double b[N+1],double c[N+1]) {
    /*
     solves Ax = d, where A is a tridiagonal matrix consisting of vectors a, b, c
     X = number of equations
     x[] = initially contains the input, d, and returns x. indexed from [0, ..., X - 1]
     a[] = subdiagonal, indexed from [1, ..., X - 1]
     b[] = main diagonal, indexed from [0, ..., X - 1]
     c[] = superdiagonal, indexed from [0, ..., X - 2]
     scratch[] = scratch space of length X, provided by caller, allowing a, b, c to be const
     not performed in this example: manual expensive common subexpression elimination
     */
	double scratch[N+1];
    scratch[0] = c[0] / b[0];
    x[0] = x[0] / b[0];

    /* loop from 1 to X - 1 inclusive */
    for (int ix = 1; ix < N+1; ix++) {
        if (ix < N){
        scratch[ix] = c[ix] / (b[ix] - a[ix] * scratch[ix - 1]);
        }
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) / (b[ix] - a[ix] * scratch[ix - 1]);
    }

    /* loop from X - 2 to 0 inclusive */
    for (int ix = N-1; ix >= 0; ix--)
        x[ix] -= scratch[ix] * x[ix + 1];
}

double ro(double potential,double Nd){
	//printf("%e %e\n\n",Nc*exp(-1*e*potential/Kt),intrinsic2/(Nc*exp(e*-1*potential/Kt)));
	return(electron*(intrinsic2/(Nc*exp(-1*electron*potential/Kt))-Nc*exp(-1*electron*potential/Kt)+Nd));
}

double jacobian_and_d(long N,double ep[N+1],double potential[N+1],double del,double Nd,double w1){
	double a[N+1], b[N+1], c[N+1], d[N+1],factor,r;
	b[0] = 1;
	a[0] = 0;
	c[0] = 0;
	d[0] = 0;
	b[N] = 1;
	a[N] = 0;
	c[N] = 0;
	d[N] = 0;


	for(long i =1;i<N;i++){
			a[i] = (ep[i]+ep[i-1])*0.5/pow(del,2);
			c[i] = (ep[i]+ep[i+1])*0.5/pow(del,2);
			if (i*del<=w1){
				factor = 0;
				r = 0; 
			}
			else{
				factor =  -1*(pow(electron,2)/Kt)*(intrinsic2/(Nc*exp(-1*electron*potential[i]/Kt))+Nc*exp(-1*electron*potential[i]/Kt));
				r = ro(potential[i],Nd);
			}
			b[i] = (-0.5*(ep[i+1]+2*ep[i]+ep[i-1])/pow(del,2))+factor; 
			d[i] = -1*a[i]*potential[i-1]-c[i]*potential[i+1]+ (0.5*(ep[i+1]+2*ep[i]+ep[i-1])/pow(del,2))*potential[i]+r;
	}
	double norm = 0;
	thomas(N,d,a,b,c);
	//printf("%lf %lf\n",d[0],d[N]);
	for(long j =0;j<N+1;j++){
		potential[j] += d[j];
		norm += d[j]*d[j];

	}
	return norm;
}

void Poisson_solver(long N,double x1[N+1],double potentialsi[N+1],double ep1[N+1],double Nd,double del1,double w1){
	int x = 1000;
	double norm1 = 0;
	while(x--){
		norm1 = jacobian_and_d(N,ep1,potentialsi,del1,Nd,w1);
	}
}
int main(){
	FILE* ptr;
	FILE* ptr1;
	FILE* ptr2;
	FILE* ptr3;
	FILE* ptr4;
	ptr4 = fopen("./Eminmax","w");
	ptr3 = fopen("./data","r");
	ptr = fopen("./poisson_sio2.txt","w");
	ptr1 = fopen("./poisson_si.txt","w");
	ptr2 = fopen("./poisson.txt","w");
	double tox;
	double tsi;
	long N;
	double Vg;
	double Nd;

	fscanf(ptr3,"%lf %lf %ld %lf %lf\n",&tox,&tsi,&N,&Vg,&Nd);
	//printf("%e %e %ld %e %e",tox,tsi,N,Vg,Nd);
	//return 0;
	double phi  = 1*Kt*log(Nc/Nd)/electron;
	double w1 = tox;
	double w2 = tsi;
	double del1 = (w1+w2)/N;
	double potentialsi[N+1];
	double E[N+1];//E_c  
	double x1[N+1];
	double ep1[N+1];
	potentialsi[0] = phi - Vg;
	potentialsi[N] = phi ;
	x1[0] = 0;
	x1[N] = w1+w2;
	ep1[0] = epsio2;
	ep1[N] = epsi;
	E[0] = 0;
	E[N] = 0;
	for(long i=1;i<N;i++){
		x1[i] = i*del1;
		if(x1[i]<=w1){
			potentialsi[i] = phi;
			ep1[i] = epsio2;
		}
		else{
			potentialsi[i] = phi;
			ep1[i] = epsi;
		}
		E[i] =0;
	}

	fprintf(ptr,"%ld %e %e\n",0,x1[0],electron*potentialsi[0]);
	Poisson_solver(N,x1,potentialsi,ep1,Nd,del1,w1);
	long int count = 0;
	int c = 0;
	for(long i =0;i<N+1;i++){
		if(x1[i]<=w1){
			E[i] = potentialsi[i]*electron+Esisio2;
			count++;
			fprintf(ptr,"%ld %e %e\n",count,x1[i],E[i]);
		}
		else{
			E[i] = potentialsi[i]*electron;
			fprintf(ptr1,"%e %e\n",x1[i],E[i]);
		}	
	}
	fprintf(ptr4,"%e %e",E[count]/electron,E[count-1]/electron);
	for(long i =0;i<N+1;i++){
		if(x1[i]<=w1){
			E[i] = potentialsi[i]*electron+Esisio2;
		}
		else{
			E[i] = potentialsi[i]*electron;
		}
		fprintf(ptr2,"%e %e\n",x1[i],E[i]);
		//printf("%e %e\n",x1[i],potentialsi[i]);
	}	
	fprintf(ptr,"%ld %e %e\n",count+1,x1[count],E[count]);
	fclose(ptr);
	fclose(ptr1);
	fclose(ptr2);
}

