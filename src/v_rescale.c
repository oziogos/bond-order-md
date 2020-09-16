// 13/04/13
////////////////////////////////////////////////////////////////////////
// Velocity rescaling using the Andersen thermostat scheme
////////////////////////////////////////////////////////////////////////

#include"MD.h"

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void v_rescale(double ux[],double uy[],double uz[],double temperature,double pP,int particles,int flagArray[],int f_particles,double mass[]);
void v_rescale(double ux[],double uy[],double uz[],double temperature,double pP,int particles,int flagArray[],int f_particles,double mass[])
{
	int k;
	double myrand1,myrand2,uxcm,uycm,uzcm,temp;
	uxcm=0;uycm=0;uzcm=0;						// center-of-mass velocity components
	for (k=0;k<particles;++k)					// particle loop
	{
		if (flagArray[k]==1)					// flag check
		{
			temp=rand();
			myrand1=(double)temp/RAND_MAX;		// draw random number
			if (myrand1<pP)						// particle if
			{
				myrand1=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);		// Box-Muller method to create gaussian random numbers
				myrand2=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);
				ux[k]=sqrt((kB*temperature)/mass[k])*sqrt(-2*log(myrand1))*cos(2*pi*myrand2);
				uy[k]=sqrt((kB*temperature)/mass[k])*sqrt(-2*log(myrand1))*sin(2*pi*myrand2);
				myrand1=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);
				myrand2=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);
				uz[k]=sqrt((kB*temperature)/mass[k])*sqrt(-2*log(myrand1))*cos(2*pi*myrand2);
			}
			uxcm=uxcm+ux[k];					// center-of-mass velocity sums
			uycm=uycm+uy[k];
			uzcm=uzcm+uz[k];
		}
	}
	// center-of-mass velocity calculation -- divide with f_particles!!
	uxcm=uxcm/f_particles;uycm=uycm/f_particles;uzcm=uzcm/f_particles;
	// remove rigid body translation for mobile particles
	for (k=0;k<particles;++k)
	{
		if (flagArray[k]==1)
		{
			ux[k]=ux[k]-uxcm;
			uy[k]=uy[k]-uycm;
			uz[k]=uz[k]-uzcm;
		}
	}
}
