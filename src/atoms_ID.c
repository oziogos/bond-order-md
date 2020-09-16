// 24/12/12
////////////////////////////////////////////////////////////////////////
// 1. Checks the species of each atom and assigns the following:
// - atomic number Z
// - mass in amu
// - atomic radius in A
// - van der Waals radius in A
// 2. Counts the number of each element in the simulation
////////////////////////////////////////////////////////////////////////
#include<stdio.h>
#include<string.h>
//#include"global_vars.h"
FILE *fplog;
//void atoms_ID(void);
//void atoms_ID(void)
// input: particles, **species, *Z, *mass, *atomicR, *vanderWaalsR, logpath[]
void atoms_ID(char *species[],int particles,int Z[],double mass[],double atomicR[],double vanderWaalsR[],char logpath[]);
void atoms_ID(char *species[],int particles,int Z[],double mass[],double atomicR[],double vanderWaalsR[],char logpath[])
{
	int k;
	// ptable is a XY type matrix. Each row corresponds to a specific
	// atomic number Z (the first row is for hydrogen, the eighth is for
	// oxygen and so on...). The first column stores a counter for the same
	// species particles while the second column stores the largest particle
	// index corresponding to each element. 
	int ptable[119][2];
	for (k=0;k<119;++k){ptable[k][0]=0;}	//initialize sum column to zero
	// particle loop for species ID. Uses strcmp to identify and stores
	// Z, mass, atomicR and vanderWaalsR for each particle.
	for (k=0;k<particles;++k)
	{
		
		if (strcmp(species[k],"H")==0)
		{
			// Hydrogen
			Z[k]=1;mass[k]=1.00794;atomicR[k]=0.31;vanderWaalsR[k]=1.2;
			ptable[Z[k]][0]=ptable[Z[k]][0]+1;
			ptable[Z[k]][1]=k;
		}
		
		if (strcmp(species[k],"C")==0)
		{
			// Carbon
			Z[k]=6;mass[k]=12.0107;atomicR[k]=0.7;vanderWaalsR[k]=1.7;
			ptable[Z[k]][0]=ptable[Z[k]][0]+1;
			ptable[Z[k]][1]=k;
		}
		
		if (strcmp(species[k],"N")==0)
		{
			// Nitrogen
			Z[k]=7;mass[k]=14.0067;atomicR[k]=0.71;vanderWaalsR[k]=1.55;
			ptable[Z[k]][0]=ptable[Z[k]][0]+1;
			ptable[Z[k]][1]=k;
		}
		
		if (strcmp(species[k],"O")==0)
		{
			// Oxygen
			Z[k]=8;mass[k]=15.9994;atomicR[k]=0.66;vanderWaalsR[k]=1.52;
			ptable[Z[k]][0]=ptable[Z[k]][0]+1;
			ptable[Z[k]][1]=k;
		}
		
		if (strcmp(species[k],"Ne")==0)
		{
			// Neon
			Z[k]=10;mass[k]=20.1797;atomicR[k]=0.58;vanderWaalsR[k]=1.54;
			ptable[Z[k]][0]=ptable[Z[k]][0]+1;
			ptable[Z[k]][1]=k;
		}
		
		if (strcmp(species[k],"Si")==0)
		{
			// Silicon
			Z[k]=14;mass[k]=28.0855;atomicR[k]=1.1;vanderWaalsR[k]=2.1;
			ptable[Z[k]][0]=ptable[Z[k]][0]+1;
			ptable[Z[k]][1]=k;
		}
		
		if (strcmp(species[k],"Ar")==0)
		{
			// Argon
			Z[k]=18;mass[k]=39.948;atomicR[k]=1.06;vanderWaalsR[k]=1.88;
			ptable[Z[k]][0]=ptable[Z[k]][0]+1;
			ptable[Z[k]][1]=k;
		}
		
		if (strcmp(species[k],"Ge")==0)
		{
			// Germanium
			Z[k]=32;mass[k]=72.63;atomicR[k]=1.22;vanderWaalsR[k]=2.11;
			ptable[Z[k]][0]=ptable[Z[k]][0]+1;
			ptable[Z[k]][1]=k;
		}
		
	}
	// console output
	printf("Atomic species:\n");
	// log
	fplog=fopen(logpath,"a");
	fprintf(fplog,"Atomic species:\n");
	fclose(fplog);
	for(k=0;k<119;++k)
	{
		if(ptable[k][0]!=0)
		{
			printf("---[%s]: %d atoms\n",species[ptable[k][1]],ptable[k][0]);
			// log
			fplog=fopen(logpath,"a");
			fprintf(fplog,"---[%s]: %d atoms\n",species[ptable[k][1]],ptable[k][0]);
			fclose(fplog);
		}
	}
}
