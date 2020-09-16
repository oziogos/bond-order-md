//
#include<stdio.h>
#include<math.h>

double LJsigma(int Zi);
double S_coeff(int Zi);

FILE *fplog;
void cells_init(int Mx_in, int My_in, int Mz_in, double Lx, double Ly, double Lz, int interaction, double LJcut, int particles, int Z[], char logpath[], int *Mx_out, int *My_out, int *Mz_out);
void cells_init(int Mx_in, int My_in, int Mz_in, double Lx, double Ly, double Lz, int interaction, double LJcut, int particles, int Z[], char logpath[], int *Mx_out, int *My_out, int *Mz_out)
{
	int i;//,j,k,l;
	//int ghostsx,ghostsy,**cells,row,column,ii,jj;
	double LJsigma_max,T3S_max,cut;
	int Mx,My,Mz;//,M2,M3;
	double Lcx,Lcy,Lcz;
	
	
	//define max cutoff
	if(interaction==0)
	{
		LJsigma_max=LJsigma(Z[0]);
		for(i=0;i<particles;++i)
		{
			if(LJsigma_max<LJsigma(Z[i]))
			{
				LJsigma_max=LJsigma(Z[i]);
			}
		}
		cut=LJcut*LJsigma_max;
	}
	if(interaction==1)
	{
		T3S_max=S_coeff(Z[0]);
		for(i=0;i<particles;++i)
		{
			if(T3S_max<S_coeff(Z[i]))
			{
				T3S_max=S_coeff(Z[i]);
			}
		}
		cut=T3S_max;
	}
	
	if(Mx_in==0 || My_in==0 || Mz_in ==0)
	{
		Mx=floor(Lx/cut);
		My=floor(Ly/cut);
		Mz=floor(Lz/cut);
		if(Mx<3){Mx=3;}
		if(My<3){My=3;}
		if(Mz<3){Mz=3;}
	
		Lcx=Lx/Mx;Lcy=Ly/My;Lcz=Lz/Mz;
	

			// console
			printf("Cell info:\n---Automatic grid\n---Number of cells: %d-%d-%d\n",Mx,My,Mz);
			printf("---Cell dimensions: %lfA %lfA %lfA\n",Lcx,Lcy,Lcz);
			// log
			fplog=fopen(logpath,"a");
			fprintf(fplog,"Cell info:\n---Automatic grid\n---Number of cells: %d-%d-%d\n",Mx,My,Mz);
			fprintf(fplog,"---Cell dimensions: %lfA %lfA %lfA\n",Lcx,Lcy,Lcz);
			fclose(fplog);
		
	}
	else
	{
		Mx=Mx_in;
		My=My_in;
		Mz=Mz_in;
		Lcx=Lx/Mx;Lcy=Ly/My;Lcz=Lz/Mz;					// cell dimensions
	

			// console
			printf("Cell info:\n---Custom grid\n---Number of cells: %d-%d-%d\n",Mx,My,Mz);
			printf("---Cell dimensions: %lfA %lfA %lfA\n",Lcx,Lcy,Lcz);
			// log
			fplog=fopen(logpath,"a");
			fprintf(fplog,"Cell info:\n---Custom grid\n---Number of cells: %d-%d-%d\n",Mx,My,Mz);
			fprintf(fplog,"---Cell dimensions: %lfA %lfA %lfA\n",Lcx,Lcy,Lcz);
			fclose(fplog);
		
			if (Lcx<cut)								// check with rcut
			{
				
				// console
				printf("---Warning! Incompatible cell size for x. rcut=%lfA > %lfA\n",cut,Lcx);
				// log
				fplog=fopen(logpath,"a");
				fprintf(fplog,"---Warning! Incompatible cell size for x. rcut=%lfA > %lfA\n",cut,Lcx);
				fclose(fplog);
			}
			if (Lcy<cut)								// check with rcut
			{
				
				// console
				printf("---Warning! Incompatible cell size for y. rcut=%lfA > %lfA\n",cut,Lcy);
				// log
				fplog=fopen(logpath,"a");
				fprintf(fplog,"---Warning! Incompatible cell size for y. rcut=%lfA > %lfA\n",cut,Lcy);
				fclose(fplog);
			}
			if (Lcz<cut)								// check with rcut
			{
				
				// console
				printf("---Warning! Incompatible cell size for z. rcut=%lfA > %lfA\n",cut,Lcz);
				// log
				fplog=fopen(logpath,"a");
				fprintf(fplog,"---Warning! Incompatible cell size for z. rcut=%lfA > %lfA\n",cut,Lcz);
				fclose(fplog);
			}
		
	}
	
	*Mx_out=Mx;
	*My_out=My;
	*Mz_out=Mz;
	
}
