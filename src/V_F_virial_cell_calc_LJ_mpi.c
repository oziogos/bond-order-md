// 13/04/13
////////////////////////////////////////////////////////////////////////
#include<math.h>
double LJsigma(int Zi);
double LJepsilon(int Zi);
void V_F_virial_cell_calc_LJ_mpi(double x[],double y[],double z[],int particles,int Z[],double Lx,double Ly,double Lz,double LJcut,double sumFx[],double sumFy[],double sumFz[],double *V_res,double *sigmaxx_res,double *sigmayy_res,double *sigmazz_res,int M3,int head[],int list[],int **neighbors,int cellpartarray[],int neighborcellpartarray[],int cellsN,int cells[]);
void V_F_virial_cell_calc_LJ_mpi(double x[],double y[],double z[],int particles,int Z[],double Lx,double Ly,double Lz,double LJcut,double sumFx[],double sumFy[],double sumFz[],double *V_res,double *sigmaxx_res,double *sigmayy_res,double *sigmazz_res,int M3,int head[],int list[],int **neighbors,int cellpartarray[],int neighborcellpartarray[],int cellsN,int cells[])
{
	int i,j;
	int cellpart,neighborcellpart,cellindex,neighborcellindex,headpointer,listpointer,locali,localj;
	double r2,r6,r8,r12,r14,xforce,yforce,zforce,forcecoeff,mycoeff,rcut2,sigma,epsilon,sigma6,sigma8,sigma12,sigma14;
	double Fxij,Fyij,Fzij;
	double V;//,virial;
	
	double volume;
	double sigmaxx,sigmayy,sigmazz;
	
	//
	int ic,cellsN_local;
	
	if(cells[cellsN-1]==-1){cellsN_local=cellsN-1;}else{cellsN_local=cellsN;}
	
	volume=Lx*Ly*Lz;
	
	V=0;//virial=0;
	sigmaxx=0;
	sigmayy=0;
	sigmazz=0;
	for (i=0;i<particles;++i){sumFx[i]=0;sumFy[i]=0;sumFz[i]=0;}
	//for (cellindex=0;cellindex<M3;++cellindex)
	for(ic=0;ic<cellsN_local;++ic)
	{
		
		cellindex=cells[ic];
		
		//
		cellpart=0;
		headpointer=head[cellindex];
		if (headpointer!=-1)
		{
			cellpartarray[cellpart]=headpointer;
			cellpart=cellpart+1;
			listpointer=list[headpointer];
			while (listpointer!=-1)
			{
				if (listpointer!=-1)
				{
					cellpartarray[cellpart]=listpointer;
					cellpart=cellpart+1;
					listpointer=list[listpointer];
				}
			}
		}
		//
		if (cellpart!=0)
		{
			for (i=0;i<cellpart;++i)
			{
				for (j=0;j<cellpart;++j)
				{
					if (i!=j)
					{
						locali=cellpartarray[i];
						localj=cellpartarray[j];
						sigma=0.5*(LJsigma(Z[locali])+LJsigma(Z[localj]));
						epsilon=sqrt(LJepsilon(Z[locali])*LJepsilon(Z[localj]));
						xforce=x[localj];
						yforce=y[localj];
						zforce=z[localj];
						if (x[localj]-x[locali]>Lx/2){xforce=x[localj]-Lx;}
						if (x[localj]-x[locali]<-Lx/2){xforce=x[localj]+Lx;}
						if (y[localj]-y[locali]>Ly/2){yforce=y[localj]-Ly;}
						if (y[localj]-y[locali]<-Ly/2){yforce=y[localj]+Ly;}
						if (z[localj]-z[locali]>Lz/2){zforce=z[localj]-Lz;}
						if (z[localj]-z[locali]<-Lz/2){zforce=z[localj]+Lz;}
						r2=(xforce-x[locali])*(xforce-x[locali])+(yforce-y[locali])*(yforce-y[locali])+(zforce-z[locali])*(zforce-z[locali]);
						rcut2=LJcut*LJcut*sigma*sigma;
						if (r2<=rcut2)// cutoff if
						{
							r6=r2*r2*r2;// r^6 calculation
							r8=r6*r2;// r^8 calculation
							r12=r6*r6;// r^12 calculation
							r14=r6*r8;// r^14 calculation
							mycoeff=48*epsilon/(sigma*sigma);// force coeff 1
							sigma6=sigma*sigma*sigma*sigma*sigma*sigma;// sigma^6
							sigma8=sigma6*sigma*sigma;// sigma^8
							sigma12=sigma6*sigma6;// sigma^12
							sigma14=sigma6*sigma8;// sigma^14
							// forces
							forcecoeff=mycoeff*(sigma14/r14-0.5*sigma8/r8);
							Fxij=forcecoeff*(x[locali]-xforce);
							Fyij=forcecoeff*(y[locali]-yforce);
							Fzij=forcecoeff*(z[locali]-zforce);
							sumFx[locali]=sumFx[locali]+Fxij;
							sumFy[locali]=sumFy[locali]+Fyij;
							sumFz[locali]=sumFz[locali]+Fzij;
							V=V+0.5*4*epsilon*(-sigma6/r6+sigma12/r12);// potential energy sum
							//virial=virial-0.5*Fxij*(xforce-x[locali])-0.5*Fyij*(yforce-y[locali])-0.5*Fzij*(zforce-z[locali]);
							sigmaxx=sigmaxx-0.5*Fxij*(xforce-x[locali])/volume;
							sigmayy=sigmayy-0.5*Fyij*(yforce-y[locali])/volume;
							sigmazz=sigmazz-0.5*Fzij*(zforce-z[locali])/volume;
						}		
					}
				}
			}
		}
		
		neighborcellpart=0;
		for (neighborcellindex=0;neighborcellindex<26;++neighborcellindex)
		{
			//
			//neighborcellpart=0;
			headpointer=head[neighbors[neighborcellindex][cellindex]];
			if (headpointer!=-1)
			{
				neighborcellpartarray[neighborcellpart]=headpointer;
				neighborcellpart=neighborcellpart+1;
				listpointer=list[headpointer];
				while (listpointer!=-1)
				{
					if (listpointer!=-1)
					{
						neighborcellpartarray[neighborcellpart]=listpointer;
						neighborcellpart=neighborcellpart+1;
						listpointer=list[listpointer];
					}
				}
			}
		}
			//
			if (neighborcellpart!=0)
			{
				for (i=0;i<cellpart;++i)
				{
					for (j=0;j<neighborcellpart;++j)
					{
							locali=cellpartarray[i];
							localj=neighborcellpartarray[j];
							sigma=0.5*(LJsigma(Z[locali])+LJsigma(Z[localj]));
							epsilon=sqrt(LJepsilon(Z[locali])*LJepsilon(Z[localj]));
							xforce=x[localj];
							yforce=y[localj];
							zforce=z[localj];
							if (x[localj]-x[locali]>Lx/2){xforce=x[localj]-Lx;}
							if (x[localj]-x[locali]<-Lx/2){xforce=x[localj]+Lx;}
							if (y[localj]-y[locali]>Ly/2){yforce=y[localj]-Ly;}
							if (y[localj]-y[locali]<-Ly/2){yforce=y[localj]+Ly;}
							if (z[localj]-z[locali]>Lz/2){zforce=z[localj]-Lz;}
							if (z[localj]-z[locali]<-Lz/2){zforce=z[localj]+Lz;}
							r2=(xforce-x[locali])*(xforce-x[locali])+(yforce-y[locali])*(yforce-y[locali])+(zforce-z[locali])*(zforce-z[locali]);
							rcut2=LJcut*LJcut*sigma*sigma;
							if (r2<=rcut2)// cutoff if
							{
								r6=r2*r2*r2;// r^6 calculation
								r8=r6*r2;// r^8 calculation
								r12=r6*r6;// r^12 calculation
								r14=r6*r8;// r^14 calculation
								mycoeff=48*epsilon/(sigma*sigma);// force coeff 1
								sigma6=sigma*sigma*sigma*sigma*sigma*sigma;// sigma^6
								sigma8=sigma6*sigma*sigma;// sigma^8
								sigma12=sigma6*sigma6;// sigma^12
								sigma14=sigma6*sigma8;// sigma^14
								// forces
								forcecoeff=mycoeff*(sigma14/r14-0.5*sigma8/r8);
								Fxij=forcecoeff*(x[locali]-xforce);
								Fyij=forcecoeff*(y[locali]-yforce);
								Fzij=forcecoeff*(z[locali]-zforce);
								sumFx[locali]=sumFx[locali]+Fxij;
								sumFy[locali]=sumFy[locali]+Fyij;
								sumFz[locali]=sumFz[locali]+Fzij;
								V=V+0.5*4*epsilon*(-sigma6/r6+sigma12/r12);// potential energy sum
								//virial=virial-0.5*Fxij*(xforce-x[locali])-0.5*Fyij*(yforce-y[locali])-0.5*Fzij*(zforce-z[locali]);
								sigmaxx=sigmaxx-0.5*Fxij*(xforce-x[locali])/volume;
								sigmayy=sigmayy-0.5*Fyij*(yforce-y[locali])/volume;
								sigmazz=sigmazz-0.5*Fzij*(zforce-z[locali])/volume;
							}
					}
				}
			}
			//
		
	}
	
	*V_res=V;
	//*virial_res=virial;
	*sigmaxx_res=sigmaxx;
	*sigmayy_res=sigmayy;
	*sigmazz_res=sigmazz;
	
}
