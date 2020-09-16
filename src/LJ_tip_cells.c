#include<math.h>
double LJsigma(int Zi);
double LJepsilon(int Zi);
void LJ_tip_cells(int tip_particles,int particles,double tip_x[],double tip_y[],double tip_z[],int tip_Z[],double x[],double y[],double z[],int Z[],double Lx,double Ly,double Lz,double LJcut,double tip_sumFx[],double tip_sumFy[],double tip_sumFz[],double sumFx[],double sumFy[],double sumFz[],double *tip_V_res,int M3,int tip_head[],int tip_list[],int head[],int list[],int **neighbors,int tip_cellpartarray[],int tip_neighborcellpartarray[],int cellpartarray[]);
void LJ_tip_cells(int tip_particles,int particles,double tip_x[],double tip_y[],double tip_z[],int tip_Z[],double x[],double y[],double z[],int Z[],double Lx,double Ly,double Lz,double LJcut,double tip_sumFx[],double tip_sumFy[],double tip_sumFz[],double sumFx[],double sumFy[],double sumFz[],double *tip_V_res,int M3,int tip_head[],int tip_list[],int head[],int list[],int **neighbors,int tip_cellpartarray[],int tip_neighborcellpartarray[],int cellpartarray[])
{
	int i,j;
	int cellpart,tip_cellpart,tip_neighborcellpart,cellindex,neighborcellindex,headpointer,listpointer,locali,localj;
	double r2,r6,r8,r12,r14,xforce,yforce,zforce,forcecoeff,mycoeff,rcut2,sigma,epsilon,sigma6,sigma8,sigma12,sigma14;
	double Fxij,Fyij,Fzij;
	double tip_V;
	
	tip_V=0;
	for (i=0;i<tip_particles;++i){tip_sumFx[i]=0;tip_sumFy[i]=0;tip_sumFz[i]=0;}
	
	for (cellindex=0;cellindex<M3;++cellindex)
	{
		// i: substrate
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
		// j: probe
		tip_cellpart=0;
		headpointer=tip_head[cellindex];
		if (headpointer!=-1)
		{
			tip_cellpartarray[tip_cellpart]=headpointer;
			tip_cellpart=tip_cellpart+1;
			listpointer=tip_list[headpointer];
			while (listpointer!=-1)
			{
				if (listpointer!=-1)
				{
					tip_cellpartarray[tip_cellpart]=listpointer;
					tip_cellpart=tip_cellpart+1;
					listpointer=tip_list[listpointer];
				}
			}
		}
		
		//
		if (cellpart!=0 && tip_cellpart!=0)
		{
			for (i=0;i<cellpart;++i)
			{
				for (j=0;j<tip_cellpart;++j)
				{
					//if (i!=j)
					//{
						locali=cellpartarray[i];
						localj=tip_cellpartarray[j];
						sigma=0.5*(LJsigma(Z[locali])+LJsigma(tip_Z[localj]));
						epsilon=sqrt(LJepsilon(Z[locali])*LJepsilon(tip_Z[localj]));
						xforce=tip_x[localj];
						yforce=tip_y[localj];
						zforce=tip_z[localj];
						if (tip_x[localj]-x[locali]>Lx/2){xforce=tip_x[localj]-Lx;}
						if (tip_x[localj]-x[locali]<-Lx/2){xforce=tip_x[localj]+Lx;}
						if (tip_y[localj]-y[locali]>Ly/2){yforce=tip_y[localj]-Ly;}
						if (tip_y[localj]-y[locali]<-Ly/2){yforce=tip_y[localj]+Ly;}
						if (tip_z[localj]-z[locali]>Lz/2){zforce=tip_z[localj]-Lz;}
						if (tip_z[localj]-z[locali]<-Lz/2){zforce=tip_z[localj]+Lz;}
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
							tip_sumFx[localj]=tip_sumFx[localj]-Fxij;
							tip_sumFy[localj]=tip_sumFy[localj]-Fyij;
							tip_sumFz[localj]=tip_sumFz[localj]-Fzij;
							tip_V=tip_V+4*epsilon*(-sigma6/r6+sigma12/r12);// potential energy sum
						}		
					//}
				}
			}
		}
		
		//
		tip_neighborcellpart=0;
		for (neighborcellindex=0;neighborcellindex<26;++neighborcellindex)
		{
			// j
			headpointer=tip_head[neighbors[neighborcellindex][cellindex]];
			if (headpointer!=-1)
			{
				tip_neighborcellpartarray[tip_neighborcellpart]=headpointer;
				tip_neighborcellpart=tip_neighborcellpart+1;
				listpointer=tip_list[headpointer];
				while (listpointer!=-1)
				{
					if (listpointer!=-1)
					{
						tip_neighborcellpartarray[tip_neighborcellpart]=listpointer;
						tip_neighborcellpart=tip_neighborcellpart+1;
						listpointer=tip_list[listpointer];
					}
				}
			}
		}
		
		//
		if (tip_neighborcellpart!=0)
		{
			for (i=0;i<cellpart;++i)
			{
				for (j=0;j<tip_neighborcellpart;++j)
				{
						locali=cellpartarray[i];
						localj=tip_neighborcellpartarray[j];
						sigma=0.5*(LJsigma(Z[locali])+LJsigma(tip_Z[localj]));
						epsilon=sqrt(LJepsilon(Z[locali])*LJepsilon(tip_Z[localj]));
						xforce=tip_x[localj];
						yforce=tip_y[localj];
						zforce=tip_z[localj];
						if (tip_x[localj]-x[locali]>Lx/2){xforce=tip_x[localj]-Lx;}
						if (tip_x[localj]-x[locali]<-Lx/2){xforce=tip_x[localj]+Lx;}
						if (tip_y[localj]-y[locali]>Ly/2){yforce=tip_y[localj]-Ly;}
						if (tip_y[localj]-y[locali]<-Ly/2){yforce=tip_y[localj]+Ly;}
						if (tip_z[localj]-z[locali]>Lz/2){zforce=tip_z[localj]-Lz;}
						if (tip_z[localj]-z[locali]<-Lz/2){zforce=tip_z[localj]+Lz;}
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
							tip_sumFx[localj]=tip_sumFx[localj]-Fxij;
							tip_sumFy[localj]=tip_sumFy[localj]-Fyij;
							tip_sumFz[localj]=tip_sumFz[localj]-Fzij;
							tip_V=tip_V+4*epsilon*(-sigma6/r6+sigma12/r12);// potential energy sum
						}
				}
			}
		}
		//
		
	}
	
	*tip_V_res=tip_V;

}
