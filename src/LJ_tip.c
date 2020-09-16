// 22/04/12
////////////////////////////////////////////////////////////////////////
#include<math.h>
//#include"global_vars.h"
double LJsigma(int Zi);
double LJepsilon(int Zi);
void LJ_tip(int tip_particles,int particles,double tip_x[],double tip_y[],double tip_z[],int tip_Z[],double x[],double y[],double z[],int Z[],double Lx,double Ly,double Lz,double LJcut,double tip_sumFx[],double tip_sumFy[],double tip_sumFz[],double sumFx[],double sumFy[],double sumFz[],double *tip_V_res);
void LJ_tip(int tip_particles,int particles,double tip_x[],double tip_y[],double tip_z[],int tip_Z[],double x[],double y[],double z[],int Z[],double Lx,double Ly,double Lz,double LJcut,double tip_sumFx[],double tip_sumFy[],double tip_sumFz[],double sumFx[],double sumFy[],double sumFz[],double *tip_V_res)
{
	int i,j;
	double r2,r6,r8,r12,r14,xforce,yforce,zforce,forcecoeff,mycoeff,rcut2,sigma,epsilon,sigma6,sigma8,sigma12,sigma14;
	double Fxvalue,Fyvalue,Fzvalue;
	double tip_V;
	
	tip_V=0;
	for (i=0;i<tip_particles;++i){tip_sumFx[i]=0;tip_sumFy[i]=0;tip_sumFz[i]=0;}
	
	for (i=0;i<particles;++i){
		for (j=0;j<tip_particles;++j){
			
			sigma=0.5*(LJsigma(Z[i])+LJsigma(tip_Z[j]));
			epsilon=sqrt(LJepsilon(Z[i])*LJepsilon(tip_Z[j]));
			xforce=tip_x[j];// PBC
			yforce=tip_y[j];
			zforce=tip_z[j];
			if (tip_x[j]-x[i]>Lx/2){xforce=tip_x[j]-Lx;}
			if (tip_x[j]-x[i]<-Lx/2){xforce=tip_x[j]+Lx;}
			if (tip_y[j]-y[i]>Ly/2){yforce=tip_y[j]-Ly;}
			if (tip_y[j]-y[i]<-Ly/2){yforce=tip_y[j]+Ly;}
			if (tip_z[j]-z[i]>Lz/2){zforce=tip_z[j]-Lz;}
			if (tip_z[j]-z[i]<-Lz/2){zforce=tip_z[j]+Lz;}
			// r^2 calculation
			r2=(xforce-x[i])*(xforce-x[i])+(yforce-y[i])*(yforce-y[i])+(zforce-z[i])*(zforce-z[i]);
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
				Fxvalue=forcecoeff*(x[i]-xforce);
				Fyvalue=forcecoeff*(y[i]-yforce);
				Fzvalue=forcecoeff*(z[i]-zforce);
				sumFx[i]=sumFx[i]+Fxvalue;
				sumFy[i]=sumFy[i]+Fyvalue;
				sumFz[i]=sumFz[i]+Fzvalue;
				tip_sumFx[j]=tip_sumFx[j]-Fxvalue;
				tip_sumFy[j]=tip_sumFy[j]-Fyvalue;
				tip_sumFz[j]=tip_sumFz[j]-Fzvalue;
				tip_V=tip_V+4*epsilon*(-sigma6/r6+sigma12/r12);// potential energy sum
			}
			
			
		}
	}
	
	*tip_V_res=tip_V;
	
}
