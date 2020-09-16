// 13/04/13
////////////////////////////////////////////////////////////////////////
#include<math.h>
//#include<stdio.h>

#define pi 3.14159265
double A_coeff(int Zi);
double B_coeff(int Zi);
double lamda_coeff(int Zi);
double mi_coeff(int Zi);
double beta_coeff(int Zi);
double eta_coeff(int Zi);
double c_coeff(int Zi);
double d_coeff(int Zi);
double h_coeff(int Zi);
double R_coeff(int Zi);
double S_coeff(int Zi);
double tersoff_xi(int Zi,int Zj);

void V_F_virial_cell_calc_tersoff_mpi(double x[],double y[],double z[],int particles,int Z[],double Lx,double Ly,double Lz,double sumFx[],double sumFy[],double sumFz[],double *V_res,double *sigmaxx_res,double *sigmayy_res,double *sigmazz_res,int M3,int head[],int list[],int **neighbors,int cellpartarray[],int neighborcellpartarray[],int allpartarray[],double dx_zeta_III[],double dy_zeta_III[],double dz_zeta_III[],int cellsN,int cells[]);
void V_F_virial_cell_calc_tersoff_mpi(double x[],double y[],double z[],int particles,int Z[],double Lx,double Ly,double Lz,double sumFx[],double sumFy[],double sumFz[],double *V_res,double *sigmaxx_res,double *sigmayy_res,double *sigmazz_res,int M3,int head[],int list[],int **neighbors,int cellpartarray[],int neighborcellpartarray[],int allpartarray[],double dx_zeta_III[],double dy_zeta_III[],double dz_zeta_III[],int cellsN,int cells[])
{
	// linked-cells variables
	int i,j,k;
	int cellpart,neighborcellpart,cellindex,neighborcellindex,headpointer,listpointer,locali,localj;
	int allcellpart,localk;	// for many body interactions
	
	// Tersoff variables
	double xi,A_mean,B_mean,lamda_mean,mi_mean,R_mean,S_mean,fCij,fR,fA,zeta,costheta,gfunction,fCik,bij;
	
	// for the derivatives
	double dx_rij,dy_rij,dz_rij,dx_fR,dy_fR,dz_fR,dx_fA,dy_fA,dz_fA,dx_fCij,dy_fCij,dz_fCij,dx_zeta,dy_zeta,dz_zeta,dx_rik,dy_rik,dz_rik,dx_fCik,dy_fCik,dz_fCik,dx_g,dy_g,dz_g,dx_costheta,dy_costheta,dz_costheta,dx_bij,dy_bij,dz_bij;
	
	// time-saving factors
	double bij_factor,fCij_factor,fCik_factor,g_factor,d_bij_factor;	
	
	// positions and distances
	double xPBCj,yPBCj,zPBCj,xPBCk,yPBCk,zPBCk,rij2,rij,rik2,rik;		
	
	// Tersoff forces
	double dx_zeta_II,dy_zeta_II,dz_zeta_II,rev,dx_bij_II,dy_bij_II,dz_bij_II;
	
	// potential energy
	double V;
	
	// stress related variables
	double sigmaxx,sigmayy,sigmazz;
	double sigmaxx_sum,sigmayy_sum,sigmazz_sum;
	double rijLx,rijLy,rijLz;
	double rikLx,rikLy,rikLz;
	double cosLx,cosLy,cosLz;
	double s_factor_1,s_factor_2;
	
	// mpi related
	int ic,cellsN_local;
	
	if(cells[cellsN-1]==-1){cellsN_local=cellsN-1;}else{cellsN_local=cellsN;}
	
	V=0;
	sigmaxx=0;
	sigmayy=0;
	sigmazz=0;
	
	for (k=0;k<particles;++k)
	{
		sumFx[k]=0;sumFy[k]=0;sumFz[k]=0;
	}
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
		allcellpart=0;
		headpointer=head[cellindex];
		if (headpointer!=-1)
		{
			allpartarray[allcellpart]=headpointer;
			allcellpart=allcellpart+1;
			listpointer=list[headpointer];
			while (listpointer!=-1)
			{
				if (listpointer!=-1)
				{
					allpartarray[allcellpart]=listpointer;
					allcellpart=allcellpart+1;
					listpointer=list[listpointer];
				}
			}
		}
		for (neighborcellindex=0;neighborcellindex<26;++neighborcellindex)
		{
			headpointer=head[neighbors[neighborcellindex][cellindex]];
			if (headpointer!=-1)
			{
				allpartarray[allcellpart]=headpointer;
				allcellpart=allcellpart+1;
				listpointer=list[headpointer];
				while (listpointer!=-1)
				{
					if (listpointer!=-1)
					{
						allpartarray[allcellpart]=listpointer;
						allcellpart=allcellpart+1;
						listpointer=list[listpointer];
					}
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
						xPBCj=x[localj];yPBCj=y[localj];zPBCj=z[localj];// PBCs
						if (x[localj]-x[locali]>Lx/2){xPBCj=x[localj]-Lx;}
						if (x[localj]-x[locali]<-Lx/2){xPBCj=x[localj]+Lx;}
						if (y[localj]-y[locali]>Ly/2){yPBCj=y[localj]-Ly;}
						if (y[localj]-y[locali]<-Ly/2){yPBCj=y[localj]+Ly;}
						if (z[localj]-z[locali]>Lz/2){zPBCj=z[localj]-Lz;}
						if (z[localj]-z[locali]<-Lz/2){zPBCj=z[localj]+Lz;}
						// i-j distance calculation
						rij2=(xPBCj-x[locali])*(xPBCj-x[locali])+(yPBCj-y[locali])*(yPBCj-y[locali])+(zPBCj-z[locali])*(zPBCj-z[locali]);
						rij=sqrt(rij2);
						//S_mean=sqrt(S_[locali]*S_[localj]);// S mean
						S_mean=sqrt(S_coeff(Z[locali])*S_coeff(Z[localj]));
						if (rij<S_mean)// long cutoff interval
						{
							//A_mean=sqrt(A_[locali]*A_[localj]);// various means
							A_mean=sqrt(A_coeff(Z[locali])*A_coeff(Z[localj]));
							//B_mean=sqrt(B_[locali]*B_[localj]);
							B_mean=sqrt(B_coeff(Z[locali])*B_coeff(Z[localj]));
							//R_mean=sqrt(R_[locali]*R_[localj]);
							R_mean=sqrt(R_coeff(Z[locali])*R_coeff(Z[localj]));
							//lamda_mean=(lamda_[locali]+lamda_[localj])/2;
							lamda_mean=(lamda_coeff(Z[locali])+lamda_coeff(Z[localj]))/2;
							//mi_mean=(mi_[locali]+mi_[localj])/2;
							mi_mean=(mi_coeff(Z[locali])+mi_coeff(Z[localj]))/2;
							fR=A_mean*exp(-lamda_mean*rij);// Morse terms
							fA=-B_mean*exp(-mi_mean*rij);
							// distance derivatives
							dx_rij=(x[locali]-xPBCj)/rij;
							dy_rij=(y[locali]-yPBCj)/rij;
							dz_rij=(z[locali]-zPBCj)/rij;
							// Morse derivatives
							dx_fR=-lamda_mean*fR*dx_rij;
							dy_fR=-lamda_mean*fR*dy_rij;
							dz_fR=-lamda_mean*fR*dz_rij;
							dx_fA=-mi_mean*fA*dx_rij;
							dy_fA=-mi_mean*fA*dy_rij;
							dz_fA=-mi_mean*fA*dz_rij;
							// medium cutoff
							if (rij<S_mean && rij>=R_mean)
							{
								fCij=0.5+0.5*cos((pi*(rij-R_mean))/(S_mean-R_mean));
								fCij_factor=((0.5*pi)/(R_mean-S_mean))*sin((pi*(rij-R_mean))/(S_mean-R_mean));// cutoff derivatives
								dx_fCij=fCij_factor*dx_rij;
								dy_fCij=fCij_factor*dy_rij;
								dz_fCij=fCij_factor*dz_rij;
							}
							// short cutoff
							else if (rij<R_mean)
							{
								fCij=1;dx_fCij=0;dy_fCij=0;dz_fCij=0;
								
								fCij_factor=0;
								
							}
							else {fCij=0;dx_fCij=0;dy_fCij=0;dz_fCij=0;}// for -Wall warning
							// @#@#@#@#@ k-loop for the bond order term @#@#@#@#@
							zeta=0;dx_zeta=0;dy_zeta=0;dz_zeta=0;
							dx_zeta_II=0;dy_zeta_II=0;dz_zeta_II=0;
							
							sigmaxx_sum=0;
							sigmayy_sum=0;
							sigmazz_sum=0;
							
							for (k=0;k<allcellpart;++k)
							{
								localk=allpartarray[k];						
						
								dx_zeta_III[localk]=0;										//<------------------
								dy_zeta_III[localk]=0;
								dz_zeta_III[localk]=0;
							
								if (localk!=locali && localk!=localj)
								{
									xPBCk=x[localk];yPBCk=y[localk];zPBCk=z[localk];// PBCs
									if (x[localk]-x[locali]>Lx/2){xPBCk=x[localk]-Lx;}
									if (x[localk]-x[locali]<-Lx/2){xPBCk=x[localk]+Lx;}
									if (y[localk]-y[locali]>Ly/2){yPBCk=y[localk]-Ly;}
									if (y[localk]-y[locali]<-Ly/2){yPBCk=y[localk]+Ly;}
									if (z[localk]-z[locali]>Lz/2){zPBCk=z[localk]-Lz;}
									if (z[localk]-z[locali]<-Lz/2){zPBCk=z[localk]+Lz;}
									// i-k distance calculation
									rik2=(xPBCk-x[locali])*(xPBCk-x[locali])+(yPBCk-y[locali])*(yPBCk-y[locali])+(zPBCk-z[locali])*(zPBCk-z[locali]);
									rik=sqrt(rik2);
									// distance derivatives
									dx_rik=(x[locali]-xPBCk)/rik;
									dy_rik=(y[locali]-yPBCk)/rik;
									dz_rik=(z[locali]-zPBCk)/rik;
									//S_mean=sqrt(S_[locali]*S_[localk]);// S mean
									S_mean=sqrt(S_coeff(Z[locali])*S_coeff(Z[localk]));
									if (rik<S_mean)// long cutoff interval
									{
										//R_mean=sqrt(R_[locali]*R_[localk]);// mean
										R_mean=sqrt(R_coeff(Z[locali])*R_coeff(Z[localk]));
										// medium cutoff
										if (rik<S_mean && rik>=R_mean)
										{
											fCik=0.5+0.5*cos((pi*(rik-R_mean))/(S_mean-R_mean));
											fCik_factor=((0.5*pi)/(R_mean-S_mean))*sin((pi*(rik-R_mean))/(S_mean-R_mean));// cutoff derivatives
											dx_fCik=fCik_factor*dx_rik;
											dy_fCik=fCik_factor*dy_rik;
											dz_fCik=fCik_factor*dz_rik;
										}
										// short cutoff
										else if (rik<R_mean)
										{
											fCik=1;dx_fCik=0;dy_fCik=0;dz_fCik=0;
											
											fCik_factor=0;
											
										}
										else {fCik=0;dx_fCik=0;dy_fCik=0;dz_fCik=0;}// for -Wall warning
										rev=1.0/(rij*rik);
										// cosine calculation
										//costheta=((x[locali]-xPBCj)*(x[locali]-xPBCk)+(y[locali]-yPBCj)*(y[locali]-yPBCk)+(z[locali]-zPBCj)*(z[locali]-zPBCk))/(rij*rik);
										costheta=rev*(x[locali]*x[locali]-x[locali]*xPBCj-x[locali]*xPBCk+xPBCj*xPBCk+y[locali]*y[locali]-y[locali]*yPBCj-y[locali]*yPBCk+yPBCj*yPBCk+z[locali]*z[locali]-z[locali]*zPBCj-z[locali]*zPBCk+zPBCj*zPBCk);
										// cosine derivatives
										//dx_costheta=(2*x[locali]-xPBCj-xPBCk)/(rij*rik)-costheta*((x[locali]-xPBCj)/rij2+(x[locali]-xPBCk)/rik2);
										//dy_costheta=(2*y[locali]-yPBCj-yPBCk)/(rij*rik)-costheta*((y[locali]-yPBCj)/rij2+(y[locali]-yPBCk)/rik2);
										//dz_costheta=(2*z[locali]-zPBCj-zPBCk)/(rij*rik)-costheta*((z[locali]-zPBCj)/rij2+(z[locali]-zPBCk)/rik2);
										dx_costheta=rev*(2*x[locali]-xPBCj-xPBCk-costheta*(rij*dx_rik+rik*dx_rij));
										dy_costheta=rev*(2*y[locali]-yPBCj-yPBCk-costheta*(rij*dy_rik+rik*dy_rij));
										dz_costheta=rev*(2*z[locali]-zPBCj-zPBCk-costheta*(rij*dz_rik+rik*dz_rij));
										// g function
										gfunction=1+(c_coeff(Z[locali])*c_coeff(Z[locali]))/(d_coeff(Z[locali])*d_coeff(Z[locali]))-(c_coeff(Z[locali])*c_coeff(Z[locali]))/(d_coeff(Z[locali])*d_coeff(Z[locali])+(h_coeff(Z[locali])-costheta)*(h_coeff(Z[locali])-costheta));
										// g function derivatives
										g_factor=-2*c_coeff(Z[locali])*c_coeff(Z[locali])*(h_coeff(Z[locali])-costheta)/((d_coeff(Z[locali])*d_coeff(Z[locali])+(h_coeff(Z[locali])-costheta)*(h_coeff(Z[locali])-costheta))*(d_coeff(Z[locali])*d_coeff(Z[locali])+(h_coeff(Z[locali])-costheta)*(h_coeff(Z[locali])-costheta)));
										dx_g=g_factor*dx_costheta;
										dy_g=g_factor*dy_costheta;
										dz_g=g_factor*dz_costheta;
										zeta=zeta+fCik*gfunction;// zeta summation - omega=1
										// zeta derivatives
										dx_zeta=dx_zeta+gfunction*dx_fCik+dx_g*fCik;
										dy_zeta=dy_zeta+gfunction*dy_fCik+dy_g*fCik;
										dz_zeta=dz_zeta+gfunction*dz_fCik+dz_g*fCik;
										
										dx_zeta_II=dx_zeta_II+fCik*g_factor*rev*(xPBCk-x[locali]+costheta*rik*dx_rij);				//<------------
										dy_zeta_II=dy_zeta_II+fCik*g_factor*rev*(yPBCk-y[locali]+costheta*rik*dy_rij);
										dz_zeta_II=dz_zeta_II+fCik*g_factor*rev*(zPBCk-z[locali]+costheta*rik*dz_rij);
										
										dx_zeta_III[localk]=-gfunction*dx_fCik+fCik*g_factor*rev*(xPBCj-x[locali]+costheta*rij*dx_rik);	//<------------
										dy_zeta_III[localk]=-gfunction*dy_fCik+fCik*g_factor*rev*(yPBCj-y[locali]+costheta*rij*dy_rik);
										dz_zeta_III[localk]=-gfunction*dz_fCik+fCik*g_factor*rev*(zPBCj-z[locali]+costheta*rij*dz_rik);
										
										rikLx=(dx_rik*(x[locali]-xPBCk))/Lx;
										rikLy=(dy_rik*(y[locali]-yPBCk))/Ly;
										rikLz=(dz_rik*(z[locali]-zPBCk))/Lz;
										
										cosLx=(2.0*dx_rij*dx_rik-costheta*(dx_rij*dx_rij+dx_rik*dx_rik))/Lx;
										cosLy=(2.0*dy_rij*dy_rik-costheta*(dy_rij*dy_rij+dy_rik*dy_rik))/Ly;
										cosLz=(2.0*dz_rij*dz_rik-costheta*(dz_rij*dz_rij+dz_rik*dz_rik))/Lz;
										
										sigmaxx_sum=sigmaxx_sum+fCik_factor*gfunction*rikLx+fCik*g_factor*cosLx;
										sigmayy_sum=sigmayy_sum+fCik_factor*gfunction*rikLy+fCik*g_factor*cosLy;
										sigmazz_sum=sigmazz_sum+fCik_factor*gfunction*rikLz+fCik*g_factor*cosLz;
										
									}// i-k main cutoff if
								}// k!=i && k!=j if
							}// k-loop
							xi=tersoff_xi(Z[locali],Z[localj]);
							bij_factor=1+pow(beta_coeff(Z[locali])*zeta,eta_coeff(Z[locali]));
							bij=xi*pow(bij_factor,-1/(2*eta_coeff(Z[locali])));
							V=V+0.5*fCij*(fR+bij*fA);
							// ******************* force calculation *******************
							if (zeta!=0)// possible problem with pow(zeta,eta_[i]-1)
							{
								d_bij_factor=-0.5*((pow(beta_coeff(Z[locali]),eta_coeff(Z[locali]))*bij)/bij_factor)*pow(zeta,eta_coeff(Z[locali])-1);
						
								dx_bij=d_bij_factor*dx_zeta;
								dy_bij=d_bij_factor*dy_zeta;
								dz_bij=d_bij_factor*dz_zeta;
								
								dx_bij_II=d_bij_factor*dx_zeta_II;				//<----------
								dy_bij_II=d_bij_factor*dy_zeta_II;
								dz_bij_II=d_bij_factor*dz_zeta_II;
								
								for (k=0;k<allcellpart;++k)
								{
									localk=allpartarray[k];
									sumFx[localk]=sumFx[localk]-0.5*fCij*fA*d_bij_factor*dx_zeta_III[localk];
									sumFy[localk]=sumFy[localk]-0.5*fCij*fA*d_bij_factor*dy_zeta_III[localk];
									sumFz[localk]=sumFz[localk]-0.5*fCij*fA*d_bij_factor*dz_zeta_III[localk];
								} 
							}
							else
							{
								
								d_bij_factor=0;
								
								dx_bij=0;										//<----------
								dy_bij=0;
								dz_bij=0;
								
								dx_bij_II=0;									//<----------
								dy_bij_II=0;
								dz_bij_II=0;
								
							}
							
							sumFx[locali]=sumFx[locali]-0.5*(dx_fCij*(fR+bij*fA)+fCij*(dx_fR+bij*dx_fA+fA*dx_bij));		//<----------
							sumFy[locali]=sumFy[locali]-0.5*(dy_fCij*(fR+bij*fA)+fCij*(dy_fR+bij*dy_fA+fA*dy_bij));
							sumFz[locali]=sumFz[locali]-0.5*(dz_fCij*(fR+bij*fA)+fCij*(dz_fR+bij*dz_fA+fA*dz_bij));
							
							sumFx[localj]=sumFx[localj]-0.5*(-dx_fCij*(fR+bij*fA)+fCij*(-dx_fR-bij*dx_fA+fA*dx_bij_II));	//<----------
							sumFy[localj]=sumFy[localj]-0.5*(-dy_fCij*(fR+bij*fA)+fCij*(-dy_fR-bij*dy_fA+fA*dy_bij_II));
							sumFz[localj]=sumFz[localj]-0.5*(-dz_fCij*(fR+bij*fA)+fCij*(-dz_fR-bij*dz_fA+fA*dz_bij_II));
							
							// stress
							
							rijLx=(dx_rij*(x[locali]-xPBCj))/Lx;
							rijLy=(dy_rij*(y[locali]-yPBCj))/Ly;
							rijLz=(dz_rij*(z[locali]-zPBCj))/Lz;
												
							s_factor_1=-lamda_mean*fR*fCij+fR*fCij_factor+bij*(-mi_mean*fA*fCij+fA*fCij_factor);
							s_factor_2=fA*fCij*d_bij_factor;
														
							sigmaxx=sigmaxx-(0.5/(Ly*Lz))*s_factor_1*rijLx-(0.5/(Ly*Lz))*s_factor_2*sigmaxx_sum;
							sigmayy=sigmayy-(0.5/(Lx*Lz))*s_factor_1*rijLy-(0.5/(Lx*Lz))*s_factor_2*sigmayy_sum;
							sigmazz=sigmazz-(0.5/(Lx*Ly))*s_factor_1*rijLz-(0.5/(Lx*Ly))*s_factor_2*sigmazz_sum;
							
						}// else for the i-j main cutoff if
					}
				}
			}
		}
		//
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
							xPBCj=x[localj];yPBCj=y[localj];zPBCj=z[localj];// PBCs
							if (x[localj]-x[locali]>Lx/2){xPBCj=x[localj]-Lx;}
							if (x[localj]-x[locali]<-Lx/2){xPBCj=x[localj]+Lx;}
							if (y[localj]-y[locali]>Ly/2){yPBCj=y[localj]-Ly;}
							if (y[localj]-y[locali]<-Ly/2){yPBCj=y[localj]+Ly;}
							if (z[localj]-z[locali]>Lz/2){zPBCj=z[localj]-Lz;}
							if (z[localj]-z[locali]<-Lz/2){zPBCj=z[localj]+Lz;}
							// i-j distance calculation
							rij2=(xPBCj-x[locali])*(xPBCj-x[locali])+(yPBCj-y[locali])*(yPBCj-y[locali])+(zPBCj-z[locali])*(zPBCj-z[locali]);
							rij=sqrt(rij2);
							//S_mean=sqrt(S_[locali]*S_[localj]);// S mean
							S_mean=sqrt(S_coeff(Z[locali])*S_coeff(Z[localj]));
							if (rij<S_mean)// long cutoff interval
							{
								//A_mean=sqrt(A_[locali]*A_[localj]);// various means
								A_mean=sqrt(A_coeff(Z[locali])*A_coeff(Z[localj]));
								//B_mean=sqrt(B_[locali]*B_[localj]);
								B_mean=sqrt(B_coeff(Z[locali])*B_coeff(Z[localj]));
								//R_mean=sqrt(R_[locali]*R_[localj]);
								R_mean=sqrt(R_coeff(Z[locali])*R_coeff(Z[localj]));
								//lamda_mean=(lamda_[locali]+lamda_[localj])/2;
								lamda_mean=(lamda_coeff(Z[locali])+lamda_coeff(Z[localj]))/2;
								//mi_mean=(mi_[locali]+mi_[localj])/2;
								mi_mean=(mi_coeff(Z[locali])+mi_coeff(Z[localj]))/2;
								fR=A_mean*exp(-lamda_mean*rij);// Morse terms
								fA=-B_mean*exp(-mi_mean*rij);
								// distance derivatives
								dx_rij=(x[locali]-xPBCj)/rij;
								dy_rij=(y[locali]-yPBCj)/rij;
								dz_rij=(z[locali]-zPBCj)/rij;
								// Morse derivatives
								dx_fR=-lamda_mean*fR*dx_rij;
								dy_fR=-lamda_mean*fR*dy_rij;
								dz_fR=-lamda_mean*fR*dz_rij;
								dx_fA=-mi_mean*fA*dx_rij;
								dy_fA=-mi_mean*fA*dy_rij;
								dz_fA=-mi_mean*fA*dz_rij;
								// medium cutoff
								if (rij<S_mean && rij>=R_mean)
								{
									fCij=0.5+0.5*cos((pi*(rij-R_mean))/(S_mean-R_mean));
									fCij_factor=((0.5*pi)/(R_mean-S_mean))*sin((pi*(rij-R_mean))/(S_mean-R_mean));// cutoff derivatives
									dx_fCij=fCij_factor*dx_rij;
									dy_fCij=fCij_factor*dy_rij;
									dz_fCij=fCij_factor*dz_rij;
								}
								// short cutoff
								else if (rij<R_mean)
								{
									fCij=1;dx_fCij=0;dy_fCij=0;dz_fCij=0;
									
									fCij_factor=0;
									
								}
								else {fCij=0;dx_fCij=0;dy_fCij=0;dz_fCij=0;}// for -Wall warning
								// @#@#@#@#@ k-loop for the bond order term @#@#@#@#@
								zeta=0;dx_zeta=0;dy_zeta=0;dz_zeta=0;
								dx_zeta_II=0;dy_zeta_II=0;dz_zeta_II=0;	
								
								sigmaxx_sum=0;
								sigmayy_sum=0;
								sigmazz_sum=0;
								
								for (k=0;k<allcellpart;++k)
								{
									localk=allpartarray[k];
									
									dx_zeta_III[localk]=0;										//<------------------
									dy_zeta_III[localk]=0;
									dz_zeta_III[localk]=0;
									
									if (localk!=locali && localk!=localj)
									{
										xPBCk=x[localk];yPBCk=y[localk];zPBCk=z[localk];// PBCs
										if (x[localk]-x[locali]>Lx/2){xPBCk=x[localk]-Lx;}
										if (x[localk]-x[locali]<-Lx/2){xPBCk=x[localk]+Lx;}
										if (y[localk]-y[locali]>Ly/2){yPBCk=y[localk]-Ly;}
										if (y[localk]-y[locali]<-Ly/2){yPBCk=y[localk]+Ly;}
										if (z[localk]-z[locali]>Lz/2){zPBCk=z[localk]-Lz;}
										if (z[localk]-z[locali]<-Lz/2){zPBCk=z[localk]+Lz;}
										// i-k distance calculation
										rik2=(xPBCk-x[locali])*(xPBCk-x[locali])+(yPBCk-y[locali])*(yPBCk-y[locali])+(zPBCk-z[locali])*(zPBCk-z[locali]);
										rik=sqrt(rik2);
										// distance derivatives
										dx_rik=(x[locali]-xPBCk)/rik;
										dy_rik=(y[locali]-yPBCk)/rik;
										dz_rik=(z[locali]-zPBCk)/rik;
										//S_mean=sqrt(S_[locali]*S_[localk]);// S mean
										S_mean=sqrt(S_coeff(Z[locali])*S_coeff(Z[localk]));
										if (rik<S_mean)// long cutoff interval
										{
											//R_mean=sqrt(R_[locali]*R_[localk]);// mean
											R_mean=sqrt(R_coeff(Z[locali])*R_coeff(Z[localk]));
											// medium cutoff
											if (rik<S_mean && rik>=R_mean)
											{
												fCik=0.5+0.5*cos((pi*(rik-R_mean))/(S_mean-R_mean));
												fCik_factor=((0.5*pi)/(R_mean-S_mean))*sin((pi*(rik-R_mean))/(S_mean-R_mean));// cutoff derivatives
												dx_fCik=fCik_factor*dx_rik;
												dy_fCik=fCik_factor*dy_rik;
												dz_fCik=fCik_factor*dz_rik;
											}
											// short cutoff
											else if (rik<R_mean)
											{
												fCik=1;dx_fCik=0;dy_fCik=0;dz_fCik=0;
												
												fCik_factor=0;
											
											}
											else {fCik=0;dx_fCik=0;dy_fCik=0;dz_fCik=0;}// for -Wall warning
											rev=1.0/(rij*rik);
											// cosine calculation
											//costheta=((x[locali]-xPBCj)*(x[locali]-xPBCk)+(y[locali]-yPBCj)*(y[locali]-yPBCk)+(z[locali]-zPBCj)*(z[locali]-zPBCk))/(rij*rik);
											costheta=rev*(x[locali]*x[locali]-x[locali]*xPBCj-x[locali]*xPBCk+xPBCj*xPBCk+y[locali]*y[locali]-y[locali]*yPBCj-y[locali]*yPBCk+yPBCj*yPBCk+z[locali]*z[locali]-z[locali]*zPBCj-z[locali]*zPBCk+zPBCj*zPBCk);
											// cosine derivatives
											//dx_costheta=(2*x[locali]-xPBCj-xPBCk)/(rij*rik)-costheta*((x[locali]-xPBCj)/rij2+(x[locali]-xPBCk)/rik2);
											//dy_costheta=(2*y[locali]-yPBCj-yPBCk)/(rij*rik)-costheta*((y[locali]-yPBCj)/rij2+(y[locali]-yPBCk)/rik2);
											//dz_costheta=(2*z[locali]-zPBCj-zPBCk)/(rij*rik)-costheta*((z[locali]-zPBCj)/rij2+(z[locali]-zPBCk)/rik2);
											dx_costheta=rev*(2*x[locali]-xPBCj-xPBCk-costheta*(rij*dx_rik+rik*dx_rij));
											dy_costheta=rev*(2*y[locali]-yPBCj-yPBCk-costheta*(rij*dy_rik+rik*dy_rij));
											dz_costheta=rev*(2*z[locali]-zPBCj-zPBCk-costheta*(rij*dz_rik+rik*dz_rij));
											// g function
											gfunction=1+(c_coeff(Z[locali])*c_coeff(Z[locali]))/(d_coeff(Z[locali])*d_coeff(Z[locali]))-(c_coeff(Z[locali])*c_coeff(Z[locali]))/(d_coeff(Z[locali])*d_coeff(Z[locali])+(h_coeff(Z[locali])-costheta)*(h_coeff(Z[locali])-costheta));
											// g function derivatives
											g_factor=-2*c_coeff(Z[locali])*c_coeff(Z[locali])*(h_coeff(Z[locali])-costheta)/((d_coeff(Z[locali])*d_coeff(Z[locali])+(h_coeff(Z[locali])-costheta)*(h_coeff(Z[locali])-costheta))*(d_coeff(Z[locali])*d_coeff(Z[locali])+(h_coeff(Z[locali])-costheta)*(h_coeff(Z[locali])-costheta)));
											dx_g=g_factor*dx_costheta;
											dy_g=g_factor*dy_costheta;
											dz_g=g_factor*dz_costheta;
											zeta=zeta+fCik*gfunction;// zeta summation - omega=1
											// zeta derivatives
											dx_zeta=dx_zeta+gfunction*dx_fCik+dx_g*fCik;
											dy_zeta=dy_zeta+gfunction*dy_fCik+dy_g*fCik;
											dz_zeta=dz_zeta+gfunction*dz_fCik+dz_g*fCik;
											
											dx_zeta_II=dx_zeta_II+fCik*g_factor*rev*(xPBCk-x[locali]+costheta*rik*dx_rij);				//<------------
											dy_zeta_II=dy_zeta_II+fCik*g_factor*rev*(yPBCk-y[locali]+costheta*rik*dy_rij);
											dz_zeta_II=dz_zeta_II+fCik*g_factor*rev*(zPBCk-z[locali]+costheta*rik*dz_rij);
											
											dx_zeta_III[localk]=-gfunction*dx_fCik+fCik*g_factor*rev*(xPBCj-x[locali]+costheta*rij*dx_rik);	//<------------
											dy_zeta_III[localk]=-gfunction*dy_fCik+fCik*g_factor*rev*(yPBCj-y[locali]+costheta*rij*dy_rik);
											dz_zeta_III[localk]=-gfunction*dz_fCik+fCik*g_factor*rev*(zPBCj-z[locali]+costheta*rij*dz_rik);
											
											rikLx=(dx_rik*(x[locali]-xPBCk))/Lx;
											rikLy=(dy_rik*(y[locali]-yPBCk))/Ly;
											rikLz=(dz_rik*(z[locali]-zPBCk))/Lz;
											
											cosLx=(2.0*dx_rij*dx_rik-costheta*(dx_rij*dx_rij+dx_rik*dx_rik))/Lx;
											cosLy=(2.0*dy_rij*dy_rik-costheta*(dy_rij*dy_rij+dy_rik*dy_rik))/Ly;
											cosLz=(2.0*dz_rij*dz_rik-costheta*(dz_rij*dz_rij+dz_rik*dz_rik))/Lz;
											
											sigmaxx_sum=sigmaxx_sum+fCik_factor*gfunction*rikLx+fCik*g_factor*cosLx;
											sigmayy_sum=sigmayy_sum+fCik_factor*gfunction*rikLy+fCik*g_factor*cosLy;
											sigmazz_sum=sigmazz_sum+fCik_factor*gfunction*rikLz+fCik*g_factor*cosLz;
											
										}// i-k main cutoff if
									}// k!=i && k!=j if
								}// k-loop
								xi=tersoff_xi(Z[locali],Z[localj]);
							bij_factor=1+pow(beta_coeff(Z[locali])*zeta,eta_coeff(Z[locali]));
							bij=xi*pow(bij_factor,-1/(2*eta_coeff(Z[locali])));
							V=V+0.5*fCij*(fR+bij*fA);
							// ******************* force calculation *******************
							if (zeta!=0)// possible problem with pow(zeta,eta_[i]-1)
							{
							
								d_bij_factor=-0.5*((pow(beta_coeff(Z[locali]),eta_coeff(Z[locali]))*bij)/bij_factor)*pow(zeta,eta_coeff(Z[locali])-1);
						
								dx_bij=d_bij_factor*dx_zeta;
								dy_bij=d_bij_factor*dy_zeta;
								dz_bij=d_bij_factor*dz_zeta;
								
								dx_bij_II=d_bij_factor*dx_zeta_II;				//<----------
								dy_bij_II=d_bij_factor*dy_zeta_II;
								dz_bij_II=d_bij_factor*dz_zeta_II;
								
								for (k=0;k<allcellpart;++k)
								{
									localk=allpartarray[k];
									sumFx[localk]=sumFx[localk]-0.5*fCij*fA*d_bij_factor*dx_zeta_III[localk];
									sumFy[localk]=sumFy[localk]-0.5*fCij*fA*d_bij_factor*dy_zeta_III[localk];
									sumFz[localk]=sumFz[localk]-0.5*fCij*fA*d_bij_factor*dz_zeta_III[localk];
								} 
							}
							else
							{
								
								d_bij_factor=0;
														
								dx_bij=0;										//<----------
								dy_bij=0;
								dz_bij=0;
								
								dx_bij_II=0;									//<----------
								dy_bij_II=0;
								dz_bij_II=0;
								
							}
							
							sumFx[locali]=sumFx[locali]-0.5*(dx_fCij*(fR+bij*fA)+fCij*(dx_fR+bij*dx_fA+fA*dx_bij));		//<----------
							sumFy[locali]=sumFy[locali]-0.5*(dy_fCij*(fR+bij*fA)+fCij*(dy_fR+bij*dy_fA+fA*dy_bij));
							sumFz[locali]=sumFz[locali]-0.5*(dz_fCij*(fR+bij*fA)+fCij*(dz_fR+bij*dz_fA+fA*dz_bij));
							
							sumFx[localj]=sumFx[localj]-0.5*(-dx_fCij*(fR+bij*fA)+fCij*(-dx_fR-bij*dx_fA+fA*dx_bij_II));	//<----------
							sumFy[localj]=sumFy[localj]-0.5*(-dy_fCij*(fR+bij*fA)+fCij*(-dy_fR-bij*dy_fA+fA*dy_bij_II));
							sumFz[localj]=sumFz[localj]-0.5*(-dz_fCij*(fR+bij*fA)+fCij*(-dz_fR-bij*dz_fA+fA*dz_bij_II));
							
							// stress
							
							rijLx=(dx_rij*(x[locali]-xPBCj))/Lx;
							rijLy=(dy_rij*(y[locali]-yPBCj))/Ly;
							rijLz=(dz_rij*(z[locali]-zPBCj))/Lz;
												
							s_factor_1=-lamda_mean*fR*fCij+fR*fCij_factor+bij*(-mi_mean*fA*fCij+fA*fCij_factor);
							s_factor_2=fA*fCij*d_bij_factor;
							
							sigmaxx=sigmaxx-(0.5/(Ly*Lz))*s_factor_1*rijLx-(0.5/(Ly*Lz))*s_factor_2*sigmaxx_sum;
							sigmayy=sigmayy-(0.5/(Lx*Lz))*s_factor_1*rijLy-(0.5/(Lx*Lz))*s_factor_2*sigmayy_sum;
							sigmazz=sigmazz-(0.5/(Lx*Ly))*s_factor_1*rijLz-(0.5/(Lx*Ly))*s_factor_2*sigmazz_sum;
							
							}// else for the i-j main cutoff if
						//}				
						
					}
				}
			}
			//
		
	}
	
	*V_res=V;
	*sigmaxx_res=sigmaxx;
	*sigmayy_res=sigmayy;
	*sigmazz_res=sigmazz;

}

