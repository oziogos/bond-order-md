// 20/04/13
////////////////////////////////////////////////////////////////////////
// Notes:
// 1. build with: mpicc -o MDp *.c *.h -lm -O3
//----------------------------------------------------------------------
#include<mpi.h>			// MPI library
// libraries
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
// OS check: important for libraries, system operations and paths
#ifdef _WIN32
#define OS 1
#include<direct.h>
#elif _WIN64
#define OS 1
#include<direct.h>
#elif __linux__
#define OS 2
#include<unistd.h>
#elif __APPLE__
#define OS 3
#include<unistd.h>
#elif __MACH__
#define OS 3
#include<unistd.h>
#else
#define OS 0
#include<unistd.h>
#endif
// include MD header
#include"MD.h"
// all-purposes file pointer
FILE *fp;
// log file pointer
FILE *fplog;

// --- Functions -------------------------------------------------------
// PES_parameters() contains the parametrizations for the implemented FFs

// Identifies atomic species and defines the atomic number Z, the atomic mass and the atomic and van der Waals radii
void atoms_ID(char *species[],int particles,int Z[],double mass[],double atomicR[],double vanderWaalsR[],char logpath[]);
// Velocity rescaling according to the Andersen thermostat scheme
void v_rescale(double ux[],double uy[],double uz[],double temperature,double pP,int particles,int flagArray[],int f_particles,double mass[]);								
// Functions for PES and forces using cell calculations 
// Initialize cells
void cells_init(int Mx_in, int My_in, int Mz_in, double Lx, double Ly, double Lz, int interaction, double LJcut, int particles, int Z[], char logpath[], int *Mx_out, int *My_out, int *Mz_out);										
void build_neighbors(int Mx, int My, int Mz, int **neighbors);
// Routine for LJ potential energy, forces and virial pressure calculations using cells  --  MPI
void V_F_virial_cell_calc_LJ_mpi(double x[],double y[],double z[],int particles,int Z[],double Lx,double Ly,double Lz,double LJcut,double sumFx[],double sumFy[],double sumFz[],double *V_res,double *sigmaxx_res,double *sigmayy_res,double *sigmazz_res,int M3,int head[],int list[],int **neighbors,int cellpartarray[],int neighborcellpartarray[],int cellsN,int cells[]);
// Routine for Tersoff potential energy and forces using cells  --  MPI
void V_F_virial_cell_calc_tersoff_mpi(double x[],double y[],double z[],int particles,int Z[],double Lx,double Ly,double Lz,double sumFx[],double sumFy[],double sumFz[],double *V_res,double *sigmaxx_res,double *sigmayy_res,double *sigmazz_res,int M3,int head[],int list[],int **neighbors,int cellpartarray[],int neighborcellpartarray[],int allpartarray[],double dx_zeta_III[],double dy_zeta_III[],double dz_zeta_III[],int cellsN,int cells[]);
// Rigid tip functions
// The tip-sample L-J interaction function
void LJ_tip(int tip_particles,int particles,double tip_x[],double tip_y[],double tip_z[],int tip_Z[],double x[],double y[],double z[],int Z[],double Lx,double Ly,double Lz,double LJcut,double tip_sumFx[],double tip_sumFy[],double tip_sumFz[],double sumFx[],double sumFy[],double sumFz[],double *tip_V_res);
void LJ_tip_cells(int tip_particles,int particles,double tip_x[],double tip_y[],double tip_z[],int tip_Z[],double x[],double y[],double z[],int Z[],double Lx,double Ly,double Lz,double LJcut,double tip_sumFx[],double tip_sumFy[],double tip_sumFz[],double sumFx[],double sumFy[],double sumFz[],double *tip_V_res,int M3,int tip_head[],int tip_list[],int head[],int list[],int **neighbors,int tip_cellpartarray[],int tip_neighborcellpartarray[],int cellpartarray[]);

// for the cutoff calculation
double LJsigma(int Zi);
double S_coeff(int Zi);

// --- MAIN ------------------------------------------------------------

int main()
{
	int rank,size;
	MPI_Status status;
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	
	// the MASTER rank
	if(rank==0)
	{
	//------------------------------------------------------------------
	// MPI related variables
	int source,dest,id;													// communication variables
	int k,l;															// indices
	int cellsN;															// number of cells per process
	int *cells;															// cells array
	int *p_book,**book;													// cells bookkeeping
	double Vrecv;														// potential energy receiving variable
	double *sumFxrecv,*sumFyrecv,*sumFzrecv;							// total force receiving arrays
	double sigmaxx_recv,sigmayy_recv,sigmazz_recv;						// stress receiving arrays
	//------------------------------------------------------------------
	
	// general
	char *cwd_check;													// getcwd check
	char start_time[cmax_length],stop_time[cmax_length];				// strings for start/stop time
	char buffer[cmax_length];											// char buffer
	char temp[cmax_length];												// char dummy
	char floating[cmax_length];											// data export mask
	char sim_folder[cmax_length];										// simulation folder
	char configpath[cmax_length];										// config path
	char initpath[cmax_length];											// init path
	char resumepath[cmax_length];										// resume path
	char tip_resume_path[cmax_length];									// tip resume path
	char logpath[cmax_length];											// simulation log path
	char dataexportpath[cmax_length];									// thermodynamic data export path 
	char thermoexportpath[cmax_length];									// path for data exports from the thermostated atoms only 
	char P_path[cmax_length];											// pressure data export path
	char csvexportpath[cmax_length];									// csv export path
	char xyzexportpath[cmax_length];									// xyz export path
	char posexportpath[cmax_length];									// coordinates export path
	char velexportpath[cmax_length];									// velocities export path
	char VACFexportpath[cmax_length];									// VACF export file path
	char CoM_path[cmax_length];											// CoM related export file path
	int i,j,n;															// indices
	char *fcheck;														// function check variable for fgets
	int int_random;														// random numbers variable
	double myrand;														// random numbers variable
	double time_spent;													// total time in sec
	
	// config
	int myseed;															// rand() seed
	int steps;															// step to write thermodynamic data and tip data
	int resume,linked_cells;											// resume and linked-list cells method flags
	int csvout,xyzout,posout,velout,console_out,xyz_animation;			// export flags/steps
	int VACFcalcstep;													// VACF calculation interval: starting from myN and back!
	double dtSI;														// timestep in fsec
	int interaction;													// interaction type for the substrate: 0-LJ, 1-T3
	int myN;															// total number of steps
	int tip;															// tip flag
	int tip_version;													// tip version
	int write_res;														// step to write resume
	
	// init - resume
	char title[cmax_length];											// simulation title
	char **species;														// atomic species
	int particles,mobile_particles,thermo_particles;					// particles
	int nres;															// resuming time index
	int *mobileFlag,*thermoFlag;										// atomic flags
	int *Z;																// atomic number Z
	double xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz;						// supercell parameters	
	double *mass;														// atomic mass
	double *vanderWaalsR,*atomicR;										// radii (for Paraview csv files)
	
	// T3
	double *dx_zeta_III,*dy_zeta_III,*dz_zeta_III;						// T3 type-III interaction contribution arrays //<-----------------------
	
	// thermostat
	int Tintervals,*Tstart,*Tstop;										// thermostat(s) initialization variables
	double berendsen_t,contactP,particleP,T;							// thermostat parameters and temperature
	double *Tvalue;														// thermostat(s) initialization variables
	double *thermostat;													// thermostat intervals
	
	// integration
	double dt;															// integration step
	double fpk,fpk_thermo,lamda;										// temperature calculation coefficients and Berendsen lamda
	double V,K,K_thermo,Temperature,Temperature_thermo,pressure;		// thermodynamic quantities
	double Kx,Ky,Kz;													// kinetic energy components
	double Pxx,Pyy,Pzz,sigmaxx,sigmayy,sigmazz;							// pressure and stress components
	double *x_1,*y_1,*z_1;												// positions
	double *ux_1,*uy_1,*uz_1;											// velocities
	double *ax_1,*ay_1,*az_1,*ax_2,*ay_2,*az_2;							// the acceleration arrays
	double *sumFx,*sumFy,*sumFz;										// total force arrays
	
	// VACF
	double VACF,VACFsum,*uxsave,*uysave,*uzsave,VACFinit;				// variables for the VACF calculation - simple version

	// cells
	int cindex;															// the linked-cells id index
	int Mx,My,Mz;														// number of cells
	int M3;																// total number of cells
	int *p;																// array pointer for contiguous allocations -- coupled with neighbors
	int *cellpartarray,*neighborcellpartarray,*head,*list;				// linked-list cells arrays
	int *allpartarray;													// for T3
	int **neighbors;													// neighbors matrix
	double Lcx,Lcy,Lcz;													// cells lengths
	double cut,LJsigma_max,T3S_max;										// cutoff variables
	
	// tip
	char tip_path[cmax_length];//,tip_title[cmax_length];				// tip path and title
	char **tip_species;													// tip atomic species
	int tip_particles;													// number of tip particles
	int *tip_Z;															// tip atomic number Z
	double tip_M,*tip_mass;												// tip mass variables
	double tip_sampleFx,tip_sampleFy,tip_sampleFz;						// tip-sample total force components
	double tip_sampleFx_1,tip_sampleFy_1,tip_sampleFz_1;				// tip-sample total force components
	double xholder,yholder,zholder,x_cm,y_cm,z_cm,x_cm_1,y_cm_1,z_cm_1;	// holder and CoM coords
	double Fxtarget,Fytarget,Fztarget;									// external tip load
	double kx,ky;														// tip-holder spring constants
	double tip_ux,tip_uy,tip_uz;										// tip velocities
	double tip_gamma_x,tip_gamma_y,tip_gamma_z;							// drag coefficients
	double tip_V;														// tip-sample interaction energy
	double *tip_x,*tip_y,*tip_z;										// tip atomic positions
	double *tip_vanderWaalsR,*tip_atomicR;								// tip atomic radii (for Paraview csv files)
	double *tip_sumFx,*tip_sumFy,*tip_sumFz;							// tip-sample interaction arrays
	
	// tip cells
	int tip_linked_cells;												// tip linked-list cells method flag
	int tip_Mx,tip_My,tip_Mz,tip_M3;									// LJ cells														
	int *tip_cellpartarray,*tip_neighborcellpartarray,*tip_head,*tip_list;	// LJ cells arrays
	int *p1;															// array pointer for contiguous allocations
	int	**tip_neighbors;												// LJ neighbors matrix	
	double tip_Lcx,tip_Lcy,tip_Lcz;										// LJ cells lengths	

	// CoM related variables
	double CoMx,CoMy,CoMz,M;											// CoM coordinates and mass
	double CoMux,CoMuy,CoMuz;											// CoM velocity components
	double CoMLx,CoMLy,CoMLz;											// CoM angular velocity components
	
	// barostat variables
	char **bflagX_array,**bflagY_array,**bflagZ_array;					// barostat on/off char flags
	int change_box;														// box rescale flag
	int realloc_flag;													// realloc arrays due to rescale flag
	int Mx_new,My_new,Mz_new;											// new cell numbers
	int Pintervals,*Pstart,*Pstop;										// barostat intervals info
	int *bflagX,*bflagY,*bflagZ;										// arrays to store the on/off info
	double mixx,miyy,mizz;												// Berendsen mi
	double mi_factor,mi_factor_2;										// speed-up factors
	double bulkM;														// bulk modulus
	double berendsen_p;													// Barostat time constant
	double *Pxxvalue,*Pyyvalue,*Pzzvalue;								// barostat intervals pressure values
	double *barostatX,*barostatY,*barostatZ;							// arrays to store the applied pressure

	int M3_0,fullN,*full;

	////////////////////////////////////////////////////////////////////		
			
	// --- Initializations ---------------------------------------------
	
	// Timing initializations
	clock_t begin,end;
	time_t rawtime;
		
	// Read the OS variable and write to console the detected OS type.
	// The OS type determines the path format when reading or writing to a file.
	// In UNIX-based systems, the directories seperator is the '\' character
	// while in Windows is the '/' character.
	if(OS==1){printf("Windows OS detected!\n");}
	else if(OS==2){printf("Linux OS detected!\n");}
	else if(OS==3){printf("MAC OS detected!\n");}
	else{printf("Cannot resolve OS...\nWill run in linux mode.\n");}
	
	// Read the working directory and save to string
	cwd_check=getcwd(sim_folder,cmax_length);
	if(cwd_check==NULL){printf("Could not initialize cwd!\n");exit(-1);}
	
	// Define the log file path
	// All the information printed on screen gets dumped in the log file as well.
	// This is realized using the fplog pointer.
	if (OS==1){sprintf(logpath,"%s\\%s",sim_folder,log_file);}
	else {sprintf(logpath,"%s/%s",sim_folder,log_file);}
	
	// Print the working directory
	// console
	printf("---Current working directory: %s\nLoading files:\n",sim_folder);
	// log
	fplog=fopen(logpath,"w+");	// "w+" flag: create the log file
	fprintf(fplog,"---Current working directory: %s\nLoading files:\n",sim_folder);
	fclose(fplog);

	// Define the path of the config file
	if (OS==1){sprintf(configpath,"%s\\%s",sim_folder,config_file);}
	else{sprintf(configpath,"%s/%s",sim_folder,config_file);}
	
	////////////////////////////////////////////////////////////////////
	
	// Load the config file and initialize the thermostat(s)
	// --------------- read the config file ----------------------------
	// file check
	fp=fopen(configpath,"r");
	if (fp==0 || fp==NULL)
	{
		// console
		printf("Could not locate %s\nThe program will now exit...\n",configpath);
		// log
		fplog=fopen(logpath,"a");
		fprintf(fplog,"Could not locate %s\nThe program will now exit...\n",configpath);
		fclose(fplog);
		exit(-3); // error code
	}		
	else
	{
		// console
		printf("%s ... OK\n",configpath);
		// log
		fplog=fopen(logpath,"a");
		fprintf(fplog,"%s ... OK\n",configpath);
		fclose(fplog);
	}
	
	fplog=fopen(logpath,"a");
	fprintf(fplog,"\n---Config file parameters:\n\n");
	
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&tip);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&tip_version);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&tip_linked_cells);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&tip_Mx);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&tip_My);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&tip_Mz);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&interaction);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&linked_cells);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&myseed);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&resume);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&csvout);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&xyzout);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&xyz_animation);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&T);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&dtSI);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&myN);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&console_out);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&steps);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&write_res);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&VACFcalcstep);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&posout);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&velout);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&Mx);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&My);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&Mz);
	fprintf(fplog,"%s",buffer);
			
	// ---------------- I/O paths --------------------------------------
	// path format depending on the OS type
	if (OS==1)
	{
		sprintf(initpath,"%s\\%s",sim_folder,init_file);				// initial configuration file
		sprintf(dataexportpath,"%s\\out\\%s",sim_folder,out_file);		// K, V, T and virial pressure
		sprintf(VACFexportpath,"%s\\out\\%s",sim_folder,vacf_file);		// VACF
		sprintf(csvexportpath,"%s\\out\\csvout\\",sim_folder);			// csv files directory
		sprintf(xyzexportpath,"%s\\out\\xyzout\\",sim_folder);			// xyz files directory
		sprintf(resumepath,"%s\\out\\%s",sim_folder,resume_file);		// resuming configuration file
		sprintf(thermoexportpath,"%s\\out\\%s",sim_folder,thermo_file);	// temperature of the thermostated atoms
		sprintf(posexportpath,"%s\\out\\posout\\",sim_folder);			// positions directory
		sprintf(velexportpath,"%s\\out\\velout\\",sim_folder);			// velocities directory
		sprintf(CoM_path,"%s\\out\\%s",sim_folder,CoM_file);			// 
		sprintf(P_path,"%s\\out\\%s",sim_folder,P_file);				// 
	}
	if (OS==2||OS==3||OS==0)
	{
		sprintf(initpath,"%s/%s",sim_folder,init_file);					// initial configuration file
		sprintf(dataexportpath,"%s/out/%s",sim_folder,out_file);		// K, V, T and virial pressure
		sprintf(VACFexportpath,"%s/out/%s",sim_folder,vacf_file);		// VACF
		sprintf(csvexportpath,"%s/out/csvout/",sim_folder);				// csv files directory
		sprintf(xyzexportpath,"%s/out/xyzout/",sim_folder);				// xyz files directory
		sprintf(resumepath,"%s/out/%s",sim_folder,resume_file);			// resuming configuration file
		sprintf(thermoexportpath,"%s/out/%s",sim_folder,thermo_file);	// temperature of the thermostated atoms
		sprintf(posexportpath,"%s/out/posout/",sim_folder);				// positions directory
		sprintf(velexportpath,"%s/out/velout/",sim_folder);				// velocities directory
		sprintf(CoM_path,"%s/out/%s",sim_folder,CoM_file);				// 
		sprintf(P_path,"%s/out/%s",sim_folder,P_file);					// 
	}
	
	//------thermostat_initialization-----------------------------------
	// thermostat intervals
	thermostat=(double*)malloc((myN+1)*sizeof(double));
	for (i=0;i<=myN;++i)
	{
		thermostat[i]=-1;			// set to -1: NVE by default
	}
	
	fprintf(fplog,"\n---Thermostat info:\n\n");
	
	fcheck=fgets(buffer,cmax_length,fp);	// black line
	
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&contactP);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&particleP);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&berendsen_t);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&Tintervals);
	fprintf(fplog,"%s",buffer);
	
	Tvalue=(double*)malloc(Tintervals*sizeof(double));
	Tstart=(int*)malloc(Tintervals*sizeof(int));
	Tstop=(int*)malloc(Tintervals*sizeof(int));
	for (i=0;i<Tintervals;++i)
	{
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%lf",&Tstart[i],&Tstop[i],&Tvalue[i]);
		fprintf(fplog,"%s",buffer);
		
		if (Tstop[i]>myN)
		{
			// console
			printf("Thermostat interval exceeds simulation interval. The program will now exit.\n");
			// log
			fprintf(fplog,"Thermostat interval exceeds simulation interval. The program will now exit.\n");
			fclose(fplog);
			exit(-43); // error code
		}

	}
	// populate the thermostat array with the appropriate temperatures
	for (i=0;i<Tintervals;++i)
	{
		for (j=Tstart[i];j<=Tstop[i];++j)
		{
			thermostat[j]=Tvalue[i];
		}
	}
	
	free(Tstart);free(Tstop);free(Tvalue);
	//------------------------------------------------------------------
	
	//------barostat_initialization-----------------------------------
	// barostat intervals
	barostatX=(double*)malloc((myN+1)*sizeof(double));
	barostatY=(double*)malloc((myN+1)*sizeof(double));
	barostatZ=(double*)malloc((myN+1)*sizeof(double));
	bflagX=(int*)malloc((myN+1)*sizeof(int));
	bflagY=(int*)malloc((myN+1)*sizeof(int));
	bflagZ=(int*)malloc((myN+1)*sizeof(int));
	for (i=0;i<=myN;++i)
	{
		bflagX[i]=0;			
		bflagY[i]=0;			
		bflagZ[i]=0;			
		barostatX[i]=0;
		barostatY[i]=0;
		barostatZ[i]=0;
	}
	
	fprintf(fplog,"\n---Barostat info:\n\n");
	
	fcheck=fgets(buffer,cmax_length,fp);	// black line
	
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&bulkM);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&berendsen_p);
	fprintf(fplog,"%s",buffer);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&Pintervals);
	fprintf(fplog,"%s",buffer);
	
	Pxxvalue=(double*)malloc(Pintervals*sizeof(double));
	Pyyvalue=(double*)malloc(Pintervals*sizeof(double));
	Pzzvalue=(double*)malloc(Pintervals*sizeof(double));
	Pstart=(int*)malloc(Pintervals*sizeof(int));
	Pstop=(int*)malloc(Pintervals*sizeof(int));
	bflagX_array=(char**)malloc(Pintervals*sizeof(char*));
	bflagY_array=(char**)malloc(Pintervals*sizeof(char*));
	bflagZ_array=(char**)malloc(Pintervals*sizeof(char*));
	for (i=0;i<Pintervals;++i)
	{
		bflagX_array[i]=(char*)malloc(2*sizeof(char));
		bflagY_array[i]=(char*)malloc(2*sizeof(char));
		bflagZ_array[i]=(char*)malloc(2*sizeof(char));
	}
	for (i=0;i<Pintervals;++i)
	{
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%lf\t%lf\t%lf\t%s\t%s\t%s",&Pstart[i],&Pstop[i],&Pxxvalue[i],&Pyyvalue[i],&Pzzvalue[i],bflagX_array[i],bflagY_array[i],bflagZ_array[i]);
		fprintf(fplog,"%s",buffer);
		
		if (Pstop[i]>myN)
		{
			// console
			printf("Barostat interval exceeds simulation interval. The program will now exit.\n");
			// log
			fprintf(fplog,"Barostat interval exceeds simulation interval. The program will now exit.\n");
			fclose(fplog);
			exit(-43); // error code
		}	
	}
	
	fclose(fplog);
	
	//------------------------------------------------------------------
	// barostat arrays
	
	// change bulkM from Mbar to reduced units
	bulkM=bulkM/P_to_Mbar;
	// calculate rescaling factor
	mi_factor=dtSI/(berendsen_p*bulkM);
	// populate the barostat array with the appropriate pressures and express in reduced units
	for (i=0;i<Pintervals;++i)
	{
		for (j=Pstart[i];j<=Pstop[i];++j)
		{
			if(strcmp(bflagX_array[i],"y")==0){bflagX[j]=1;}else{bflagX[j]=0;}
			if(strcmp(bflagY_array[i],"y")==0){bflagY[j]=1;}else{bflagY[j]=0;}
			if(strcmp(bflagZ_array[i],"y")==0){bflagZ[j]=1;}else{bflagZ[j]=0;}
			
			barostatX[j]=Pxxvalue[i]/P_to_Mbar;
			barostatY[j]=Pyyvalue[i]/P_to_Mbar;
			barostatZ[j]=Pzzvalue[i]/P_to_Mbar;
		}
	}
	
	free(Pstart);free(Pstop);free(Pxxvalue);free(Pyyvalue);free(Pzzvalue);
	for(i=0;i<Pintervals;++i){free(bflagX_array[i]);free(bflagY_array[i]);free(bflagZ_array[i]);}
	free(bflagX_array);free(bflagY_array);free(bflagZ_array);
	
	//------------------------------------------------------------------
	// read tip data
	if(tip==1)
	{
		fplog=fopen(logpath,"a");
		
		fprintf(fplog,"\n---Tip info info:\n\n");
		
		fcheck=fgets(buffer,cmax_length,fp);	// black line
		
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&Fxtarget);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&Fytarget);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&Fztarget);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&tip_ux);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&tip_uy);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&tip_uz);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&tip_gamma_x);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&tip_gamma_y);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&tip_gamma_z);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&kx);
		fprintf(fplog,"%s",buffer);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&ky);
		fprintf(fplog,"%s",buffer);
		
		fclose(fplog);
		
	}
	//------------------------------------------------------------------
	fclose(fp);
	
	////////////////////////////////////////////////////////////////////
	
	// File check for the initial configuration file (for fresh simulations)
	// +Initialize nres
	if (resume==0)
	{
		//
		nres=0;
		//
		fp=fopen(initpath,"r");							
		if (fp==0 || fp==NULL)
		{
			// console
			printf("Could not locate %s\nThe program will now exit...\n",initpath);
			// log
			fplog=fopen(logpath,"a");
			fprintf(fplog,"Could not locate %s\nThe program will now exit...\n",initpath);
			fclose(fplog);
			exit(-5); // error code
		}
		else
		{
			// console
			printf("%s ... OK\n",initpath);
			// log
			fplog=fopen(logpath,"a");
			fprintf(fplog,"\n%s ... OK\n",initpath);
			fclose(fplog);
		}
		fclose(fp);
	}

	// Print the interaction type
	// Lennard-Jones: 0 -- the version flag is not utilized
	// T3: 1 -- version: 1
	if (interaction==0)
	{
		// console
		printf("---Running Lennard-Jones\n");
		// log
		fplog=fopen(logpath,"a");
		fprintf(fplog,"---Running Lennard-Jones\n");
		fclose(fplog);
	}				
	if (interaction==1)
	{
		// console
		printf("---Running Tersoff T3\n");
		// log
		fplog=fopen(logpath,"a");
		fprintf(fplog,"---Running Tersoff T3\n");
		fclose(fplog);
	}

	// Create the output directories
	// console
	printf("Creating directories...\n");
	// log
	fplog=fopen(logpath,"a");
	fprintf(fplog,"Creating directories...\n");
	fclose(fplog);
	// out directory
	if (OS==1){sprintf(buffer,"mkdir %s\\out\\",sim_folder);}
	else{sprintf(buffer,"mkdir %s/out/",sim_folder);}
	system(buffer);
	// create directory for the csv files
	sprintf(buffer,"mkdir %s",csvexportpath);system(buffer);
	// create directory for the xyz files
	sprintf(buffer,"mkdir %s",xyzexportpath);system(buffer);
	// create directory for the positions dump files
	if (posout!=0)
	{
		sprintf(buffer,"mkdir %s",posexportpath);
		system(buffer);
	}
	// create directory for the velocity dump files
	if (velout!=0)
	{
		sprintf(buffer,"mkdir %s",velexportpath);
		system(buffer);
	}

	// Set the seed for the random number generator
	// If the seed parameter from the config file is zero, then a random
	// seed is selected via the system clock.
	if(myseed!=0){srand(myseed);}else{srand(time(0));}

	// Scale physical quantities
	dt=dtSI/timescale;	// time step (reduced)

	// Read the configuration data
	// If resume is equal to zero, a fresh simulation is carried out reading
	// the initial configuration file.
	// If resume is equal to one, a resuming simulation is carried out reading
	// the resuming configuration file.
	if (resume==0)
	{
		////////////////////////////////////////////////////////////////
		fp=fopen(initpath,"r");							// open the initial configuration file
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%s",temp,title);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&particles);
		// read supercell parameters	
		fcheck=fgets(buffer,cmax_length,fp);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&xmin);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&xmax);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&ymin);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&ymax);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&zmin);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&zmax);
		// preallocations
		species=(char**)malloc(particles*sizeof(char*));
		for (i=0;i<particles;++i)
		{
			species[i]=(char*)malloc(species_length*sizeof(char));
		}
		x_1=(double*)malloc(particles*sizeof(double));
		y_1=(double*)malloc(particles*sizeof(double));
		z_1=(double*)malloc(particles*sizeof(double));
		mobileFlag=(int*)malloc(particles*sizeof(int));
		thermoFlag=(int*)malloc(particles*sizeof(int));
		// read atomic data
		mobile_particles=0;		// counter for mobile particles
		thermo_particles=0;		// counter for thermostated particles
		fcheck=fgets(buffer,cmax_length,fp);
		for (i=0;i<particles;++i)
		{
			fcheck=fgets(buffer,cmax_length,fp);
			sscanf(buffer,"%s\t%lf\t%lf\t%lf\t%d\t%d",species[i],&x_1[i],&y_1[i],&z_1[i],&mobileFlag[i],&thermoFlag[i]);
			mobile_particles=mobile_particles+mobileFlag[i];		// count mobile particles
			thermo_particles=thermo_particles+thermoFlag[i];		// count thermostated particles
		}
		// end of atomic data reading
		fclose(fp);													// close the initial configuration file
		////////////////////////////////////////////////////////////////
	}	
	if (resume==1)
	{
		// console
		printf("---Resuming simulation from: %s\n",resumepath);
		// log
		fplog=fopen(logpath,"a");
		fprintf(fplog,"---Resuming simulation from: %s\n",resumepath);
		fclose(fplog);
		// read the resuming configuration file
		////////////////////////////////////////////////////////////////
		fp=fopen(resumepath,"r");						
		// file check
		if (fp==0 || fp==NULL)
		{
			// console
			printf("Could not locate %s\nThe program will now exit...\n",resumepath);
			// log
			fplog=fopen(logpath,"a");
			fprintf(fplog,"Could not locate %s\nThe program will now exit...\n",resumepath);
			fclose(fplog);
			exit(-5); // error code
		}
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%s",temp,title);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d",temp,&particles);
		// read supercell parameters	
		fcheck=fgets(buffer,cmax_length,fp);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&xmin);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&xmax);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&ymin);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&ymax);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&zmin);
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",temp,&zmax);
		// preallocations
		species=(char**)malloc(particles*sizeof(char*));
		for (i=0;i<particles;++i)
		{
			species[i]=(char*)malloc(species_length*sizeof(char));
		}
		x_1=(double*)malloc(particles*sizeof(double));
		y_1=(double*)malloc(particles*sizeof(double));
		z_1=(double*)malloc(particles*sizeof(double));
		// preallocate velocity arrays!
		ux_1=(double*)malloc(particles*sizeof(double));	
		uy_1=(double*)malloc(particles*sizeof(double));
		uz_1=(double*)malloc(particles*sizeof(double));
		mobileFlag=(int*)malloc(particles*sizeof(int));
		thermoFlag=(int*)malloc(particles*sizeof(int));
		// read atomic data
		mobile_particles=0;							// counters
		thermo_particles=0;
		fcheck=fgets(buffer,cmax_length,fp);
		for (i=0;i<particles;++i)
		{
			fcheck=fgets(buffer,cmax_length,fp);
			sscanf(buffer,"%s\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf",species[i],&x_1[i],&y_1[i],&z_1[i],&mobileFlag[i],&thermoFlag[i],&ux_1[i],&uy_1[i],&uz_1[i]);
			mobile_particles=mobile_particles+mobileFlag[i];
			thermo_particles=thermo_particles+thermoFlag[i];
		}
		// end of atomic data reading
		fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&nres);
		fclose(fp);									// close the resuming configuration file
		////////////////////////////////////////////////////////////////
	}
	
	// Load function to identify the atomic species
	// The atoms_ID function assigns also the appropriate values to the
	// Z, mass, atomicR and vanderWaalsR arrays.
	// write info
	if(tip==1){printf("Substrate:\n");fplog=fopen(logpath,"a");fprintf(fplog,"Substrate:\n");fclose(fplog);}
	// preallocations
	Z=(int*)malloc(particles*sizeof(int));
	mass=(double*)malloc(particles*sizeof(double));
	atomicR=(double*)malloc(particles*sizeof(double));
	vanderWaalsR=(double*)malloc(particles*sizeof(double));
	// function call
	atoms_ID(species,particles,Z,mass,atomicR,vanderWaalsR,logpath);
	
	// T3 preallocations
	if (interaction==1)
	{
		dx_zeta_III=(double*)malloc(particles*sizeof(double));			//<-----------------------
		dy_zeta_III=(double*)malloc(particles*sizeof(double));
		dz_zeta_III=(double*)malloc(particles*sizeof(double));		
	}	

	// Preallocations
	// particles dependant
	// accelerations
	ax_1=(double*)malloc(particles*sizeof(double));
	ay_1=(double*)malloc(particles*sizeof(double));
	az_1=(double*)malloc(particles*sizeof(double));
	ax_2=(double*)malloc(particles*sizeof(double));
	ay_2=(double*)malloc(particles*sizeof(double));
	az_2=(double*)malloc(particles*sizeof(double));
	// total forces
	sumFx=(double*)malloc(particles*sizeof(double));
	sumFy=(double*)malloc(particles*sizeof(double));
	sumFz=(double*)malloc(particles*sizeof(double));
		
	// 1
	// *** total forces: receive buffers
	sumFxrecv=(double*)malloc(particles*sizeof(double));
	sumFyrecv=(double*)malloc(particles*sizeof(double));
	sumFzrecv=(double*)malloc(particles*sizeof(double));
	// ***
	
	// Coefficients for temperature calculations with fixed mass center
	fpk=kB*(3*mobile_particles-3);	
	fpk_thermo=kB*(3*thermo_particles-3);

	// Supercell parameters
	// Simulation box size
	Lx=xmax-xmin;Ly=ymax-ymin;Lz=zmax-zmin;
	// print supercell info
	// console
	printf("Supecell dimensions:\n---Lx=%lfA\n---Ly=%lfA\n---Lz=%lfA\n",Lx,Ly,Lz);
	// log
	fplog=fopen(logpath,"a");
	fprintf(fplog,"Supecell dimensions:\n---Lx=%lfA\n---Ly=%lfA\n---Lz=%lfA\n",Lx,Ly,Lz);
	fclose(fplog);
	
	// Gaussian thermalization for initial velocities
	// If resume==1, the preallocations and assignments for ux_1,uy_1,uz_1 are already done
	if (resume==0)
	{
		ux_1=(double*)malloc(particles*sizeof(double));
		uy_1=(double*)malloc(particles*sizeof(double));
		uz_1=(double*)malloc(particles*sizeof(double));
		// use Andersen thermostat to draw Gaussian velocities
		v_rescale(ux_1,uy_1,uz_1,T,1,particles,mobileFlag,mobile_particles,mass);
		// fixed particles: set velocities equal to zero
		for (i=0;i<particles;++i)
		{
		// need to set equal to zero because v_rescale() SKIPS particles with mobileFlag==0
		// if resume==1, the velocities are zero in the resume file
			if (mobileFlag[i]==0){ux_1[i]=0;uy_1[i]=0;uz_1[i]=0;}
		}
	}
	
	// VACF initialization
	if (VACFcalcstep==0)
	{
		// console
		printf("---Caution: No VACF calculations\n");
		// log
		fplog=fopen(logpath,"a");
		fprintf(fplog,"---Caution: No VACF calculations\n");
		fclose(fplog);
	}
	if (VACFcalcstep>myN)
	{
		// console
		printf("---Warning!! No VACF calculations - interval exceeds simulation time\n");
		// log
		fplog=fopen(logpath,"a");
		fprintf(fplog,"---Warning!! No VACF calculations - interval exceeds simulation time\n");
		fclose(fplog);
	}
	// the VACF reference time point must be inside the Verlet integration loop;
	// hence the 'VACFcalcstep>0' condition
	if (VACFcalcstep>0 && VACFcalcstep<=myN)
	{
		// array for the VACF values
		//VACF=(double*)malloc(VACFcalcstep*sizeof(double));
		// arrays to save the starting velocities	
		uxsave=(double*)malloc(particles*sizeof(double));		
		uysave=(double*)malloc(particles*sizeof(double));
		uzsave=(double*)malloc(particles*sizeof(double));
	}
	
	// Linked cells initialization
	
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
		cut=LJcutoff*LJsigma_max;
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
	
	//
		
	cells_init(Mx,My,Mz,Lx,Ly,Lz,interaction,LJcutoff,particles,Z,logpath,&Mx,&My,&Mz);
	M3=Mx*My*Mz;
	Lcx=Lx/Mx;Lcy=Ly/My;Lcz=Lz/Mz;
			
	// preallocations
	cellpartarray=(int*)malloc(particles*sizeof(int));// cells
	neighborcellpartarray=(int*)malloc(particles*sizeof(int));
	allpartarray=(int*)malloc(particles*sizeof(int));
	head=(int*)malloc(M3*sizeof(int));// head
	list=(int*)malloc(particles*sizeof(int));// list
	neighbors=(int**)malloc(26*sizeof(int*));
	//for (i=0;i<26;++i){neighbors[i]=(int*)malloc(Mx*My*Mz*sizeof(int));}
	p=(int*)malloc(26*M3*sizeof(int));
	for(i=0;i<26;++i)
	{
		neighbors[i]=p+i*M3;
	}
	build_neighbors(Mx,My,Mz,neighbors);
	/*
	// 2
	// *** preallocate the bookkeeping matrix
	cellsN=(int)floor(M3/size)+1;
	cells=(int*)malloc(cellsN*sizeof(int));	// for the master process
	book=(int**)malloc(size*sizeof(int*));
	p_book=(int*)malloc(size*cellsN*sizeof(int));
	for(i=0;i<size;++i)
	{
		book[i]=p_book+i*cellsN;
	}
	// ***
	// *** populate the bookkeeping matrix
	// initial bookkeeping -- static!!!!!!!!!!!!
	i=0;l=0;
	while(i<M3)
	{
		for(k=0;k<size;++k)
		{
			book[k][cellsN-1]=-1;	// <-----
			if(i<M3)
			{
				book[k][l]=i;
				i=i+1;
			}
		}
		l=l+1;
	}
	// ***
	// *** cells for the master process
		i=0;
		for(j=0;j<cellsN;++j){cells[j]=book[i][j];}
	// ***
	*/		
	// 3
	// *** send initialization data
	for(i=1;i<size;++i)
	{		
		dest=i;
		id=0;
		MPI_Send(&particles, 1, MPI_INT, dest, id, MPI_COMM_WORLD);
		id=1;
		MPI_Send(&M3, 1, MPI_INT, dest, id, MPI_COMM_WORLD);
		id=2;
		MPI_Send(&Z[0], particles, MPI_INT, dest, id, MPI_COMM_WORLD);
		id=3;
		MPI_Send(&neighbors[0][0], 26*M3, MPI_INT, dest, id, MPI_COMM_WORLD);
		id=4;
		MPI_Send(&Lx, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
		id=5;
		MPI_Send(&Ly, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
		id=6;
		MPI_Send(&Lz, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
		
		//id=13;
		//MPI_Send(&book[i][0], cellsN, MPI_INT, dest, id, MPI_COMM_WORLD);	
		
		id=12;
		MPI_Send(&myN, 1, MPI_INT, dest, id, MPI_COMM_WORLD);
		id=19;
		MPI_Send(&interaction, 1, MPI_INT, dest, id, MPI_COMM_WORLD);
		
	}
	// ***

	// Timing
	// starting time
	time(&rawtime);
	begin = clock();
	// write time to string
	sprintf(start_time,"\nSimulation started at: %s",ctime(&rawtime));
	
//(A) --- Tip data import, preallocations and assignments --------------
	if(tip==1)
	{
		if(OS==1){
		sprintf(tip_resume_path,"%s\\out\\%s",sim_folder,tip_resume_file);}
		else{
		sprintf(tip_resume_path,"%s/out/%s",sim_folder,tip_resume_file);}
		
		// Display simulation mode
		// console
		printf("\n---Rigid probe simulation mode\n");
		if(tip_version==1){printf("[version: 1]---Dynamic probe + dynamic substrate\n");}
		if(tip_version==2){printf("[version: 2]---Dynamic probe + frozen substrate\n");}
		if(tip_version==3){printf("[version: 3]---Frozen probe + dynamic substrate\n");}
		// log
		fplog=fopen(logpath,"a");
		fprintf(fplog,"\n---Rigid probe simulation mode\n");
		if(tip_version==1){fprintf(fplog,"[version: 1]---Dynamic probe + dynamic substrate\n");}
		if(tip_version==2){fprintf(fplog,"[version: 2]---Dynamic probe + frozen substrate\n");}
		if(tip_version==3){fprintf(fplog,"[version: 3]---Frozen probe + dynamic substrate\n");}
		fclose(fplog);
		
		// Import tip data
		if(resume==0)	// fresh simulation
		{
			// console
			printf("Loading probe files:\n");
			// log
			fplog=fopen(logpath,"a");
			fprintf(fplog,"Loading probe files:\n");
			fclose(fplog);
			// Tip coordinates
			if(OS==1){sprintf(tip_path,"%s\\%s",sim_folder,tip_file);}
			else{sprintf(tip_path,"%s/%s",sim_folder,tip_file);}
			fp=fopen(tip_path,"r");
			// File check for coords file
			if (fp==0 || fp==NULL)
			{
				// console
				printf("Could not locate %s\nThe program will now exit...\n",tip_path);
				// log
				fplog=fopen(logpath,"a");
				fprintf(fplog,"Could not locate %s\nThe program will now exit...\n",tip_path);
				fclose(fplog);
				exit(-100); // error code
			}
			else
			{
				// console
				printf("%s ... OK\n",tip_path);
				// log
				fplog=fopen(logpath,"a");
				fprintf(fplog,"%s ... OK\n",tip_path);
				fclose(fplog);
			}
			// Read tip coords
			fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&tip_particles);
			fcheck=fgets(buffer,cmax_length,fp);
			tip_x=(double*)malloc(tip_particles*sizeof(double));
			tip_y=(double*)malloc(tip_particles*sizeof(double));
			tip_z=(double*)malloc(tip_particles*sizeof(double));
			tip_species=(char**)malloc(tip_particles*sizeof(char*));
			for (i=0;i<tip_particles;++i)
			{
				tip_species[i]=(char*)malloc(species_length*sizeof(char));
			}
			for(i=0;i<tip_particles;++i)
			{
				fcheck=fgets(buffer,cmax_length,fp);
				sscanf(buffer,"%s\t%lf\t%lf\t%lf",tip_species[i],&tip_x[i],&tip_y[i],&tip_z[i]);
				if(tip_x[i]<xmin){printf("tip_x[%d]<xmin - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_x[i]>xmax){printf("tip_x[%d]>xmax - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_y[i]<ymin){printf("tip_y[%d]<ymin - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_y[i]>ymax){printf("tip_y[%d]>ymax - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_z[i]<zmin){printf("tip_z[%d]<zmin - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_z[i]>zmax){printf("tip_z[%d]>zmax - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				
			}
			fclose(fp);
		}
		
		// Units
		// ...
		//
		// input is in nN
		Fxtarget=Fxtarget/1.60217;
		Fytarget=Fytarget/1.60217;
		Fztarget=Fztarget/1.60217;
		// input is in N/m
		kx=kx/16.0217;
		ky=ky/16.0217;
		// input is in m/s
		tip_ux=tip_ux/9822.67;
		tip_uy=tip_uy/9822.67;
		tip_uz=tip_uz/9822.67;
		// input is in nN*s/m
		tip_gamma_x=(tip_gamma_x*1000000)/(16.0217*10.1805);
		tip_gamma_y=(tip_gamma_y*1000000)/(16.0217*10.1805);
		tip_gamma_z=(tip_gamma_z*1000000)/(16.0217*10.1805);
		
		// Resuming simulation
		if(resume==1)	
		{
			// console
			printf("---Resuming probe simulation\n");
			// log
			fplog=fopen(logpath,"a");
			fprintf(fplog,"---Resuming probe simulation\n");
			fclose(fplog);
			// File check for the tip coords file
			// (contains also the velocity components in reduced units)
			//if(OS==1){sprintf(tip_path,"%s\\out\\%s",sim_folder,tip_resume_file);}
			//else{sprintf(tip_path,"%s/out/%s",sim_folder,tip_resume_file);}
			fp=fopen(tip_resume_path,"r");
			if (fp==0 || fp==NULL)
			{
				// console
				printf("Could not locate %s\nThe program will now exit...\n",tip_path);
				// log
				fplog=fopen(logpath,"a");
				fprintf(fplog,"Could not locate %s\nThe program will now exit...\n",tip_path);
				fclose(fplog);
				exit(-100); // error code
			}
			else
			{
				// console
				printf("%s ... OK\n",tip_path);
				// log
				fplog=fopen(logpath,"a");
				fprintf(fplog,"%s ... OK\n",tip_path);
				fclose(fplog);
			}
			// Read tip coords and velocity
			fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&tip_particles);
			fcheck=fgets(buffer,cmax_length,fp);
			tip_x=(double*)malloc(tip_particles*sizeof(double));
			tip_y=(double*)malloc(tip_particles*sizeof(double));
			tip_z=(double*)malloc(tip_particles*sizeof(double));
			tip_species=(char**)malloc(tip_particles*sizeof(char*));
			for (i=0;i<tip_particles;++i)
			{
				tip_species[i]=(char*)malloc(species_length*sizeof(char));
			}
			for(i=0;i<tip_particles;++i)
			{
				fcheck=fgets(buffer,cmax_length,fp);
				sscanf(buffer,"%s\t%lf\t%lf\t%lf",tip_species[i],&tip_x[i],&tip_y[i],&tip_z[i]);
				if(tip_x[i]<xmin){printf("tip_x[%d]<xmin - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_x[i]>xmax){printf("tip_x[%d]>xmax - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_y[i]<ymin){printf("tip_y[%d]<ymin - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_y[i]>ymax){printf("tip_y[%d]>ymax - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_z[i]<zmin){printf("tip_z[%d]<zmin - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
				if(tip_z[i]>zmax){printf("tip_z[%d]>zmax - Redeploy the tip or change simulation box size\n",i+1);exit(-100);}
			
			}
			// Velocity components
			fcheck=fgets(buffer,cmax_length,fp);
			fcheck=fgets(buffer,cmax_length,fp);
			sscanf(buffer,"%lf\t%lf\t%lf",&tip_ux,&tip_uy,&tip_uz);
			// holder coords
			fcheck=fgets(buffer,cmax_length,fp);
			fcheck=fgets(buffer,cmax_length,fp);
			sscanf(buffer,"%lf\t%lf\t%lf",&xholder,&yholder,&zholder);
						
			fclose(fp);
		}
		
		// Preallocations
		tip_Z=(int*)malloc(tip_particles*sizeof(int));
		tip_sumFx=(double*)malloc(tip_particles*sizeof(double));
		tip_sumFy=(double*)malloc(tip_particles*sizeof(double));
		tip_sumFz=(double*)malloc(tip_particles*sizeof(double));
		tip_mass=(double*)malloc(tip_particles*sizeof(double));
		tip_atomicR=(double*)malloc(tip_particles*sizeof(double));
		tip_vanderWaalsR=(double*)malloc(tip_particles*sizeof(double));
		
		// Use the atoms_ID function to identify the tip atoms
		printf("Probe:\n");
		fplog=fopen(logpath,"a");fprintf(fplog,"Probe:\n");fclose(fplog);
		atoms_ID(tip_species,tip_particles,tip_Z,tip_mass,tip_atomicR,tip_vanderWaalsR,logpath);

		// calculate tip mass
		tip_M=0;
		for (i=0;i<tip_particles;++i){tip_M=tip_M+tip_mass[i];}
		
		// calculate CoM
		x_cm=0;y_cm=0;z_cm=0;
		for(i=0;i<tip_particles;++i)
		{
			x_cm=x_cm+tip_x[i]*tip_mass[i];
			y_cm=y_cm+tip_y[i]*tip_mass[i];
			z_cm=z_cm+tip_z[i]*tip_mass[i];
		}
		x_cm=x_cm/tip_M;
		y_cm=y_cm/tip_M;
		z_cm=z_cm/tip_M;
		
		// calculate holder coords
		if(resume==0)
		{
			xholder=x_cm;
			yholder=y_cm;
			zholder=z_cm;
		}
		
		// tip cells
		if(tip_linked_cells==1)
		{
			
			cells_init(tip_Mx,tip_My,tip_Mz,Lx,Ly,Lz,0,LJcutoff,tip_particles,tip_Z,logpath,&tip_Mx,&tip_My,&tip_Mz);
			tip_M3=tip_Mx*tip_My*tip_Mz;
			tip_Lcx=Lx/tip_Mx;tip_Lcy=Ly/tip_My;tip_Lcz=Lz/tip_Mz;
			// preallocations
			tip_cellpartarray=(int*)malloc(tip_particles*sizeof(int));// cells
			tip_neighborcellpartarray=(int*)malloc(tip_particles*sizeof(int));
			tip_head=(int*)malloc(tip_M3*sizeof(int));// head
			tip_list=(int*)malloc(tip_particles*sizeof(int));// list
			tip_neighbors=(int**)malloc(26*sizeof(int*));
			p1=(int*)malloc(26*tip_M3*sizeof(int));
			for(i=0;i<26;++i)
			{
				tip_neighbors[i]=p1+i*tip_M3;
			}
			build_neighbors(tip_Mx,tip_My,tip_Mz,tip_neighbors);
		}
			
	}
//(A END) --------------------------------------------------------------
	
	// --- Initial state calculations ----------------------------------

	n=0;																// "time" index
	
	if(tip==0 || (tip==1 && (tip_version==1 || tip_version==3)))
	{
		// create head and list arrays
		for (i=0;i<M3;++i){head[i]=-1;}
		for (i=0;i<particles;++i){list[i]=-1;}
		for (i=0;i<particles;++i)
		{
			cindex=floor((x_1[i]-xmin)/Lcx)+floor((y_1[i]-ymin)/Lcy)*Mx+floor((z_1[i]-zmin)/Lcz)*Mx*My;
			list[i]=head[cindex];
			head[cindex]=i;	
		}
		
		full=(int*)malloc(M3*sizeof(int));	
		fullN=-1;
		for(i=0;i<M3;++i)
		{
			if(head[i]!=-1)
			{
				fullN=fullN+1;
				full[fullN]=i;
			}
		}
		
		M3_0=fullN+1;
	
		cellsN=(int)floor(M3_0/size)+1;
		cells=(int*)malloc(cellsN*sizeof(int));	// for the master process
		book=(int**)malloc(size*sizeof(int*));
		p_book=(int*)malloc(size*cellsN*sizeof(int));
		for(i=0;i<size;++i)
		{
			book[i]=p_book+i*cellsN;
		}
		
		for(i=0;i<size;++i){for(j=0;j<cellsN;++j){book[i][j]=-1;}}
			
		i=0;l=0;
		while(i<M3_0)
		{
			for(k=0;k<size;++k)
			{
				if(i<M3_0)
				{
					book[k][l]=full[i];
					i=i+1;
				}
			}
			l=l+1;
		}
		
		// *** cells for the master process
		i=0;
		for(j=0;j<cellsN;++j){cells[j]=book[i][j];}
		// ***

					
		// send data
		for(i=1;i<size;++i)
		{		
			dest=i;
			id=7;
			MPI_Send(&x_1[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			id=8;
			MPI_Send(&y_1[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			id=9;
			MPI_Send(&z_1[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			id=10;
			MPI_Send(&head[0], M3, MPI_INT, dest, id, MPI_COMM_WORLD);
			id=11;
			MPI_Send(&list[0], particles, MPI_INT, dest, id, MPI_COMM_WORLD);

			id=90;
			MPI_Send(&cellsN, 1, MPI_INT, dest, id, MPI_COMM_WORLD);

			id=13;
			MPI_Send(&book[i][0], cellsN, MPI_INT, dest, id, MPI_COMM_WORLD);	
			
		}
		// master process calculations		
		if(interaction==1)
		{
			V_F_virial_cell_calc_tersoff_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,allpartarray,dx_zeta_III,dy_zeta_III,dz_zeta_III,cellsN,cells);
		}
		if(interaction==0)
		{
			V_F_virial_cell_calc_LJ_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,LJcutoff,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,cellsN,cells);
		}
		// receive and reduce to the final quantities
		for(i=1;i<size;++i)
		{
			source=i;
			id=14;
			MPI_Recv(&Vrecv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			V=V+Vrecv;
			
			id=20;
			MPI_Recv(&sigmaxx_recv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			id=21;
			MPI_Recv(&sigmayy_recv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			id=22;
			MPI_Recv(&sigmazz_recv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			
			sigmaxx=sigmaxx+sigmaxx_recv;
			sigmayy=sigmayy+sigmayy_recv;
			sigmazz=sigmazz+sigmazz_recv;
			
			id=16;
			MPI_Recv(&sumFxrecv[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			id=17;
			MPI_Recv(&sumFyrecv[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			id=18;
			MPI_Recv(&sumFzrecv[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			for(j=0;j<particles;++j)
			{
				sumFx[j]=sumFx[j]+sumFxrecv[j];
				sumFy[j]=sumFy[j]+sumFyrecv[j];
				sumFz[j]=sumFz[j]+sumFzrecv[j];
			}
		}	
	}
	else
	{
		V=0;
		for (i=0;i<particles;++i){sumFx[i]=0;sumFy[i]=0;sumFz[i]=0;}
	}

//(B) --- Tip-sample interaction ---------------------------------------
	if (tip==1)
	{
		if(tip_linked_cells==1)
		{

			// create head and list arrays for the substrate
			for (i=0;i<tip_M3;++i){head[i]=-1;}
			for (i=0;i<particles;++i){list[i]=-1;}
			for (i=0;i<particles;++i)
			{
				cindex=floor((x_1[i]-xmin)/tip_Lcx)+floor((y_1[i]-ymin)/tip_Lcy)*tip_Mx+floor((z_1[i]-zmin)/tip_Lcz)*tip_Mx*tip_My;
				list[i]=head[cindex];
				head[cindex]=i;	
			}
			
			// create head and list arrays for the probe
			for (i=0;i<tip_M3;++i){tip_head[i]=-1;}
			for (i=0;i<tip_particles;++i){tip_list[i]=-1;}
			for (i=0;i<tip_particles;++i)
			{

				cindex=floor((tip_x[i]-xmin)/tip_Lcx)+floor((tip_y[i]-ymin)/tip_Lcy)*tip_Mx+floor((tip_z[i]-zmin)/tip_Lcz)*tip_Mx*tip_My;
				tip_list[i]=tip_head[cindex];
				tip_head[cindex]=i;	
				
			}
			
			LJ_tip_cells(tip_particles,particles,tip_x,tip_y,tip_z,tip_Z,x_1,y_1,z_1,Z,Lx,Ly,Lz,LJcutoff,tip_sumFx,tip_sumFy,tip_sumFz,sumFx,sumFy,sumFz,&tip_V,tip_M3,tip_head,tip_list,head,list,tip_neighbors,tip_cellpartarray,tip_neighborcellpartarray,cellpartarray);
			
		}
		else
		{
			LJ_tip(tip_particles,particles,tip_x,tip_y,tip_z,tip_Z,x_1,y_1,z_1,Z,Lx,Ly,Lz,LJcutoff,tip_sumFx,tip_sumFy,tip_sumFz,sumFx,sumFy,sumFz,&tip_V);
		}
		
		V=V+tip_V;		// add to the potential energy the tip interaction energy
		tip_sampleFx=0;			// tip total force components
		tip_sampleFy=0;
		tip_sampleFz=0;
		
		for(i=0;i<tip_particles;++i)	// calculate tip force
		{
			tip_sampleFx=tip_sampleFx+tip_sumFx[i];
			tip_sampleFy=tip_sampleFy+tip_sumFy[i];
			tip_sampleFz=tip_sampleFz+tip_sumFz[i];
		}
			
		// export tip data	
		if(OS==1){sprintf(tip_path,"%s\\out\\%s",sim_folder,tip_data_file);}
		else{sprintf(tip_path,"%s/out/%s",sim_folder,tip_data_file);}
		fp=fopen(tip_path,"w+");
		//fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",n+nres,tip_sampleFx+Fxtarget,tip_sampleFy+Fytarget,tip_sampleFz+Fztarget,x_cm_1,y_cm_1,z_cm_1);
		fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
		n+nres,
		tip_sampleFx+Fxtarget-tip_gamma_x*tip_ux+kx*(xholder-x_cm),
		tip_sampleFy+Fytarget-tip_gamma_y*tip_uy+ky*(yholder-y_cm),
		tip_sampleFz+Fztarget-tip_gamma_z*tip_uz,
		x_cm,
		y_cm,
		z_cm);
		fclose(fp);
	}
	
//(B END) --------------------------------------------------------------
	
	// Accelerations, kinetic energy, temperature and virial pressure calculations
	//K=0;										// set to zero for the summation
	Kx=0;										// set to zero for the summation
	Ky=0;										// set to zero for the summation
	Kz=0;										// set to zero for the summation
	for (i=0;i<particles;++i)
	{
		ax_1[i]=sumFx[i]/mass[i];				// acceleration calculations
		ay_1[i]=sumFy[i]/mass[i];
		az_1[i]=sumFz[i]/mass[i];
		// fixed particles: set accelerations equal to zero
		if (mobileFlag[i]==0){ax_1[i]=0;ay_1[i]=0;az_1[i]=0;}
		// kinetic energy summation			
		Kx=Kx+0.5*mass[i]*ux_1[i]*ux_1[i];	
		Ky=Ky+0.5*mass[i]*uy_1[i]*uy_1[i];	
		Kz=Kz+0.5*mass[i]*uz_1[i]*uz_1[i];	
		//K=K+0.5*mass[i]*(ux_1[i]*ux_1[i]+uy_1[i]*uy_1[i]+uz_1[i]*uz_1[i]);	
	}
	K=Kx+Ky+Kz;
	Temperature=(2*K)/fpk;				// temperature
	
	Pxx=(2.0*Kx)/(Lx*Ly*Lz)+sigmaxx;
	Pyy=(2.0*Ky)/(Lx*Ly*Lz)+sigmayy;
	Pzz=(2.0*Kz)/(Lx*Ly*Lz)+sigmazz;
	//Pxx=sigmaxx;
	//Pyy=sigmayy;
	//Pzz=sigmazz;
	pressure=(Pxx+Pyy+Pzz)/3;
	
	// kinetic energy and temperature of the thermostated particles
	K_thermo=0;
	for(i=0;i<particles;++i){if (thermoFlag[i]==1){K_thermo=K_thermo+0.5*mass[i]*(ux_1[i]*ux_1[i]+uy_1[i]*uy_1[i]+uz_1[i]*uz_1[i]);}}
	Temperature_thermo=(2*K_thermo)/fpk_thermo;
	
	// Coordinates export
	////////////////////////////////////////////////////////////////////
	sprintf(buffer,"%s%s%d.csv",csvexportpath,csv_file,n+nres);
	fp=fopen(buffer,"w+");
	fprintf(fp,"X,Y,Z,radius,Flag,atomicZ\n");
	sprintf(floating,"%%%slf,%%%slf,%%%slf,%%.3lf,%%d,%%d\n",coords_acc,coords_acc,coords_acc);
	if(interaction==0)
	{
		for (i=0;i<particles;++i)
		{
			fprintf(fp,floating,x_1[i],y_1[i],z_1[i],vanderWaalsR[i],mobileFlag[i]+thermoFlag[i],Z[i]);
		}
		if(tip==1)
		{
			for (i=0;i<tip_particles;++i)
			{
				fprintf(fp,floating,tip_x[i],tip_y[i],tip_z[i],tip_vanderWaalsR[i],-1,tip_Z[i]);
			}
		}
	}
	if(interaction==1)
	{
		for (i=0;i<particles;++i)
		{
			fprintf(fp,floating,x_1[i],y_1[i],z_1[i],atomicR[i],mobileFlag[i]+thermoFlag[i],Z[i]);
		}
		if(tip==1)
		{
			for (i=0;i<tip_particles;++i)
			{
				fprintf(fp,floating,tip_x[i],tip_y[i],tip_z[i],tip_atomicR[i],-1,tip_Z[i]);
			}
		}
	}
	fclose(fp);
	////////////////////////////////////////////////////////////////////
	sprintf(buffer,"%s%s%d.xyz",xyzexportpath,xyz_file,n+nres);
	fp=fopen(buffer,"w+");
	if(tip==0){fprintf(fp,"%d\n%s\n",particles,title);}					// tip
	if(tip==1){fprintf(fp,"%d\n%s\n",particles+tip_particles,title);}	// tip
	sprintf(floating,"%%s %%%slf %%%slf %%%slf\n",coords_acc,coords_acc,coords_acc);
	for (i=0;i<particles;++i)
	{
		fprintf(fp,floating,species[i],x_1[i],y_1[i],z_1[i]);
	}
	if(tip==1)															// tip
	{
		for (i=0;i<tip_particles;++i)
		{
		fprintf(fp,floating,tip_species[i],tip_x[i],tip_y[i],tip_z[i]);
		}
	}
	fclose(fp);
	////////////////////////////////////////////////////////////////////
	// Jmol animation
	if (xyz_animation!=0)
	{
		//Jmol_animation(x_1,y_1,z_1,"w+");
		sprintf(buffer,"%s%s",xyzexportpath,animation_file);
		fp=fopen(buffer,"w+");
		if(tip==0){fprintf(fp,"%d\n%s\n",particles,title);}					// tip
		if(tip==1){fprintf(fp,"%d\n%s\n",particles+tip_particles,title);}	// tip
		sprintf(floating,"%%s %%%slf %%%slf %%%slf\n",coords_acc,coords_acc,coords_acc);
		for (i=0;i<particles;++i)
		{
			fprintf(fp,floating,species[i],x_1[i],y_1[i],z_1[i]);
		}
		if(tip==1)															// tip
		{
			for (i=0;i<tip_particles;++i)
			{
			fprintf(fp,floating,tip_species[i],tip_x[i],tip_y[i],tip_z[i]);
			}
		}
		fclose(fp);
	}	
	////////////////////////////////////////////////////////////////////
	// Velocities and positions export to dat files
	if (velout!=0)
	{
		sprintf(buffer,"%s%s%d.dat",velexportpath,vel_file,n+nres);
		fp=fopen(buffer,"w+");
		sprintf(floating,"%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc);
		for (i=0;i<particles;++i)
		{
			fprintf(fp,floating,ux_1[i],uy_1[i],uz_1[i]);
		}
		fclose(fp);
	}
	if (posout!=0)
	{
		sprintf(buffer,"%s%s%d.dat",posexportpath,pos_file,n+nres);
		fp=fopen(buffer,"w+");
		sprintf(floating,"%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc);
		for (i=0;i<particles;++i)
		{
			fprintf(fp,floating,x_1[i],y_1[i],z_1[i]);
		}
		fclose(fp);
	}	
	////////////////////////////////////////////////////////////////////
	
	// Create files in order to append data in them
	
	//-----
	
	fp=fopen(dataexportpath,"w+");
	//sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc,acc,acc,acc,acc);
	//fprintf(fp,floating,n+nres,K,V,Temperature,pressure*P_to_Mbar,sigmaxx*P_to_Mbar,sigmayy*P_to_Mbar,sigmazz*P_to_Mbar);
	sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc);
	fprintf(fp,floating,n+nres,K,V,Temperature);
	fclose(fp);
	
	//-----
	/*
	fp=fopen(thermoexportpath,"w+");
	sprintf(floating,"%%d\t%%%slf\t%%%slf\n",acc,acc);
	fprintf(fp,floating,n+nres,K_thermo,Temperature_thermo);
	fclose(fp);
	*/
	//-----
	
	//fp=fopen(CoM_path,"w+");
	//sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc,acc,acc,acc,acc,acc,acc);
	//sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc);
	
	//------------------------------------------------------------------
	// CoM position, velocity and angular velocity
	
	M=0;
	CoMx=0;
	CoMy=0;
	CoMz=0;
	
	//CoMux=0;
	//CoMuy=0;
	//CoMuz=0;
	
	//CoMLx=0;
	//CoMLy=0;
	//CoMLz=0;
	
	for(i=0;i<particles;++i)
	{
		if(mobileFlag[i]==1)
		{
			CoMx=CoMx+x_1[i]*mass[i];
			CoMy=CoMy+y_1[i]*mass[i];
			CoMz=CoMz+z_1[i]*mass[i];
			M=M+mass[i];
			
		//	CoMux=CoMux+ux_1[i];
		//	CoMuy=CoMuy+uy_1[i];
		//	CoMuz=CoMuz+uz_1[i];
			
		//	CoMLx=CoMLx+uz_1[i]*y_1[i]-uy_1[i]*z_1[i];
		//	CoMLy=CoMLy-uz_1[i]*x_1[i]+ux_1[i]*z_1[i];
		//	CoMLz=CoMLz+uy_1[i]*x_1[i]-ux_1[i]*y_1[i];
		}
	}
	CoMx=CoMx/M;
	CoMy=CoMy/M;
	CoMz=CoMz/M;
	
	//------------------------------------------------------------------
	
	//fprintf(fp,floating,n+nres,CoMx,CoMy,CoMz,CoMux,CoMuy,CoMuz,CoMLx,CoMLy,CoMLz);
	//fprintf(fp,floating,n+nres,CoMx,CoMy,CoMz);
	
	//fclose(fp);
	
	//-----
	
	if (VACFcalcstep>0  && VACFcalcstep<=myN)
	{
		fp=fopen(VACFexportpath,"w+");
		sprintf(floating,"%%d\t%%%slf\n",acc);
		fprintf(fp,floating,0,1.0);
		fclose(fp);
	}
	
	//-----
	
	// write pressure
	fp=fopen(P_path,"w+");
	sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc,acc);
	fprintf(fp,floating,n+nres,pressure*P_to_Mbar,Pxx*P_to_Mbar,Pyy*P_to_Mbar,Pzz*P_to_Mbar);
	fclose(fp);	
	
	//-----
	
	// --- Velocity-Verlet integration scheme --------------------------
	
	// Console output header
	if (console_out!=0)
	{
		// console
		printf("\nProgress:\n---units for K,V,E: eV\n---unit for T: Kelvin\n---unit for P: Mbar\n---unit for length: Angstrom\n\n");
		// log
		fplog=fopen(logpath,"a");
		fprintf(fplog,"\nProgress:\n---units for K,V,E: eV\n---unit for T: Kelvin\n---unit for P: Mbar\n---unit for length: Angstrom\n\n");
		fclose(fplog);
	}	
	
	// Time integration
	for (n=1;n<=myN;++n)												// time loop
	{
	
		// time integration
		if(tip==0 || (tip==1 && tip_version==1) || (tip==1 && tip_version==3))
		{
			for (i=0;i<particles;++i)
			{
				// Velocity-Verlet scheme
				x_1[i]=x_1[i]+dt*ux_1[i]+0.5*dt*dt*ax_1[i];
				y_1[i]=y_1[i]+dt*uy_1[i]+0.5*dt*dt*ay_1[i];
				z_1[i]=z_1[i]+dt*uz_1[i]+0.5*dt*dt*az_1[i];
			//	x_1[i]+=dt*ux_1[i]+0.5*dt*dt*ax_1[i];
			//	y_1[i]+=dt*uy_1[i]+0.5*dt*dt*ay_1[i];
			//	z_1[i]+=dt*uz_1[i]+0.5*dt*dt*az_1[i];
				if (x_1[i]>xmax){x_1[i]=x_1[i]-Lx;}						// PBCs
				if (x_1[i]<xmin){x_1[i]=x_1[i]+Lx;}
				if (y_1[i]>ymax){y_1[i]=y_1[i]-Ly;}
				if (y_1[i]<ymin){y_1[i]=y_1[i]+Ly;}
				if (z_1[i]>zmax){z_1[i]=z_1[i]-Lz;}
				if (z_1[i]<zmin){z_1[i]=z_1[i]+Lz;}
				
			}
		}
		
//(C) --- Tip movement -------------------------------------------------
		if(tip==1 && (tip_version==1 || tip_version==2))
		{
			x_cm=0;y_cm=0;z_cm=0;
			for(i=0;i<tip_particles;++i)
			{
				x_cm=x_cm+tip_x[i]*tip_mass[i];
				y_cm=y_cm+tip_y[i]*tip_mass[i];
				z_cm=z_cm+tip_z[i]*tip_mass[i];
			}
			x_cm=x_cm/tip_M;
			y_cm=y_cm/tip_M;
			z_cm=z_cm/tip_M;
			
			for(i=0;i<tip_particles;++i)
			{
				tip_x[i]=tip_x[i]+dt*tip_ux+0.5*dt*dt*(tip_sampleFx+Fxtarget-tip_gamma_x*tip_ux+kx*(xholder-x_cm))/tip_M;
				tip_y[i]=tip_y[i]+dt*tip_uy+0.5*dt*dt*(tip_sampleFy+Fytarget-tip_gamma_y*tip_uy+ky*(yholder-y_cm))/tip_M;
				tip_z[i]=tip_z[i]+dt*tip_uz+0.5*dt*dt*(tip_sampleFz+Fztarget-tip_gamma_z*tip_uz)/tip_M;
				
				if (tip_x[i]>xmax){tip_x[i]=tip_x[i]-Lx;}						// PBCs
				if (tip_x[i]<xmin){tip_x[i]=tip_x[i]+Lx;}
				if (tip_y[i]>ymax){tip_y[i]=tip_y[i]-Ly;}
				if (tip_y[i]<ymin){tip_y[i]=tip_y[i]+Ly;}
				if (tip_z[i]>zmax){tip_z[i]=tip_z[i]-Lz;}
				if (tip_z[i]<zmin){tip_z[i]=tip_z[i]+Lz;}
				
			}
			
			x_cm_1=0;y_cm_1=0;z_cm_1=0;
			for(i=0;i<tip_particles;++i)
			{
				x_cm_1=x_cm_1+tip_x[i]*tip_mass[i];
				y_cm_1=y_cm_1+tip_y[i]*tip_mass[i];
				z_cm_1=z_cm_1+tip_z[i]*tip_mass[i];
			}
			x_cm_1=x_cm_1/tip_M;
			y_cm_1=y_cm_1/tip_M;
			z_cm_1=z_cm_1/tip_M;
			
		}
//(C END) --------------------------------------------------------------

		// Potential energy, force and virial calculations
		if(tip==0 || (tip==1 && (tip_version==1 || tip_version==3)))
		{
			// create head and list arrays
			for (i=0;i<M3;++i){head[i]=-1;}
			for (i=0;i<particles;++i){list[i]=-1;}
			for (i=0;i<particles;++i)
			{
				cindex=floor((x_1[i]-xmin)/Lcx)+floor((y_1[i]-ymin)/Lcy)*Mx+floor((z_1[i]-zmin)/Lcz)*Mx*My;
				list[i]=head[cindex];
				head[cindex]=i;	
			}
			
			// repopulate book and send parts to all processes
			
			//full=(int*)malloc(M3*sizeof(int));	
			fullN=-1;
			for(i=0;i<M3;++i)
			{
				if(head[i]!=-1)
				{
					fullN=fullN+1;
					full[fullN]=i;
				}
			}
			
			M3_0=fullN+1;
		
			cellsN=(int)floor(M3_0/size)+1;
			
			/*
			book=(int**)malloc(size*sizeof(int*));
			p_book=(int*)malloc(size*cellsN*sizeof(int));
			for(i=0;i<size;++i)
			{
				book[i]=p_book+i*cellsN;
			}
			*/
			
			cells=realloc(cells,cellsN*sizeof(int));
			p_book=realloc(p_book,size*cellsN*sizeof(int));
			for(i=0;i<size;++i)
			{
				book[i]=p_book+i*cellsN;
			}
			
			
			for(i=0;i<size;++i){for(j=0;j<cellsN;++j){book[i][j]=-1;}}
				
			i=0;l=0;
			while(i<M3_0)
			{
				for(k=0;k<size;++k)
				{
					if(i<M3_0)
					{
						book[k][l]=full[i];
						i=i+1;
					}
				}
				l=l+1;
			}
			
			// *** cells for the master process
			i=0;
			for(j=0;j<cellsN;++j){cells[j]=book[i][j];}
			// ***
			
				
			// send data
			for(i=1;i<size;++i)
			{		
				dest=i;
				id=7;
				MPI_Send(&x_1[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				id=8;
				MPI_Send(&y_1[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				id=9;
				MPI_Send(&z_1[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				id=10;
				MPI_Send(&head[0], M3, MPI_INT, dest, id, MPI_COMM_WORLD);
				id=11;
				MPI_Send(&list[0], particles, MPI_INT, dest, id, MPI_COMM_WORLD);
				
				id=90;
				MPI_Send(&cellsN, 1, MPI_INT, dest, id, MPI_COMM_WORLD);

				id=13;
				MPI_Send(&book[i][0], cellsN, MPI_INT, dest, id, MPI_COMM_WORLD);	
				
			}
			// master process calculations		
			if(interaction==1)
			{
				V_F_virial_cell_calc_tersoff_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,allpartarray,dx_zeta_III,dy_zeta_III,dz_zeta_III,cellsN,cells);
			}
			if(interaction==0)
			{
				V_F_virial_cell_calc_LJ_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,LJcutoff,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,cellsN,cells);
			}
			// receive and reduce to the final quantities
			for(i=1;i<size;++i)
			{
				source=i;
				id=14;
				MPI_Recv(&Vrecv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				V=V+Vrecv;
				
				id=20;
				MPI_Recv(&sigmaxx_recv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=21;
				MPI_Recv(&sigmayy_recv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=22;
				MPI_Recv(&sigmazz_recv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				
				sigmaxx=sigmaxx+sigmaxx_recv;
				sigmayy=sigmayy+sigmayy_recv;
				sigmazz=sigmazz+sigmazz_recv;
					
				id=16;
				MPI_Recv(&sumFxrecv[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=17;
				MPI_Recv(&sumFyrecv[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=18;
				MPI_Recv(&sumFzrecv[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				for(j=0;j<particles;++j)
				{
					sumFx[j]=sumFx[j]+sumFxrecv[j];
					sumFy[j]=sumFy[j]+sumFyrecv[j];
					sumFz[j]=sumFz[j]+sumFzrecv[j];
				}
			}

		}
		else
		{
			V=0;
		}
		
//(D) --- Tip-sample interaction ---------------------------------------
		if (tip==1)
		{
			// interaction
			if(tip_linked_cells==1)
			{
				// create head and list arrays for the substrate
				for (i=0;i<tip_M3;++i){head[i]=-1;}
				for (i=0;i<particles;++i){list[i]=-1;}
				for (i=0;i<particles;++i)
				{
					cindex=floor((x_1[i]-xmin)/tip_Lcx)+floor((y_1[i]-ymin)/tip_Lcy)*tip_Mx+floor((z_1[i]-zmin)/tip_Lcz)*tip_Mx*tip_My;
					list[i]=head[cindex];
					head[cindex]=i;	
				}
				
				// create head and list arrays for the probe
				for (i=0;i<tip_M3;++i){tip_head[i]=-1;}
				for (i=0;i<tip_particles;++i){tip_list[i]=-1;}
				for (i=0;i<tip_particles;++i)
				{
					cindex=floor((tip_x[i]-xmin)/tip_Lcx)+floor((tip_y[i]-ymin)/tip_Lcy)*tip_Mx+floor((tip_z[i]-zmin)/tip_Lcz)*tip_Mx*tip_My;
					tip_list[i]=tip_head[cindex];
					tip_head[cindex]=i;	
				}
				
				LJ_tip_cells(tip_particles,particles,tip_x,tip_y,tip_z,tip_Z,x_1,y_1,z_1,Z,Lx,Ly,Lz,LJcutoff,tip_sumFx,tip_sumFy,tip_sumFz,sumFx,sumFy,sumFz,&tip_V,tip_M3,tip_head,tip_list,head,list,tip_neighbors,tip_cellpartarray,tip_neighborcellpartarray,cellpartarray);
				
			}
			else
			{
				LJ_tip(tip_particles,particles,tip_x,tip_y,tip_z,tip_Z,x_1,y_1,z_1,Z,Lx,Ly,Lz,LJcutoff,tip_sumFx,tip_sumFy,tip_sumFz,sumFx,sumFy,sumFz,&tip_V);
			}
			
			V=V+tip_V;
			// forces
			tip_sampleFx_1=0;
			tip_sampleFy_1=0;
			tip_sampleFz_1=0;
			for(i=0;i<tip_particles;++i)
			{
				tip_sampleFx_1=tip_sampleFx_1+tip_sumFx[i];
				tip_sampleFy_1=tip_sampleFy_1+tip_sumFy[i];
				tip_sampleFz_1=tip_sampleFz_1+tip_sumFz[i];
			}
			// tip velocity
			tip_ux=(tip_ux+0.5*dt*(tip_sampleFx+tip_sampleFx_1+2*Fxtarget-tip_gamma_x*tip_ux+kx*(xholder-x_cm)+kx*(xholder-x_cm_1))/tip_M)/(1.0+0.5*dt*tip_gamma_x/tip_M);
			tip_uy=(tip_uy+0.5*dt*(tip_sampleFy+tip_sampleFy_1+2*Fytarget-tip_gamma_y*tip_uy+ky*(yholder-y_cm)+ky*(yholder-y_cm_1))/tip_M)/(1.0+0.5*dt*tip_gamma_y/tip_M);
			tip_uz=(tip_uz+0.5*dt*(tip_sampleFz+tip_sampleFz_1+2*Fztarget-tip_gamma_z*tip_uz)/tip_M)/(1.0+0.5*dt*tip_gamma_z/tip_M);
			if(tip_version==3)	// rigid tip: zero velocity
			{
				tip_ux=0;tip_uy=0;tip_uz=0;
			}
			//update
			tip_sampleFz=tip_sampleFz_1;
			tip_sampleFy=tip_sampleFy_1;
			tip_sampleFx=tip_sampleFx_1;
			
			// export tip data
			//if(OS==1){sprintf(tip_path,"%s\\out\\%s",sim_folder,tip_data_file);}
			//else{sprintf(tip_path,"%s/out/%s",sim_folder,tip_data_file);}
			if(n%steps==0)
			{
				fp=fopen(tip_path,"a");
				//fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",n+nres,tip_sampleFx+Fxtarget,tip_sampleFy+Fytarget,tip_sampleFz+Fztarget,x_cm_1,y_cm_1,z_cm_1);
				fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
				n+nres,
				tip_sampleFx_1+Fxtarget-tip_gamma_x*tip_ux+kx*(xholder-x_cm_1),
				tip_sampleFy_1+Fytarget-tip_gamma_y*tip_uy+ky*(yholder-y_cm_1),
				tip_sampleFz_1+Fztarget-tip_gamma_z*tip_uz,
				x_cm_1,
				y_cm_1,
				z_cm_1);
				fclose(fp);
			}
			
		}

//(D END) --------------------------------------------------------------
		// rigid system
		if(tip==1&&tip_version==2)
		{
			for(i=0;i<particles;++i)
			{
				ux_1[i]=0;uy_1[i]=0;uz_1[i]=0;
				ax_2[i]=0;ay_2[i]=0;az_2[i]=0;
			}
			K_thermo=0;
			K=0;
		}
		// moving system
		else
		{
			// Thermostat check
			if (thermostat[n]==-1)											// NVE
			{
				// total force loop and calculations
				//K=0;
				Kx=0;
				Ky=0;
				Kz=0;
				K_thermo=0;										// set to zero for the summation
				for (i=0;i<particles;++i)
				{
					ax_2[i]=sumFx[i]/mass[i];								// acceleration calculation
					ay_2[i]=sumFy[i]/mass[i];
					az_2[i]=sumFz[i]/mass[i];
					if (mobileFlag[i]==0){ax_2[i]=0;ay_2[i]=0;az_2[i]=0;}	// fixed particles: set accelerations equal to zero
					ux_1[i]=ux_1[i]+0.5*dt*ax_1[i]+0.5*dt*ax_2[i];			// velocity calculations
					uy_1[i]=uy_1[i]+0.5*dt*ay_1[i]+0.5*dt*ay_2[i];
					uz_1[i]=uz_1[i]+0.5*dt*az_1[i]+0.5*dt*az_2[i];
				//	ux_1[i]+=0.5*dt*ax_1[i]+0.5*dt*ax_2[i];			// velocity calculations
				//	uy_1[i]+=0.5*dt*ay_1[i]+0.5*dt*ay_2[i];
				//	uz_1[i]+=0.5*dt*az_1[i]+0.5*dt*az_2[i];
					// kinetic energy summation
					//K=K+0.5*mass[i]*(ux_1[i]*ux_1[i]+uy_1[i]*uy_1[i]+uz_1[i]*uz_1[i]);
					Kx=Kx+0.5*mass[i]*ux_1[i]*ux_1[i];
					Ky=Ky+0.5*mass[i]*uy_1[i]*uy_1[i];
					Kz=Kz+0.5*mass[i]*uz_1[i]*uz_1[i];
					// kinetic energy of the thermostated particles
					if (thermoFlag[i]==1){K_thermo=K_thermo+0.5*mass[i]*(ux_1[i]*ux_1[i]+uy_1[i]*uy_1[i]+uz_1[i]*uz_1[i]);}
				}
				K=Kx+Ky+Kz;
			}
			else															// NVT
			{
				// total force loop and calculations
				for (i=0;i<particles;++i)
				{
					ax_2[i]=sumFx[i]/mass[i];								// acceleration calculation
					ay_2[i]=sumFy[i]/mass[i];
					az_2[i]=sumFz[i]/mass[i];
					if (mobileFlag[i]==0){ax_2[i]=0;ay_2[i]=0;az_2[i]=0;}	// fixed particles: set accelerations equal to zero
					ux_1[i]=ux_1[i]+0.5*dt*ax_1[i]+0.5*dt*ax_2[i];			// velocity calculations
					uy_1[i]=uy_1[i]+0.5*dt*ay_1[i]+0.5*dt*ay_2[i];
					uz_1[i]=uz_1[i]+0.5*dt*az_1[i]+0.5*dt*az_2[i];	
				//	ux_1[i]+=0.5*dt*ax_1[i]+0.5*dt*ax_2[i];			// velocity calculations
				//	uy_1[i]+=0.5*dt*ay_1[i]+0.5*dt*ay_2[i];
				//	uz_1[i]+=0.5*dt*az_1[i]+0.5*dt*az_2[i];	
				}
				// thermostat
				// **Andersen thermostat
				if (berendsen_t==0)											
				{
					int_random=rand();
					myrand=(double)int_random/RAND_MAX;
					// thermostat application
					if (myrand<contactP)
					{
						v_rescale(ux_1,uy_1,uz_1,thermostat[n],particleP,particles,thermoFlag,thermo_particles,mass);
					}	
				}
				// **Berendsen thermostat
				else											
				{
					// kinetic energy calculation - set to zero for the summation
					K_thermo=0;				
					for (i=0;i<particles;++i)
					{
						// kinetic energy summation
						if (thermoFlag[i]==1){K_thermo=K_thermo+0.5*mass[i]*(ux_1[i]*ux_1[i]+uy_1[i]*uy_1[i]+uz_1[i]*uz_1[i]);}
					}
					Temperature_thermo=(2*K_thermo)/fpk_thermo;
					lamda=sqrt(1.0+(dtSI/berendsen_t)*(thermostat[n]/Temperature_thermo-1.0));
					for (i=0;i<particles;++i)
					{
						if (thermoFlag[i]==1)
						{
							ux_1[i]=lamda*ux_1[i];
							uy_1[i]=lamda*uy_1[i];
							uz_1[i]=lamda*uz_1[i];
						}
					}
				}
				
				// kinetic energy calculation - set to zero for the summation
				//K=0;						
				Kx=0;						
				Ky=0;						
				Kz=0;						
				for (i=0;i<particles;++i)
				{
					//K=K+0.5*mass[i]*(ux_1[i]*ux_1[i]+uy_1[i]*uy_1[i]+uz_1[i]*uz_1[i]);
					Kx=Kx+0.5*mass[i]*ux_1[i]*ux_1[i];
					Ky=Ky+0.5*mass[i]*uy_1[i]*uy_1[i];
					Kz=Kz+0.5*mass[i]*uz_1[i]*uz_1[i];
				}
				K=Kx+Ky+Kz;
				// kinetic energy and temperature of the thermostated particles
				K_thermo=0;				
				for(i=0;i<particles;++i)
				{
					if (thermoFlag[i]==1){K_thermo=K_thermo+0.5*mass[i]*(ux_1[i]*ux_1[i]+uy_1[i]*uy_1[i]+uz_1[i]*uz_1[i]);}
				}

			}
		}
		
		// Virial pressure and temperature
		Temperature=(2*K)/fpk;
		Temperature_thermo=(2*K_thermo)/fpk_thermo;
		
		Pxx=(2.0*Kx)/(Lx*Ly*Lz)+sigmaxx;
		Pyy=(2.0*Ky)/(Lx*Ly*Lz)+sigmayy;
		Pzz=(2.0*Kz)/(Lx*Ly*Lz)+sigmazz;
		//Pxx=sigmaxx;
		//Pyy=sigmayy;
		//Pzz=sigmazz;
		pressure=(Pxx+Pyy+Pzz)/3;
		
		// -- A section --
		
		change_box=0;
		mixx=1.0;
		miyy=1.0;
		mizz=1.0;
		if(bflagX[n]==1)
		{
			mi_factor_2=1.0+mi_factor*(Pxx-barostatX[n]);
			mixx=pow(mi_factor_2,1.0/3.0);
			xmin=mixx*xmin;
			xmax=mixx*xmax;
			Lx=mixx*Lx;
			change_box=1;
		}
		if(bflagY[n]==1)
		{
			mi_factor_2=1.0+mi_factor*(Pyy-barostatY[n]);
			miyy=pow(mi_factor_2,1.0/3.0);
			ymin=miyy*ymin;
			ymax=miyy*ymax;
			Ly=miyy*Ly;
			change_box=1;
		}
		if(bflagZ[n]==1)
		{
			mi_factor_2=1.0+mi_factor*(Pzz-barostatZ[n]);
			mizz=pow(mi_factor_2,1.0/3.0);
			zmin=mizz*zmin;
			zmax=mizz*zmax;
			Lz=mizz*Lz;
			change_box=1;
		}
		
		// communicate: change_box
		for(i=1;i<size;++i)
		{		
			dest=i;
			id=40;
			MPI_Send(&change_box, 1, MPI_INT, dest, id, MPI_COMM_WORLD);			
		}
		
		if(change_box==1)
		{
			
			// communicate mi
			for(i=1;i<size;++i)
			{		
				dest=i;
				id=41;
				MPI_Send(&mixx, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);			
				id=42;
				MPI_Send(&miyy, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);			
				id=43;
				MPI_Send(&mizz, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);	
				id=4;
				MPI_Send(&Lx, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				id=5;
				MPI_Send(&Ly, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				id=6;
				MPI_Send(&Lz, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);		
			}
			
			for(i=0;i<particles;++i)
			{
				x_1[i]=mixx*x_1[i];
				y_1[i]=miyy*y_1[i];
				z_1[i]=mizz*z_1[i];
			}
			
			// cells	
			Mx_new=floor(Lx/cut);
			My_new=floor(Ly/cut);
			Mz_new=floor(Lz/cut);
			if(Mx_new<3){Mx_new=3;}
			if(My_new<3){My_new=3;}
			if(Mz_new<3){Mz_new=3;}
			
			if(Mx_new!=Mx || My_new!=My || Mz_new!=Mz)
			{
				Mx=Mx_new;
				My=My_new;
				Mz=Mz_new;
				M3=Mx*My*Mz;
				Lcx=Lx/Mx;Lcy=Ly/My;Lcz=Lz/Mz;
				head=realloc(head,M3*sizeof(int));
				p=realloc(p,26*M3*sizeof(int));
				for(i=0;i<26;++i)
				{
					neighbors[i]=p+i*M3;
				}
				build_neighbors(Mx,My,Mz,neighbors);
				/*
				cellsN=(int)floor(M3/size)+1;
				cells=realloc(cells,cellsN*sizeof(int));
				p_book=realloc(p_book,size*cellsN*sizeof(int));
				for(i=0;i<size;++i)
				{
					book[i]=p_book+i*cellsN;
				}
				i=0;l=0;
				while(i<M3)
				{
					for(k=0;k<size;++k)
					{
						book[k][cellsN-1]=-1;	// <-----
						if(i<M3)
						{
							book[k][l]=i;
							i=i+1;
						}
					}
					l=l+1;
				}
				i=0;
				for(j=0;j<cellsN;++j){cells[j]=book[i][j];}
				*/
				for (i=0;i<M3;++i){head[i]=-1;}
				for (i=0;i<particles;++i){list[i]=-1;}
				for (i=0;i<particles;++i)
				{
					cindex=floor((x_1[i]-xmin)/Lcx)+floor((y_1[i]-ymin)/Lcy)*Mx+floor((z_1[i]-zmin)/Lcz)*Mx*My;
					list[i]=head[cindex];
					head[cindex]=i;	
				}
				full=realloc(full,M3*sizeof(int));	
				fullN=-1;
				for(i=0;i<M3;++i)
				{
					if(head[i]!=-1)
					{
						fullN=fullN+1;
						full[fullN]=i;
					}
				}
				
				M3_0=fullN+1;
			
				cellsN=(int)floor(M3_0/size)+1;
				
				cells=realloc(cells,cellsN*sizeof(int));
				p_book=realloc(p_book,size*cellsN*sizeof(int));
				for(i=0;i<size;++i)
				{
					book[i]=p_book+i*cellsN;
				}
				
				for(i=0;i<size;++i){for(j=0;j<cellsN;++j){book[i][j]=-1;}}
				
				i=0;l=0;
				while(i<M3_0)
				{
					for(k=0;k<size;++k)
					{
						if(i<M3_0)
						{
							book[k][l]=full[i];
							i=i+1;
						}
					}
					l=l+1;
				}
				
				i=0;
				for(j=0;j<cellsN;++j){cells[j]=book[i][j];}
				
				realloc_flag=1;
				for(i=1;i<size;++i)
				{		
					dest=i;
					id=44;
					MPI_Send(&realloc_flag, 1, MPI_INT, dest, id, MPI_COMM_WORLD);			
		//		}
		//		for(i=1;i<size;++i)
		//		{		
					dest=i;
					id=45;
					MPI_Send(&M3, 1, MPI_INT, dest, id, MPI_COMM_WORLD);
					id=46;
					MPI_Send(&neighbors[0][0], 26*M3, MPI_INT, dest, id, MPI_COMM_WORLD);
					
					id=90;
					MPI_Send(&cellsN, 1, MPI_INT, dest, id, MPI_COMM_WORLD);
					
					id=47;
					MPI_Send(&book[i][0], cellsN, MPI_INT, dest, id, MPI_COMM_WORLD);		
					id=48;
					MPI_Send(&head[0], M3, MPI_INT, dest, id, MPI_COMM_WORLD);
					id=49;
					MPI_Send(&list[0], particles, MPI_INT, dest, id, MPI_COMM_WORLD);				
				}
				// communicate: M3, neighbors, book, head, list
				
			}
			else
			{
				Lcx=Lx/Mx;Lcy=Ly/My;Lcz=Lz/Mz;
				for (i=0;i<M3;++i){head[i]=-1;}
				for (i=0;i<particles;++i){list[i]=-1;}
				for (i=0;i<particles;++i)
				{
					cindex=floor((x_1[i]-xmin)/Lcx)+floor((y_1[i]-ymin)/Lcy)*Mx+floor((z_1[i]-zmin)/Lcz)*Mx*My;
					list[i]=head[cindex];
					head[cindex]=i;	
				}
				
				full=realloc(full,M3*sizeof(int));	
				fullN=-1;
				for(i=0;i<M3;++i)
				{
					if(head[i]!=-1)
					{
						fullN=fullN+1;
						full[fullN]=i;
					}
				}
				
				M3_0=fullN+1;
			
				cellsN=(int)floor(M3_0/size)+1;
				
				cells=realloc(cells,cellsN*sizeof(int));
				p_book=realloc(p_book,size*cellsN*sizeof(int));
				for(i=0;i<size;++i)
				{
					book[i]=p_book+i*cellsN;
				}
				
				for(i=0;i<size;++i){for(j=0;j<cellsN;++j){book[i][j]=-1;}}
				
				i=0;l=0;
				while(i<M3_0)
				{
					for(k=0;k<size;++k)
					{
						if(i<M3_0)
						{
							book[k][l]=full[i];
							i=i+1;
						}
					}
					l=l+1;
				}
				
				i=0;
				for(j=0;j<cellsN;++j){cells[j]=book[i][j];}
				
				realloc_flag=0;
				for(i=1;i<size;++i)
				{		
					dest=i;
					id=44;
					MPI_Send(&realloc_flag, 1, MPI_INT, dest, id, MPI_COMM_WORLD);			
			//	}
			//	for(i=1;i<size;++i)
			//	{				
					id=60;
					MPI_Send(&head[0], M3, MPI_INT, dest, id, MPI_COMM_WORLD);
					id=61;
					MPI_Send(&list[0], particles, MPI_INT, dest, id, MPI_COMM_WORLD);	
					
					id=90;
					MPI_Send(&cellsN, 1, MPI_INT, dest, id, MPI_COMM_WORLD);
					
					id=47;
					MPI_Send(&book[i][0], cellsN, MPI_INT, dest, id, MPI_COMM_WORLD);
								
				}
				// communicate: head, list
				
			}
			
			// MASTER calculations
			// receive from SLAVE nodes and reduce (V,F,sigma)
			
			// master process calculations		
			if(interaction==1)
			{
				V_F_virial_cell_calc_tersoff_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,allpartarray,dx_zeta_III,dy_zeta_III,dz_zeta_III,cellsN,cells);
			}
			if(interaction==0)
			{
				V_F_virial_cell_calc_LJ_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,LJcutoff,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,cellsN,cells);
			}
			
			// receive and reduce to the final quantities
			for(i=1;i<size;++i)
			{				
				source=i;
				id=50;
				MPI_Recv(&Vrecv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				V=V+Vrecv;
								
				id=20;
				MPI_Recv(&sigmaxx_recv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=21;
				MPI_Recv(&sigmayy_recv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=22;
				MPI_Recv(&sigmazz_recv, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				
				sigmaxx=sigmaxx+sigmaxx_recv;
				sigmayy=sigmayy+sigmayy_recv;
				sigmazz=sigmazz+sigmazz_recv;
	
				id=16;
				MPI_Recv(&sumFxrecv[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=17;
				MPI_Recv(&sumFyrecv[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=18;
				MPI_Recv(&sumFzrecv[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				for(j=0;j<particles;++j)
				{
					sumFx[j]=sumFx[j]+sumFxrecv[j];
					sumFy[j]=sumFy[j]+sumFyrecv[j];
					sumFz[j]=sumFz[j]+sumFzrecv[j];
				}
			}
						
			// calculate a, P
			for (i=0;i<particles;++i)
			{
				ax_2[i]=sumFx[i]/mass[i];								// acceleration calculation
				ay_2[i]=sumFy[i]/mass[i];
				az_2[i]=sumFz[i]/mass[i];
				if (mobileFlag[i]==0){ax_2[i]=0;ay_2[i]=0;az_2[i]=0;}	
			}
			Pxx=(2.0*Kx)/(Lx*Ly*Lz)+sigmaxx;
			Pyy=(2.0*Ky)/(Lx*Ly*Lz)+sigmayy;
			Pzz=(2.0*Kz)/(Lx*Ly*Lz)+sigmazz;
			
			pressure=(Pxx+Pyy+Pzz)/3;
		}

		// -- end of A section --

		// CoM position, velocity and angular velocity initializations
		M=0;
		CoMx=0;
		CoMy=0;
		CoMz=0;
				
		//CoMux=0;
		//CoMuy=0;
		//CoMuz=0;
		
		//CoMLx=0;
		//CoMLy=0;
		//CoMLz=0;
		
		// Array update
		for (i=0;i<particles;++i)
		{
			ax_1[i]=ax_2[i];ay_1[i]=ay_2[i];az_1[i]=az_2[i];
			
			//------------------------------------------------------------------
			// CoM 
			
			if(mobileFlag[i]==1)
			{
				CoMx=CoMx+x_1[i]*mass[i];
				CoMy=CoMy+y_1[i]*mass[i];
				CoMz=CoMz+z_1[i]*mass[i];
				M=M+mass[i];
				
				
				//CoMux=CoMux+ux_1[i];
				//CoMuy=CoMuy+uy_1[i];
				//CoMuz=CoMuz+uz_1[i];
				
				//CoMLx=CoMLx+uz_1[i]*y_1[i]-uy_1[i]*z_1[i];
				//CoMLy=CoMLy-uz_1[i]*x_1[i]+ux_1[i]*z_1[i];
				//CoMLz=CoMLz+uy_1[i]*x_1[i]-ux_1[i]*y_1[i];
				
			}
			
			//------------------------------------------------------------------
			
		}
		CoMx=CoMx/M;
		CoMy=CoMy/M;
		CoMz=CoMz/M;
		
		// VACF calculation
		// check time step to save the starting point velocities
		if (n==myN-VACFcalcstep+1 && VACFcalcstep>0 && VACFcalcstep<=myN)	
		{
			VACFinit=0;
			// loop over particles
			for (i=0;i<particles;++i)					
			{
				// check if the particle is mobile or fixed
				if (mobileFlag[i]==1)					
				{
					// save the reference velocity
					uxsave[i]=ux_1[i];uysave[i]=uy_1[i];uzsave[i]=uz_1[i];		
					VACFinit=VACFinit+uxsave[i]*uxsave[i]+uysave[i]*uysave[i]+uzsave[i]*uzsave[i];
				}
			}
			// set to 1 for normalized VACF
			//VACF[n-myN+VACFcalcstep-1]=1;							
		}
		// check time step to calculate VACF
		if (n>myN-VACFcalcstep+1 && VACFcalcstep>0 && VACFcalcstep<=myN)			
		{
			VACFsum=0;
			// sum over all particles								
			for (i=0;i<particles;++i)
			{
				// check if the particle is mobile or fixed
				if (mobileFlag[i]==1)							
				{
					VACFsum=VACFsum+uxsave[i]*ux_1[i]+uysave[i]*uy_1[i]+uzsave[i]*uz_1[i];
				}
			}
			//VACF[n-myN+VACFcalcstep-1]=VACFsum/VACFinit;
			VACF=VACFsum/VACFinit;
			
			fp=fopen(VACFexportpath,"a");
			sprintf(floating,"%%d\t%%%slf\n",acc);
			fprintf(fp,floating,n-myN+VACFcalcstep-1,VACF);
			fclose(fp);
			
		}
		
		// write data
		//-----
		if (n % steps == 0)
		{
			fp=fopen(dataexportpath,"a");
			//sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc,acc,acc,acc,acc);
			//fprintf(fp,floating,n+nres,K,V,Temperature,pressure*P_to_Mbar,sigmaxx*P_to_Mbar,sigmayy*P_to_Mbar,sigmazz*P_to_Mbar);
			sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc);
			fprintf(fp,floating,n+nres,K,V,Temperature);
			fclose(fp);
			/*
			fp=fopen(thermoexportpath,"a");
			sprintf(floating,"%%d\t%%%slf\t%%%slf\n",acc,acc);
			fprintf(fp,floating,n+nres,K_thermo,Temperature_thermo);
			fclose(fp);
			*/
			//fp=fopen(CoM_path,"a");
			//sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc,acc,acc,acc,acc,acc,acc);
			//fprintf(fp,floating,n+nres,CoMx,CoMy,CoMz,CoMux,CoMuy,CoMuz,CoMLx,CoMLy,CoMLz);
			//sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc);
			//fprintf(fp,floating,n+nres,CoMx,CoMy,CoMz);
		
			//fclose(fp);
			
			// write pressure
			fp=fopen(P_path,"a");
			sprintf(floating,"%%d\t%%%slf\t%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc,acc);
			fprintf(fp,floating,n+nres,pressure*P_to_Mbar,Pxx*P_to_Mbar,Pyy*P_to_Mbar,Pzz*P_to_Mbar);
			fclose(fp);	
			
		}
		//-----
		
		// Console output: progress monitoring
		if (console_out!=0)
		{
			
			fplog=fopen(logpath,"a");
				
			if (n %  console_out== 0)
			{
				if(thermostat[n]==-1){
					
					if(bflagX[n]==0 && bflagY[n]==0 && bflagZ[n]==0)
					{					
						// console
						printf("[NVE] ");
						// log
						fprintf(fplog,"[NVE] ");
					}
					else
					{
						// console
						printf("[Ber_P] ");
						// log
						fprintf(fplog,"[Ber_P] ");
					}
					}
				else if(berendsen_t==0){
					if(bflagX[n]==0 && bflagY[n]==0 && bflagZ[n]==0)
					{
						// console
						printf("[NVT-Andersen] ");
						// log
						fprintf(fplog,"[NVT-Andersen] ");
					}
					else
					{
						// console
						printf("[NPT: Ander_T-Ber_P] ");
						// log
						fprintf(fplog,"[NPT: Ander_T-Ber_P] ");
					}
					}
				else if(berendsen_t!=0){
					if(bflagX[n]==0 && bflagY[n]==0 && bflagZ[n]==0)
					{
						// console
						printf("[NVT-Berendsen] ");
						// log
						fprintf(fplog,"[NVT-Berendsen] ");
					}
					else
					{
						// console
						printf("[NPT] ");
						// log
						fprintf(fplog,"[NPT] ");
					}
					}
					// console
					printf("(%d) %d steps out of %d. ",n+nres,n,myN);
					// log
					fprintf(fplog,"(%d) %d steps out of %d. ",n+nres,n,myN);
					// console
					printf("K=%.4lf -- V=%.4lf -- T=%.4lf -- E=%.4lf\n",K,V,Temperature,K+V);
					// log
					fprintf(fplog,"K=%.4lf -- V=%.4lf -- T=%.4lf -- E=%.4lf\n",K,V,Temperature,K+V);
					//
					if(bflagX[n]!=0 || bflagY[n]!=0 || bflagZ[n]!=0)
					{
					// console
					printf("(%dx%dx%d) -- Pxx=%.4lf -- Pyy=%.4lf -- Pzz=%.4lf -- Lx=%.4lf -- Ly=%.4lf -- Lz=%.4lf\n",Mx,My,Mz,Pxx*P_to_Mbar,Pyy*P_to_Mbar,Pzz*P_to_Mbar,Lx,Ly,Lz);
					// log
					fprintf(fplog,"(%dx%dx%d) -- Pxx=%.4lf -- Pyy=%.4lf -- Pzz=%.4lf -- Lx=%.4lf -- Ly=%.4lf -- Lz=%.4lf\n",Mx,My,Mz,Pxx*P_to_Mbar,Pyy*P_to_Mbar,Pzz*P_to_Mbar,Lx,Ly,Lz);
						
					
					}
			}
			
		fclose(fplog);	
			
		}		
		
		// Coordinates export
		////////////////////////////////////////////////////////////////
		if (csvout>0 && (n % csvout==0))
		{
			sprintf(buffer,"%s%s%d.csv",csvexportpath,csv_file,n+nres);
			fp=fopen(buffer,"w+");
			fprintf(fp,"X,Y,Z,radius,Flag,atomicZ\n");
			sprintf(floating,"%%%slf,%%%slf,%%%slf,%%.3lf,%%d,%%d\n",coords_acc,coords_acc,coords_acc);
			if(interaction==0)
			{
				for (i=0;i<particles;++i)
				{
					fprintf(fp,floating,x_1[i],y_1[i],z_1[i],vanderWaalsR[i],mobileFlag[i]+thermoFlag[i],Z[i]);
				}
				if(tip==1)
				{
					for (i=0;i<tip_particles;++i)
					{
						fprintf(fp,floating,tip_x[i],tip_y[i],tip_z[i],tip_vanderWaalsR[i],-1,tip_Z[i]);
					}
				}
			}
			if(interaction==1)
			{
				for (i=0;i<particles;++i)
				{
					fprintf(fp,floating,x_1[i],y_1[i],z_1[i],atomicR[i],mobileFlag[i]+thermoFlag[i],Z[i]);
				}
				if(tip==1)
				{
					for (i=0;i<tip_particles;++i)
					{
						fprintf(fp,floating,tip_x[i],tip_y[i],tip_z[i],tip_atomicR[i],-1,tip_Z[i]);
					}
				}
			}
			fclose(fp);
		}
		////////////////////////////////////////////////////////////////
		if (xyzout>0 && (n % xyzout==0))
		{
			sprintf(buffer,"%s%s%d.xyz",xyzexportpath,xyz_file,n+nres);
			fp=fopen(buffer,"w+");
			if(tip==0){fprintf(fp,"%d\n%s\n",particles,title);}					// tip
			if(tip==1){fprintf(fp,"%d\n%s\n",particles+tip_particles,title);}	// tip
			sprintf(floating,"%%s %%%slf %%%slf %%%slf\n",coords_acc,coords_acc,coords_acc);
			for (i=0;i<particles;++i)
			{
				fprintf(fp,floating,species[i],x_1[i],y_1[i],z_1[i]);
			}
			if(tip==1)															// tip
			{
				for (i=0;i<tip_particles;++i)
				{
				fprintf(fp,floating,tip_species[i],tip_x[i],tip_y[i],tip_z[i]);
				}
			}
			fclose(fp);
		}
		////////////////////////////////////////////////////////////////
		// Jmol animation
		if(xyz_animation!=0 && (n % xyz_animation==0))
		{
			sprintf(buffer,"%s%s",xyzexportpath,animation_file);
			fp=fopen(buffer,"a");
			if(tip==0){fprintf(fp,"%d\n%s\n",particles,title);}					// tip
			if(tip==1){fprintf(fp,"%d\n%s\n",particles+tip_particles,title);}	// tip
			sprintf(floating,"%%s %%%slf %%%slf %%%slf\n",coords_acc,coords_acc,coords_acc);
			for (i=0;i<particles;++i)
			{
				fprintf(fp,floating,species[i],x_1[i],y_1[i],z_1[i]);
			}
			if(tip==1)															// tip
			{
				for (i=0;i<tip_particles;++i)
				{
				fprintf(fp,floating,tip_species[i],tip_x[i],tip_y[i],tip_z[i]);
				}
			}
			fclose(fp);
		}	
		////////////////////////////////////////////////////////////////
		// Velocities and positions export to dat files
		if (velout>0 && (n % velout==0))
		{
			sprintf(buffer,"%s%s%d.dat",velexportpath,vel_file,n+nres);
			fp=fopen(buffer,"w+");
			sprintf(floating,"%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc);
			for (i=0;i<particles;++i)
			{
				fprintf(fp,floating,ux_1[i],uy_1[i],uz_1[i]);
			}
			fclose(fp);
		}
		if (posout>0 && (n % posout==0))
		{
			sprintf(buffer,"%s%s%d.dat",posexportpath,pos_file,n+nres);
			fp=fopen(buffer,"w+");
			sprintf(floating,"%%%slf\t%%%slf\t%%%slf\n",acc,acc,acc);
			for (i=0;i<particles;++i)
			{
				fprintf(fp,floating,x_1[i],y_1[i],z_1[i]);
			}
			fclose(fp);
		}
		////////////////////////////////////////////////////////////////
		if(n%write_res==0)
		{
			fp=fopen(resumepath,"w+");// create the resume file
			// start writting the header
			fprintf(fp,"title\t%s\n",title);// write title
			fprintf(fp,"particles\t%d\n",particles);// write particles
			// write supercell parameters
			fprintf(fp,"supercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",xmin,xmax,ymin,ymax,zmin,zmax);
			// end of header writting
			// write atomic data
			fprintf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tux\tuy\tuz\n");
			sprintf(floating,"%%s\t%%%slf\t%%%slf\t%%%slf\t%%d\t%%d\t%%%slf\t%%%slf\t%%%slf\n",resume_acc,resume_acc,resume_acc,resume_acc,resume_acc,resume_acc);
			for (i=0;i<particles;++i)
			{
				fprintf(fp,floating,species[i],x_1[i],y_1[i],z_1[i],mobileFlag[i],thermoFlag[i],ux_1[i],uy_1[i],uz_1[i]);
			}
			//fprintf(fp,"%d\n",n-1+nres);
			fprintf(fp,"%d\n",n+nres);
			fclose(fp);
			// end of atomic data writing
			
			if(tip==1)
			{
				// tip resume file
				//if(OS==1){sprintf(tip_path,"%s\\out\\%s",sim_folder,tip_resume_file);}
				//else{sprintf(tip_path,"%s/out/%s",sim_folder,tip_resume_file);}
				fp=fopen(tip_resume_path,"w+");
				//fprintf(fp,"%d\n%s\n",tip_particles,tip_title);
				fprintf(fp,"%d\n\n",tip_particles);
				for(i=0;i<tip_particles;++i){fprintf(fp,"%s\t%lf\t%lf\t%lf\n",tip_species[i],tip_x[i],tip_y[i],tip_z[i]);}
				fprintf(fp,"mass_center_velocities:\n%lf\t%lf\t%lf\n",tip_ux,tip_uy,tip_uz);
				fprintf(fp,"holder_coords:\n%lf\t%lf\t%lf\n",xholder,yholder,zholder);
				fclose(fp);
			}
			
		}
		////////////////////////////////////////////////////////////////
	}
	
	// -----------------------------------------------------------------
	// --------------- END OF VERLET -----------------------------------
	// -----------------------------------------------------------------
	
	// Console output: timing info
	// print time to console
	// console
	printf("%s\n",start_time);
	// log
	fplog=fopen(logpath,"a");
	fprintf(fplog,"%s\n",start_time);
	fclose(fplog);										
	time (&rawtime);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	sprintf(stop_time, "Simulation ended at: %s",ctime(&rawtime));
	// stop time
	// console
	printf("%s\nTotal time: %lf\n\n",stop_time,time_spent);
	// log
	fplog=fopen(logpath,"a");
	fprintf(fplog,"%s\nTotal time: %lf\n\n",stop_time,time_spent);
	fclose(fplog);
	
	// Data exports 

	fp=fopen(resumepath,"w+");// create the resume file
	// start writting the header
	fprintf(fp,"title\t%s\n",title);// write title
	fprintf(fp,"particles\t%d\n",particles);// write particles
	// write supercell parameters
	fprintf(fp,"supercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",xmin,xmax,ymin,ymax,zmin,zmax);
	// end of header writting
	// write atomic data
	fprintf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tux\tuy\tuz\n");
	sprintf(floating,"%%s\t%%%slf\t%%%slf\t%%%slf\t%%d\t%%d\t%%%slf\t%%%slf\t%%%slf\n",resume_acc,resume_acc,resume_acc,resume_acc,resume_acc,resume_acc);
	for (i=0;i<particles;++i)
	{
		fprintf(fp,floating,species[i],x_1[i],y_1[i],z_1[i],mobileFlag[i],thermoFlag[i],ux_1[i],uy_1[i],uz_1[i]);
	}
	fprintf(fp,"%d\n",n-1+nres);
	fclose(fp);
	// end of atomic data writing
	
	// export coords
	////////////////////////////////////////////////////////////////////
	sprintf(buffer,"%s%s%d.csv",csvexportpath,csv_file,n-1+nres);
	fp=fopen(buffer,"w+");
	fprintf(fp,"X,Y,Z,radius,Flag,atomicZ\n");
	sprintf(floating,"%%%slf,%%%slf,%%%slf,%%.3lf,%%d,%%d\n",coords_acc,coords_acc,coords_acc);
	if(interaction==0)
	{
		for (i=0;i<particles;++i)
		{
			fprintf(fp,floating,x_1[i],y_1[i],z_1[i],vanderWaalsR[i],mobileFlag[i]+thermoFlag[i],Z[i]);
		}
		if(tip==1)
		{
			for (i=0;i<tip_particles;++i)
			{
				fprintf(fp,floating,tip_x[i],tip_y[i],tip_z[i],tip_vanderWaalsR[i],-1,tip_Z[i]);
			}
		}
	}
	if(interaction==1)
	{
		for (i=0;i<particles;++i)
		{
			fprintf(fp,floating,x_1[i],y_1[i],z_1[i],atomicR[i],mobileFlag[i]+thermoFlag[i],Z[i]);
		}
		if(tip==1)
		{
			for (i=0;i<tip_particles;++i)
			{
				fprintf(fp,floating,tip_x[i],tip_y[i],tip_z[i],tip_atomicR[i],-1,tip_Z[i]);
			}
		}
	}
	fclose(fp);
	////////////////////////////////////////////////////////////////////
	sprintf(buffer,"%s%s%d.xyz",xyzexportpath,xyz_file,n-1+nres);
	fp=fopen(buffer,"w+");
	if(tip==0){fprintf(fp,"%d\n%s\n",particles,title);}					// tip
	if(tip==1){fprintf(fp,"%d\n%s\n",particles+tip_particles,title);}	// tip
	sprintf(floating,"%%s %%%slf %%%slf %%%slf\n",coords_acc,coords_acc,coords_acc);
	for (i=0;i<particles;++i)
	{
		fprintf(fp,floating,species[i],x_1[i],y_1[i],z_1[i]);
	}
	if(tip==1)															// tip
	{
		for (i=0;i<tip_particles;++i)
		{
		fprintf(fp,floating,tip_species[i],tip_x[i],tip_y[i],tip_z[i]);
		}
	}
	fclose(fp);
	////////////////////////////////////////////////////////////////////
	/*
	// VACF export
	if (VACFcalcstep>0  && VACFcalcstep<=myN)
	{
		fp=fopen(VACFexportpath,"w+");
		sprintf(floating,"%%d\t%%%slf\n",acc);
		for (i=0;i<VACFcalcstep;++i){fprintf(fp,floating,i,VACF[i]);}
		fclose(fp);
	}
	*/
	
//(E) --- Export tip positions for resuming simulations ----------------
	if(tip==1)
	{
		// tip resume file
		//if(OS==1){sprintf(tip_path,"%s\\out\\%s",sim_folder,tip_resume_file);}
		//else{sprintf(tip_path,"%s/out/%s",sim_folder,tip_resume_file);}
		fp=fopen(tip_resume_path,"w+");
		//fprintf(fp,"%d\n%s\n",tip_particles,tip_title);
		fprintf(fp,"%d\n\n",tip_particles);
		for(i=0;i<tip_particles;++i){fprintf(fp,"%s\t%lf\t%lf\t%lf\n",tip_species[i],tip_x[i],tip_y[i],tip_z[i]);}
		fprintf(fp,"mass_center_velocities:\n%lf\t%lf\t%lf\n",tip_ux,tip_uy,tip_uz);
		fprintf(fp,"holder_coords:\n%lf\t%lf\t%lf\n",xholder,yholder,zholder);
		fclose(fp);
	}
//(E END) --------------------------------------------------------------
	
	// Free allocated memory
	free(x_1);free(y_1);free(z_1);
    free(ux_1);free(uy_1);free(uz_1);
	free(ax_1);free(ay_1);free(az_1);
    free(ax_2);free(ay_2);free(az_2);
	free(thermostat);
    free(mass);free(vanderWaalsR);free(atomicR);free(Z);free(mobileFlag);free(thermoFlag);
	for (i=0;i<particles;++i)
	{
		free(species[i]);
	}
	free(species);
	free(sumFx);free(sumFy);free(sumFz);
	if (VACFcalcstep>0 && VACFcalcstep<=myN)
	{
	//	free(VACF);
		free(uxsave);free(uysave);free(uzsave);
	}
	if (linked_cells==1)
	{
		//for (i=0;i<26;++i){free(neighbors[i]);}
		free(p);
		free(neighbors);
		free(cellpartarray);free(neighborcellpartarray);free(head);free(list);
		free(allpartarray);
	}
	
//(F) --- Free tip-related arrays --------------------------------------
	if(tip==1)
	{
		free(tip_x);free(tip_y);free(tip_z);
		free(tip_mass);
		for (i=0;i<tip_particles;++i)
		{
			free(tip_species[i]);
		}
		free(tip_species);
		free(tip_Z);
		free(tip_sumFx);free(tip_sumFy);free(tip_sumFz);
		free(tip_atomicR);free(tip_vanderWaalsR);
	}
	
	if(tip==1 && tip_linked_cells==1)
	{
		free(tip_cellpartarray);free(tip_neighborcellpartarray);free(tip_head);free(tip_list);
		free(tip_neighbors);free(p1);
	}
//(F END) --------------------------------------------------------------
	
	if(interaction==1){free(dx_zeta_III);free(dy_zeta_III);free(dz_zeta_III);}		//<-------------------
	
	//
	free(p_book);free(book);
	
	//
	free(cells);
	free(sumFxrecv);
	free(sumFyrecv);
	free(sumFzrecv);
	
	//
	free(barostatX);
	free(barostatY);
	free(barostatZ);
	free(bflagX);
	free(bflagY);
	free(bflagZ);
	
	free(full);
	
	}
	else
	{
		
		int i;															// index
		int n,myN;														// time index and integration steps
		int particles;													// number of particles
		int M3;															// number of cells	
		int *Z;															// atomic number array
		int *p,**neighbors;												// neighbors matrix arrays
		double Lx,Ly,Lz;												// box lengths
		int *head,*list;												// head and list arrays
		double *x_1,*y_1,*z_1;											// positions
		double *sumFx,*sumFy,*sumFz;									// total force arrays
		double V;														// potential energy
		int *cellpartarray,*neighborcellpartarray,*allpartarray;		// cell related arrays
		double *dx_zeta_III,*dy_zeta_III,*dz_zeta_III;					// T3 zeta arrays
		int interaction;												// interaction type
		double sigmaxx,sigmayy,sigmazz;									// internal stress components
		// mpi related
		int dest,source,id;
		int cellsN;
		int *cells;
		// barostat related
		int change_box;
		double mixx,miyy,mizz;
		int realloc_flag;
		
		// receive
		source=0;
		id=0;
		MPI_Recv(&particles, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		
		Z=(int*)malloc(particles*sizeof(int));
		x_1=(double*)malloc(particles*sizeof(double));
		y_1=(double*)malloc(particles*sizeof(double));
		z_1=(double*)malloc(particles*sizeof(double));
		// total forces
		sumFx=(double*)malloc(particles*sizeof(double));
		sumFy=(double*)malloc(particles*sizeof(double));
		sumFz=(double*)malloc(particles*sizeof(double));
		
		id=1;
		MPI_Recv(&M3, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		
		cellpartarray=(int*)malloc(particles*sizeof(int));// cells
		neighborcellpartarray=(int*)malloc(particles*sizeof(int));
		// ------*
		allpartarray=(int*)malloc(particles*sizeof(int));
		dx_zeta_III=(double*)malloc(particles*sizeof(double));			//<-----------------------
		dy_zeta_III=(double*)malloc(particles*sizeof(double));
		dz_zeta_III=(double*)malloc(particles*sizeof(double));
		// ------*
		head=(int*)malloc(M3*sizeof(int));// head
		list=(int*)malloc(particles*sizeof(int));// list
		neighbors=(int**)malloc(26*sizeof(int*));
		p=(int*)malloc(26*M3*sizeof(int));
		for(i=0;i<26;++i)
		{
			neighbors[i]=p+i*M3;
		}
		
		//cellsN=(int)floor(M3/size)+1;
		//cells=(int*)malloc(cellsN*sizeof(int));	
		
		id=2;
		MPI_Recv(&Z[0], particles, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		
		id=3;
		MPI_Recv(&neighbors[0][0], 26*M3, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		
		id=4;
		MPI_Recv(&Lx, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
		id=5;
		MPI_Recv(&Ly, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
		id=6;
		MPI_Recv(&Lz, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
		
		//id=13;
		//MPI_Recv(&cells[0], cellsN, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		
		//
		id=12;
		MPI_Recv(&myN, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		id=19;
		MPI_Recv(&interaction, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		
		//------------------
		
		//--------------***********
		
		id=7;
		MPI_Recv(&x_1[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
		id=8;
		MPI_Recv(&y_1[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
		id=9;
		MPI_Recv(&z_1[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
		id=10;
		MPI_Recv(&head[0], M3, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		id=11;
		MPI_Recv(&list[0], particles, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		
		id=90;
		MPI_Recv(&cellsN, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		cells=(int*)malloc(cellsN*sizeof(int));
		
		id=13;
		MPI_Recv(&cells[0], cellsN, MPI_INT, source, id, MPI_COMM_WORLD, &status);
		
		if(interaction==1)
		{
			V_F_virial_cell_calc_tersoff_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,allpartarray,dx_zeta_III,dy_zeta_III,dz_zeta_III,cellsN,cells);
		}
		
		if(interaction==0)
		{
			V_F_virial_cell_calc_LJ_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,LJcutoff,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,cellsN,cells);
		}
		
		dest=0;
		id=14;
		MPI_Send(&V, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
		
		dest=0;
		id=20;
		MPI_Send(&sigmaxx, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
		id=21;
		MPI_Send(&sigmayy, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
		id=22;
		MPI_Send(&sigmazz, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				
		dest=0;
		id=16;
		MPI_Send(&sumFx[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
		id=17;
		MPI_Send(&sumFy[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
		id=18;
		MPI_Send(&sumFz[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			
		for (n=1;n<=myN;++n)												// time loop
		{
			
			// receive bookkeeping data
			
			id=7;
			MPI_Recv(&x_1[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			id=8;
			MPI_Recv(&y_1[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			id=9;
			MPI_Recv(&z_1[0], particles, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
			id=10;
			MPI_Recv(&head[0], M3, MPI_INT, source, id, MPI_COMM_WORLD, &status);
			id=11;
			MPI_Recv(&list[0], particles, MPI_INT, source, id, MPI_COMM_WORLD, &status);
			
			id=90;
			MPI_Recv(&cellsN, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
			//cells=(int*)malloc(cellsN*sizeof(int));
			cells=realloc(cells,cellsN*sizeof(int));
			
			id=13;
			MPI_Recv(&cells[0], cellsN, MPI_INT, source, id, MPI_COMM_WORLD, &status);
			
			
			
			if(interaction==1)
			{
				V_F_virial_cell_calc_tersoff_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,allpartarray,dx_zeta_III,dy_zeta_III,dz_zeta_III,cellsN,cells);
			}
			
			if(interaction==0)
			{
				V_F_virial_cell_calc_LJ_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,LJcutoff,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,cellsN,cells);
			}
			
			dest=0;
			id=14;
			MPI_Send(&V, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			
			dest=0;
			id=20;
			MPI_Send(&sigmaxx, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			id=21;
			MPI_Send(&sigmayy, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			id=22;
			MPI_Send(&sigmazz, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			
			dest=0;
			id=16;
			MPI_Send(&sumFx[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			id=17;
			MPI_Send(&sumFy[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			id=18;
			MPI_Send(&sumFz[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			
			// --- A section ---
			
			// receive change_box
			id=40;
			MPI_Recv(&change_box, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
			
			if(change_box==1)
			{
				
				// receive mi
				id=41;
				MPI_Recv(&mixx, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=42;
				MPI_Recv(&miyy, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=43;
				MPI_Recv(&mizz, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				
				id=4;
				MPI_Recv(&Lx, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=5;
				MPI_Recv(&Ly, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				id=6;
				MPI_Recv(&Lz, 1, MPI_DOUBLE, source, id, MPI_COMM_WORLD, &status);
				
				// rescale
				for(i=0;i<particles;++i)
				{
					x_1[i]=mixx*x_1[i];
					y_1[i]=miyy*y_1[i];
					z_1[i]=mizz*z_1[i];
				}
				
				id=44;
				MPI_Recv(&realloc_flag, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
								
				if(realloc_flag==1)
				{
					id=45;
					MPI_Recv(&M3, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					head=realloc(head,M3*sizeof(int));
					p=realloc(p,26*M3*sizeof(int));
					for(i=0;i<26;++i)
					{
						neighbors[i]=p+i*M3;
					}
					//cellsN=(int)floor(M3/size)+1;
					//cells=realloc(cells,cellsN*sizeof(int));
					id=46;
					MPI_Recv(&neighbors[0][0], 26*M3, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					
					id=90;
					MPI_Recv(&cellsN, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					cells=realloc(cells,cellsN*sizeof(int));
					
					id=47;
					MPI_Recv(&cells[0], cellsN, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					id=48;
					MPI_Recv(&head[0], M3, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					id=49;
					MPI_Recv(&list[0], particles, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					
				}
				else
				{		
					id=60;
					MPI_Recv(&head[0], M3, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					id=61;
					MPI_Recv(&list[0], particles, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					
					id=90;
					MPI_Recv(&cellsN, 1, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					cells=realloc(cells,cellsN*sizeof(int));
					
					id=47;
					MPI_Recv(&cells[0], cellsN, MPI_INT, source, id, MPI_COMM_WORLD, &status);
					
				}
								
				if(interaction==1)
				{
					V_F_virial_cell_calc_tersoff_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,allpartarray,dx_zeta_III,dy_zeta_III,dz_zeta_III,cellsN,cells);
					
				}
				
				if(interaction==0)
				{
					V_F_virial_cell_calc_LJ_mpi(x_1,y_1,z_1,particles,Z,Lx,Ly,Lz,LJcutoff,sumFx,sumFy,sumFz,&V,&sigmaxx,&sigmayy,&sigmazz,M3,head,list,neighbors,cellpartarray,neighborcellpartarray,cellsN,cells);
				}
								
				id=50;
				MPI_Send(&V, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				
				id=20;
				MPI_Send(&sigmaxx, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				id=21;
				MPI_Send(&sigmayy, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				id=22;
				MPI_Send(&sigmazz, 1, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
			
				id=16;
				MPI_Send(&sumFx[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				id=17;
				MPI_Send(&sumFy[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				id=18;
				MPI_Send(&sumFz[0], particles, MPI_DOUBLE, dest, id, MPI_COMM_WORLD);
				
			}
			// --- end of A section ---
			
		}
		
		//--------------***********
		
		free(cellpartarray);free(neighborcellpartarray);
		free(sumFx);free(sumFy);free(sumFz);
		free(x_1);free(y_1);free(z_1);
		free(head);free(list);
		free(Z);
		free(p);
		free(neighbors);
		//----------------*
		free(allpartarray);
		free(dx_zeta_III);
		free(dy_zeta_III);
		free(dz_zeta_III);
		//----------------*
		//
		free(cells);
		
	}
	
	MPI_Finalize();
		
	return 0;
}
