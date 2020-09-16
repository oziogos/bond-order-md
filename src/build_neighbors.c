//
#include<stdlib.h>

void build_neighbors(int Mx, int My, int Mz, int **neighbors);
void build_neighbors(int Mx, int My, int Mz, int **neighbors)
{
	int i,j,k,l;
	int ghostsx,ghostsy,**cells,row,column,ii,jj;
	int M2,M3;
	
	// cells
	M2=Mx*My;
	M3=M2*Mz;
	
	// ghost cells
	ghostsx=Mx+2;
	ghostsy=My+2;
	cells=(int**)malloc(ghostsy*sizeof(int*));
	for (i=0;i<ghostsy;++i)
	{
		cells[i]=(int*)malloc(ghostsx*sizeof(int));
	}
	
	for (i=0;i<ghostsy;++i){for (j=0;j<ghostsx;++j){cells[i][j]=0;}}

	// index rule: subtract 1 from static values and keep running indices the same! reminder: place counters at the bottom!!
	l=0;
	for (i=My-1;i>=0;--i)
	{
		for (j=0;j<=Mx-1;++j)
		{
			l=l+1;
			cells[i+1][j+1]=l;
			
		}
	}

	for (j=1;j<=ghostsx-2;++j)
	{
		cells[0][j]=cells[ghostsy-2][j];
		cells[ghostsy-1][j]=cells[1][j];
	}
	
	for (i=0;i<=ghostsy-1;++i)
	{
		cells[i][ghostsx-1]=cells[i][1];
		cells[i][0]=cells[i][ghostsx-2];
	}

	for (i=0;i<26;++i){for (j=0;j<M3;++j){neighbors[i][j]=0;}}

	l=0;
	for (i=My-1;i>=0;--i)
	{
		for (j=0;j<=Mx-1;++j)
		{
			row=i+2;column=j+2;
			k=0;
			for (jj=0;jj<=2;++jj)
			{
				for (ii=0;ii<=2;++ii)
				{
					if (ii!=1||jj!=1)
					{
						neighbors[k][l]=cells[ii-2+row][jj-2+column];
						k=k+1;
					}
				}
			}
			l=l+1;
		}
	}

	for (j=M2;j<=M3-1;++j)
	{
		for (i=0;i<=7;++i)
		{
			neighbors[i][j]=neighbors[i][j-M2]+M2;
		}
	}

	for (j=0;j<=M3-1;++j)
	{
		for (i=0;i<=7;++i)
		{
			if (neighbors[i][j]+M2<=M3)
			{
				neighbors[i+8][j]=neighbors[i][j]+M2;
			}
			else
			{
				neighbors[i+8][j]=neighbors[i][j]-(Mz-1)*M2;
			}
			if (neighbors[i][j]-M2>0)
			{
				neighbors[i+16][j]=neighbors[i][j]-M2;
			}
			else
			{
				neighbors[i+16][j]=neighbors[i][j]+(Mz-1)*M2;
			}
		}
		if (j+1+M2<=M3)
		{
			neighbors[24][j]=j+1+M2;
		}
		else
		{
			neighbors[24][j]=j+1-(Mz-1)*M2;
		}
		if (j+1-M2>0)
		{
			neighbors[25][j]=j+1-M2;
		}
		else
		{
			neighbors[25][j]=j+1+(Mz-1)*M2;
		}
	}

	for (i=0;i<26;++i)
	{
		for (j=0;j<M3;++j)
		{
			neighbors[i][j]=neighbors[i][j]-1;
		}

	}

	for (i=0;i<ghostsy;++i){free(cells[i]);}
	free(cells);
}
