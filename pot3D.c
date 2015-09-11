/* HEAT3D: Versao C paralelizada
 * ARQUIVO: pot3D.c
 *
 * DESCRICAO:
 *	  
 * AUTOR: Wellington Marinho
 *        Aluno de Iniciacao Cientifica - Faculdade de Engenharia Mecanica (UFU)
 *        Laboratorio de Transferencia de Calor e Massa e Dinamica de Fluidos (LTCM)
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 0.00001
#define CMD printf("%s\n",comand);

/******* Definicoes para reduzir o codigo ********/
#define Tn(i,j,k) grid[(i)*NY*NZ + (j)*NZ + (k)] 
#define Tp(i,j,k) tmp[(i)*NY*NZ + (j)*NZ + (k)]
#define B(i,j,k) b[(i)*NY*NZ + (j)*NZ + (k)]
#define for_(i,n) for(i=0; i<n; i++)

//#define HEAT3D
//#define _STAT 		
#define _SIM
//#define PARALLEL

#ifdef PARALLEL
#include "mpi.h"
#endif

/**************************** INICIALIZACAO DAS FUNCOES ***************************/

FILE *open_file(char name[20], char op[3]);
void print_mat(const double *mat, unsigned m, unsigned n, unsigned p);
void fprint_mat(FILE *fp, const double *mat, unsigned m, unsigned n, unsigned p);
void fprint_header_vtk(FILE *fvtk, const char *title, const char *format, const char *geometry);
void fprint_points(FILE *fpoint, double *px, double *py, double *pz, unsigned nx, unsigned ny, unsigned nz);
void fprint_vector(FILE *fpoint, double *p, unsigned n);

/**********************************************************************************/

int elipse(double x, double y, double z, double xc, double yc, double zc, double a, double b, double c, double r);

int cilinder(double x, double y, double z, double xc, double yc, double zi, double zf, double a, double b, double r);

int bar(double x, double y, double z, double xi, double yi, double zi, double xf, double yf, double zf);	
	
/**********************************************************************************/

int main (int argc, char *argv[])
{

#ifdef PARALLEL
	int id, procs;
	unsigned long int size_msg;
	MPI_Status status;
#endif

unsigned xi, yi, zi, xf, yf, zf, xc, yc, zc;
int t0=0, t1=0, partial_time=0, total_time=0;
unsigned i, j, k, iter, nx, x1, x2, NX, NY, NZ, nz;
unsigned long int ntot;
float w, wi, wf, dw;

double Ap, An, As, Ae, Aw, At, Ab;
double *grid, *tmp, *b, *px, *py, *pz;
double x,y,z;
double alpha, rho, c, K, t, tf, dt;
double dx, dy, dz;


double  maxdelta, maximum, delta;
FILE *init;
	
#ifdef _STAT
	FILE *arq, *f_w;
	char op, ch;
	int  best_iter=32000, best_time=0;
	double best_w=0.0, b_delta=0.0;
#endif

#ifdef _SIM
	FILE *data, *r_w;
	char op, str[3], c0 = '0', c1='0';
	int ctrl = 0;
#endif

        char file[20], filename[20], dir[30], dir_w[30], Dir[30], _dir[30], comand[20], comand1[10];
        char format[6];
        char f_omega[6] = "omega";



/**************************** INICIALIZACOES DO MPI ******************************/
#ifdef PARALLEL
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

#endif /*PARALLEL*/	 

/*************************** INICIALIZACOES *************************************/

init = open_file("init_data","r");

fscanf(init,"%*s%f %*s%f %*s%f %*s%f %*s%lf %*s%lf %*s%lf %*s%lf %*s%lf %*s%lf",&w,&wi,&wf,&dw,&rho,&c,&K,&t,&tf,&dt);
fscanf(init,"%*s%lf %*s%d %*s%lf %*s%d %*s%lf %*s%d",&x,&NX,&y,&NY,&z,&NZ);

if(NZ == 1) nz=3;
else nz = NZ;
ntot = NX*NY*NZ;

dx = (double)x/(NX-1);
dy = (double)y/(NY-1);
dz = (double)z/(nz-1);

Aw = Ae = 1.00/(dx*dx);
An = As = 1.00/(dy*dy);
At = Ab = 1.00/(dz*dz);
Ap = An + As + Ae + Aw + At + Ab;

//alpha = rho*c/K*dt;

//Ap = Ap + alpha;

#ifdef PARALLEL
	size_msg = NY*NZ;
 	nx = NX/procs;
	if((NX%2)!=0){
                printf("\nERRO: \nA dimensao em x deve ser divisivel pelo numero de processadores\n");
                exit(1);
        }

	x1 = id*nx;
	x2 = (id+1)*nx - 1;
	if (id == 0) x1++;
        if (id == procs-1) x2--;
#else
	x1 = 1;
	x2 = NX-2;
#endif /*PARALLEL*/

#ifdef _STAT

#ifdef PARALLEL
   if(id == 0){
#endif /*PARALLEL*/ 	   
	
	strcpy(format,".dat");

#ifdef HEAT3D
#else
	system("clear");
	printf("Definicao: _STAT ");
#endif
	printf("\nAvalia performance do programa variando o parametro omega (w)\n ");
	printf("\nNome do Diretorio: ");
        scanf("%s",&dir);
	
	strcpy(comand,"cd ");
	strcat(comand,dir);
#ifdef HEAT3D
	 strcpy(comand,"mkdir ");
         strcat(comand,dir);
//         CMD
//	 system(comand);	
	 strcat(comand, "/statistic ");
 //        strcat(comand,dir);
	 CMD
	 system(comand);	 
#else
	if(system(comand)== 0)
		printf("\nDiretorio ja existe\n");
	else{
	strcpy(comand,"mkdir ");
	strcat(comand,dir);
	CMD
	system(comand); 
	strcat(comand, "/statistic ");
        strcat(comand,dir);
	strcat(comand, "/simulation");	
	CMD
	system(comand); 
	}
#endif
	printf("Nome do Arquivo (.dat):  ");
	scanf("%s",&filename);
	
	strcat(filename,format);
	strcat(dir,"/statistic/");

	strcpy(dir_w,dir);
	
	strcat(dir,filename);
	strcat(dir_w,"omega");

	arq = open_file(dir,"w");
	f_w = open_file(dir_w,"w+");
	
	printf("\nRESULTADOS\ndim %d %d %d\n",NX,NY,NZ);
#ifdef PARALLEL
	printf("Numero de processadores: %d\n", procs);
	printf("Cada processador calculou %d nos\n", (x2-x1+2)*NY*NZ);
#endif 
	
	
	printf("\nomega    iteracoes    epsilon    tempo \n");

	fprintf(arq,"\nRESULTADOS\ndim %d %d %d\n",NX,NY,NZ);
#ifdef PARALLEL
        fprintf(arq,"Numero de processadores: %d\n", nprocs);
        fprintf(arq,"Cada processador calculou %d nos\n", (x2-x1+2)*NY*NZ);              
#endif 
	fprintf(arq,"\nomega    iteracoes    epsilon    tempo \n");

#ifdef PARALLEL	
}
#endif
#ifdef PARALLEL
        printf("Numero de processadores: %d\n", procs);
	printf("Processador <%d> processou %d nos", id, (x2-x1+2)*NY*NZ);
#endif
		
/************************ ALOCACAO DINAMICA DE MEMORIA ***************************/

	grid = malloc(sizeof(double)*ntot);
	tmp = malloc(sizeof(double)*ntot);
	b = malloc(sizeof(double)*ntot);
#ifdef PARALLEL
   if(id == 0){
#endif	   
	px = malloc(sizeof(double)*NX);
	py = malloc(sizeof(double)*NY);
	pz = malloc(sizeof(double)*NZ);
#ifdef PARALLEL
   }
#endif

for(w = wi; w < wf; w += dw) {
   iter = 0;
  	
#endif /*_STAT*/

#ifdef _SIM

#ifdef PARALLEL
   if(id == 0){
#endif
	strcpy(format,".vtk");
	system("clear");
	printf("Definicao: _SIM \n");
	printf("Resolve a equacao de Laplace gerando arquivo formato vtk\n ");

        printf("\nNome do Diretorio: ");
        scanf("%s",&dir);

	strcpy(dir_w,dir);
	strcpy(Dir,dir);
	strcat(dir_w,"/statistic/");
	strcpy(comand,"cd ");

	strcat(comand,dir);
	strcat(comand, "/simulation");
	if(system(comand)==0);
	else{
		strcpy(comand,"mkdir ");
        	strcat(comand,dir);
		CMD
		system(comand);
        	strcat(comand, "/simulation");
		CMD
        	system(comand);
	}			
	printf("Nome do Arquivo (.vtk):  ");
	scanf("%s",&file);

	strcpy(filename,file);
	str[0] = c0;
	str[1] = c1;
	str[2] = '\0';
        strcat(filename,str);	
       	strcat(filename,format);
        strcat(dir,"/simulation/");

        strcat(dir,filename);
	c1++;

/**********************************************************************************/
#ifdef HEAT3D
	        system("./STAT_heat3D");
#endif
/**********************************************************************************/
		
	
	data = open_file(dir,"w");
	
        strcpy(comand,"cd ");
        strcat(comand,dir_w);
	    
   	if(system(comand)==0)
        {
         strcat(dir_w,"omega");
	 
	 r_w = open_file(dir_w,"r");
	 
           fscanf(r_w,"%f",&w);
           ctrl=1;
          }
	
	
	printf("\nProcessando ...\n");
	printf("\nRESULTADOS\n");
	printf("dim %d %d %d\nomega    iteracoes    epsilon    tempo \n",NX,NY,NZ);
	
#ifdef PARALLEL
   }
#endif   
/************************ ALOCACAO DINAMICA DE MEMORIA ***************************/

	grid = malloc(sizeof(double)*ntot);
	tmp = malloc(sizeof(double)*ntot);
	b = malloc(sizeof(double)*ntot);
#ifdef PARALLEL
   if(id == 0){
#endif
	px = malloc(sizeof(double)*NX);
	py = malloc(sizeof(double)*NY);
	pz = malloc(sizeof(double)*NZ);
#ifdef PARALLEL
   }
#endif   
	
/******************************* INICIALIZACOES **********************************/
#ifdef PARALLEL
      if(id == 0){
#endif
//Gera grid
	  for_(i,NX) 
	     px[i] = i*dx;
	  for_(i,NY) 
	      py[i] = i*dy; 
          for_(i,NZ) 
              pz[i] = i*dz; 
	

#ifdef PARALLEL
   }
#endif      
	iter = 0;
#endif  /*_SIM*/

    for_(i,NX)
       for_(j,NY)	    
	  for_(k,NZ)
             Tn(i,j,k) = Tp(i,j,k)=0.0;	  

/************************** Condicoes de contorno ********************************/
   xc = NX/2;	//COORDENADAS DE CENTRO
   yc = NY/2;	//COORDENADAS DE CENTRO
   zc = NZ/2;   //COORDENADAS DE CENTRO
   xi =  xc - 0.1*NX;
   xf =  xc + 0.1*NX;
   yi =  yc - 0.1*NY;
   yf =  yc + 0.1*NY;
   zi =  zc - 0.1*NZ;
   zf =  zc + 0.1*NZ;
	
   
   for_(i,NX)
       for_(j,NY)
          for_(k,NZ){
//	     if(i==0)
//            	Tn(i,j,k)=Tp(i,j,k)=10;
//	    if(i==NX-1)	
// 	      	Tn(NX-1,j,k) = Tp(NX-1,j,k) = 10;
  	
/*	     if(j==0)
	       Tn(i,j,k)=Tp(i,j,k)=100;
	     if(j==NY-1)
	       Tn(i,NY-1,k) = Tp(i,j,NY-1) = 100;
*/
	    /* if(k==0)
               Tn(i,j,k)=Tp(i,j,k)=100;
             if(k==NZ-1)
               Tn(i,j,NZ-1) = Tp(i,j,NZ-1) = 100;
	*/		  
//             if(bar(i,j,k,xi,0,0,xi,NY-1,NZ-1))		  
//		     Tn(i,j,k) = Tp(i,j,k)=-10.0;
		
//	     if(bar(i,j,k,xf,0,0,xf,NY-1,NZ-1))
//	    if(bar(i,j,k,xi,yi,zi,xf,yf,zf))
//	     if(cilinder(i,j,k,xc,yc,zi,zf,1,1,0.1*NZ))
	     if(elipse(i,j,k,xc,yc+20,zc,2,2,1,0.1*NX))
//             Tn(i,j,k) = Tp(i,j,k)=25.0;
//	     if(elipse(i,j,k,xc,yc,zc,1,1,1,0.1*NX))
	     Tn(i,j,k) = Tp(i,j,k)=10.0;

	     
	  } 
//print_mat(grid,NX,NY,NZ);
    
//  t0=time(0);
					 
/****************************** CALCULOS *****************************************/
#ifdef _SIM
for(; t<=tf; t+=dt){    
#endif

 for_(i,NX)
    for_(j,NY)
       for_(k,NZ)
           B(i,j,k)=0.0;//alpha*Tn(i,j,k);
 
    for (;;) {
	maxdelta=0.0;	       
        for (i=x1; i<=x2; i++)
            for (j=1; j <= NY-2; j++)
               for (k=1; k<=nz-2; k++) {
//		 if(bar(i,j,k,xc,yi,0,xc,yf,NZ-1))
//		 if(bar(i,j,k,xi,0,0,xi,NY-1,NZ-1))
//		 if(cilinder(i,j,k,xc,yc,zi,zf,1,1,4))
   		 if(elipse(i,j,k,xc,yc+20,zc,1,1,1,0.1*NX))
//			 continue;
//		  if(elipse(i,j,k,xc,yc,zc,1,1,1,0.1*NX))
//		    if(bar(i,j,k,xi,yi,zi,xf,yf,zf))
//		 if(bar(i,j,k,xf,0,0,xf,NY-1,NZ-1))
			  continue;
		       
	         Tp(i,j,k) = ( Aw*Tn(i-1,j,k) + Ae*Tn(i+1,j,k) +
  			      As*Tn(i,j-1,k) + An*Tn(i,j+1,k) + 
			      Ab*Tn(i,j,k-1) + At*Tn(i,j,k+1) + B(i,j,k))*w/Ap +
	            	      (1.0 - w)*Tp(i,j,k);
		  delta = Tp(i,j,k) - Tn(i,j,k);
		  if(delta > maxdelta)
			   maxdelta = fabs(delta);
		  Tn(i,j,k) = Tp(i,j,k);
                 } 
    
/**************************** PASSAGEM DE MENSAGENS ******************************/ 

#ifdef PARALLEL
       if (id != procs-1){
          MPI_Send(&grid[x2], size_msg, MPI_DOUBLE, id+1, 100, MPI_COMM_WORLD);
          MPI_Recv(&grid[x2+1], size_msg, MPI_DOUBLE, id+1, 200, MPI_COMM_WORLD,&status);
		}	

        if (id != 0) {
          MPI_Recv(&grid[x1-1], size_msg, MPI_DOUBLE, id-1, 100, MPI_COMM_WORLD, &status);
	  MPI_Send(&grid[x1], size_msg, MPI_DOUBLE, id-1, 200, MPI_COMM_WORLD);
	   }
#endif    
	
iter++;
#ifdef PARALLEL	     
	MPI_Allreduce(&maxdelta, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
	maximum = maxdelta; 
#endif
	if (maximum < EPS) break;
	}
    
	
#ifdef PARALLEL
  if (id != 0)
        MPI_Send(&grid[x1], ((x2-x1)+2)*NY*NZ, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);     

  else {for(i=1; i< procs; i++)
        MPI_Recv(&grid[i*(x2+1)], ((x2-x1)+2)*NY*NZ, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status);
   }
#endif

#ifdef PARALLEL
   if(id == 0){
#endif


   partial_time = t1-t0;
   printf("%5.3f%13d%11.6f%8.3lf\n",w, iter, maximum, t);
   
#ifdef _STAT       	
        fprintf(arq,"%5.3f%13d%11.6f%8ds\n",w, iter, maximum, partial_time);

	if(iter < best_iter) {
	  best_iter =iter;
	  best_w = w;
	  b_delta = maxdelta;
	  best_time = partial_time;
	}
   }
 	printf("\n\nMELHOR RESULTADO\n");
#ifdef PARALLEL
	printf("Numero de processadores %d\n", nprocs);
#endif
        printf("dim %d %d %d\nomega    iteracoes    epsilon    tempo	 ", NX, NY, NZ);
	printf("\n%5.3f%13d%11.6f%8ds\n",best_w, best_iter, b_delta, best_time);
	fprintf(f_w,"%f",best_w);		 
	fclose(arq);
	fclose(f_w);
#endif /*_STAT*/      	       

/********************* GERA ARQUIVO EM FORMATO VTK **********************/
#ifdef _SIM

	printf("Gerando arquivo vtk ...\n");
	fprint_header_vtk(data, "teste vtk", "ASCII", "STRUCTURED_POINTS");
	fprintf(data,"DIMENSIONS %d %d %d\n", NX, NY, NZ);
//	fprintf(data, "POINTS %d float\n", NX*NY*NZ);
//	fprint_points(data, px, py, pz, NX, NY, NZ);
        fprintf(data,"SPACING %f %f %f\n", dx, dy, dz);
        fprintf(data, "ORIGIN 0 0 0\n");
	/*
	fprintf(data, "\nX_COORDINATES %d float\n", NX);
	fprint_vector(data, px, NX);
	fprintf(data, "\nY_COORDINATES %d float\n", NY);
	fprint_vector(data, py, NY);
	fprintf(data, "\nZ_COORDINATES %d float\n", NZ);
	fprint_vector(data, pz, NZ);
	*/
	fprintf(data, "\n\nPOINT_DATA %d\n", ntot);
	fprintf(data, "SCALARS temperatura float 1\n");
	fprintf(data, "LOOKUP_TABLE default\n");
//	fprintf(data, "\nFIELD escalar 1\n");
//      fprintf(data, "temperatura 1 %d 1\n",ntot);
	fprint_mat(data,grid,NX,NY,NZ);
       
	fclose(data);
	if(ctrl)
	fclose(r_w);
       // printf("\nt=%f\n",t);
	if(t<tf)
	{  
		strcpy(filename,file);
        	str[0] = c0;
		str[1] = c1;
		strcat(filename,str);
		strcat(filename,format);
		strcpy(_dir,Dir);
		strcat(_dir,"/simulation/");
  		strcat(_dir,filename);
		if(c1=='9'){
			c1 = '0';
			c0++;
		}
   	        c1++;
	        data = open_file(_dir,"w");
		iter=0;
	}
     
   }
#endif
 fclose(init);
#ifdef PARALLEL
	}
	MPI_Finalize();
#endif

return 0; 
}
		    
 /****************************** FUNCOES  ****************************************/

FILE *open_file(char name[30], char op[3] )
{
	FILE *fp;
	if(!(fp =fopen(name,op)))
	{
		printf("\nERRO: Arquivo \"%s\" nao pode ser aberto\n",name);
		exit(1);
	}
	return fp;	
}
void print_mat(const double *mat, unsigned m, unsigned n, unsigned p)
{
	unsigned i, j, k;
	for(i=0; i<m; i++) {
	   for(j=0; j<n; j++){
	      for(k=0; k<p; k++)
	        printf("%8.4f ",mat[i*n*p + j*p + k]);
	   printf("\n");
	}printf("\n");}
}

void fprint_mat(FILE *fp, const double *mat, unsigned m, unsigned n, unsigned p)
{
	unsigned i, j, k;
	for(k=0; k<p; k++){ 
	   for(j=0; j<n; j++){
	      for(i=0; i<m; i++)	   
	         fprintf(fp,"%8.4f ",mat[i*n*p + j*p+ k]);
           fprintf(fp,"\n");
	   }fprintf(fp,"\n");}   
}
void fprint_header_vtk(FILE *fvtk, const char *title, const char *format, const char *geometry)
{

	fprintf(fvtk,"# vtk DataFile Version 2.0\n");
	fprintf(fvtk,"%s\n",title);
	fprintf(fvtk,"%s\n",format);
	fprintf(fvtk,"DATASET %s\n",geometry);
}
void fprint_points(FILE *fpoint, double *px, double *py, double *pz, unsigned nx, unsigned ny, unsigned nz)
{
	unsigned i, j, k;
	for(i=0; i<nx; i++)
		for(j=0; j<ny; j++)
			for(k=0; k<nz; k++)
				fprintf(fpoint,"%f %f %f \n", px[k], py[j], pz[i]);
			
}
void fprint_vector(FILE *fpoint, double *p, unsigned n)
{
        unsigned i;
        for(i=0; i<n; i++)
          fprintf(fpoint,"%f ", p[i]);

}
/***********************************************************************************/
int elipse(double x, double y, double z, double xc, double yc, double zc, double a, double b, double c, double r)
{
if(((x-xc)*(x-xc)/a + (y-yc)*(y-yc)/b + (z-zc)*(z-zc)/c) <= r*r)
	return 1;
else
	return 0;
}
int cilinder(double x, double y, double z, double xc, double yc, double zi, double zf, double a, double b, double r)
{
if(( ((x-xc)*(x-xc)/a + (y-yc)*(y-yc)/b) <= r*r) && (z>=zi) && (z<=zf))
	        return 1;
else
	        return 0;

}
int bar(double x, double y, double z, double xi, double yi, double zi, double xf, double yf, double zf)
{
if((x>=xi) && (x<=xf) && (y>=yi) && (y<=yf) && (z>=zi) && (z<=zf)) 
	return 1;
else
	return 0;
		       
}
/***********************************************************************************/



