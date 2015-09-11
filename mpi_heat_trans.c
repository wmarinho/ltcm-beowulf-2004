/* HEAT3D: Versao C paralelizada
 * ARQUIVO: mpi_heat_trans.c
 *
 * DESCRICAO:Este e' um programa geral para a soluçao da equacao de Laplace aplicado 
 *	      `a transferencia de calor em um dominio tridimendional.
 *            O metodo utilizado e'o  SOR (variante do metodo Gauss-Seidel) e
 *            pode ser compilado usando duas definicoes diferentes _STAT ou
 *            _SIM e utilizando processamento paralelo. No primeiro caso, 
 *	      a equacao de Laplace e'resolvida variando
 *            o valor do parametro de relaxacao omega (w) e e' computado o numero
 *            de iteracoes para cada valor. O omega em que ocorre o menor numero
 *            de iteracoes e'armazenado em arquivo. A definicao _SIM resolve a
 *            equacao de Laplace utilizando o valor de omega retornado dor _STAT.
 *	      	

 *	  
 *
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 0.00001
//#define N 10	
//#define NX 10
//#define NY 10
//#define NZ 10
//#define H 0.05
#define CMD printf("%s\n",comand);
#define DEBUG printf("\nDEBUG\n");

/******* Definicoes para reduzir o codigo ********/
#define Tn(i,j,k) grid[(i)*NY*NZ + (j)*NZ + (k)] 
#define T(i,j,k) tmp[(i)*NY*NZ + (j)*NZ + (k)]
#define B(i,j,k) b[(i)*NY*NZ + (j)*NZ + (k)]
#define for_(i,n) for(i=0; i<n; i++)

//#define HEAT3D
//#define _STAT 		
#define _SIM
#define PARALLEL

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

int main (int argc, char *argv[])
{

#ifdef PARALLEL
	int id, procs;
	unsigned long int size_msg;
	MPI_Status status;
#endif //PARALLEL

int xi, yi, zi, xf, yf, zf;
int t0=0, t1=0, partial_time=0, total_time=0;
double start=0.0, end, time = 0.0;
unsigned i, j, k, iter, nx, x1, x2, NX, NY, NZ,N;
unsigned long int ntot;
//float w=1.91, wi = 1.00, wf = 2.00, dw = 0.1;
float w, wi, wf, dw;

double Ap, An, As, Ae, Aw, At, Ab;
double *grid, *tmp, *b, *px, *py, *pz;
//double alpha, rho = 2.5, c = 0.092, K = 0.093, t=0.0, tf=10.0, dt=2.0;
//double dx=0.05, dy=0.05, dz=0.05;
double alpha, rho, c, K, t, tf, dt;
double dx, dy, dz;


double  maxdelta, maximum, delta;
FILE *init;
	
#ifdef _STAT
	FILE *arq, *f_w;
	char op, ch;
	int  best_iter=32000, best_time=0;
	double best_w=0.0, b_delta=0.0;
#endif //_STAT       

#ifdef _SIM
	FILE *data, *r_w;
	char op, str[2], count = '0';
	int ctrl = 0;
#endif  //_SIM


// if(id == 0){
 //       char file[20], filename[20], dir[20], dir_w[20],comand[20];

//        char format[6];
 //       char f_omega[6] = "omega";
 //       }
//#else

        char file[20], filename[20], dir[30], dir_w[30], Dir[30], _dir[30], comand[20], comand1[10];
        char format[6];
        char f_omega[6] = "omega";




/**************************** INICIALIZACOES DO MPI ******************************/
#ifdef PARALLEL
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

/*   if(id == 0){
	char file[20], filename[20], dir[20], dir_w[20],comand[20];
	
	char format[6];
	char f_omega[6] = "omega";
	}
#else

	char file[20], filename[20], dir[30], dir_w[30], Dir[30], _dir[30], comand[20], comand1[10];
	char format[6];
	char f_omega[6] = "omega";
*/
#endif /*PARALLEL*/	 

/*************************** INICIALIZACOES *************************************/

init = open_file("init_data","r");
fscanf(init,"%*s%f %*s%f %*s%f %*s%f %*s%lf %*s%lf %*s%lf %*s%lf %*s%lf %*s%lf",&w,&wi,&wf,&dw,&rho,&c,&K,&t,&tf,&dt);
fscanf(init,"%*s%d%*s%lf %*s%d %*s%lf %*s%d %*s%lf %*s%d",&N,&dx,&NX,&dy,&NY,&dz,&NZ);


if(fclose(init))
	printf("\nERRO: Arquivo de inicializacao nao foi fechado\n");

ntot = NX*NY*NZ;

//dx = (double)NX/N;
//dy = (double)NY/N;
//dz = (double)NZ/N;

Aw = Ae = 1.00/(dx*dx);
An = As = 1.00/(dy*dy);
At = Ab = 1.00/(dz*dz);
Ap = An + As + Ae + Aw + At + Ab;

alpha = rho*c/K*dt;

Ap = Ap + alpha;

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
#endif //HEAT3D
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
#endif //HEAT3D
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
#endif //PARALLEL
	
	
	printf("\nomega    iteracoes    epsilon    tempo \n");

	fprintf(arq,"\nRESULTADOS\ndim %d %d %d\n",NX,NY,NZ);
#ifdef PARALLEL
        fprintf(arq,"Numero de processadores: %d\n", procs);
        fprintf(arq,"Cada processador calculou %d nos\n", (x2-x1+2)*NY*NZ);              
#endif //PARALLEL
	fprintf(arq,"\nomega    iteracoes    epsilon    tempo \n");

#ifdef PARALLEL	
}
#endif //PARALLEL
#ifdef PARALLEL
        printf("\nNumero de processadores: %d\n", procs);
	printf("Processador <%d> processou %d nos\n", id, (x2-x1+2)*NY*NZ);
#endif //PARALLEL
		
/************************ ALOCACAO DINAMICA DE MEMORIA ***************************/

	grid = malloc(sizeof(double)*ntot);
	tmp = malloc(sizeof(double)*ntot);
	b = malloc(sizeof(double)*ntot);
#ifdef PARALLEL
   if(id == 0){
#endif //PARALLEL	   
	px = malloc(sizeof(double)*N);
	py = malloc(sizeof(double)*N);
	pz = malloc(sizeof(double)*N);
#ifdef PARALLEL
   }
#endif //PARALLEL

for(w = wi; w < wf; w += dw) {
   iter = 0;
  	
#endif /*_STAT*/

#ifdef _SIM

#ifdef PARALLEL
   if(id == 0){
#endif //PARALLEL
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
	str[0] = count;
	str[1]='\0';
        strcat(filename,str);	
       	strcat(filename,format);
        strcat(dir,"/simulation/");

        strcat(dir,filename);
	count++;

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
	
	
//	printf("\nProcessando ...\n");
	printf("\nRESULTADOS\n");
#ifdef PARALLEL
	printf("Programa rodando sob %d processo(s)\n", procs);
#endif
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
	px = malloc(sizeof(double)*N);
	py = malloc(sizeof(double)*N);
	pz = malloc(sizeof(double)*N);
#ifdef PARALLEL
   }
#endif   

/******************************* INICIALIZACOES **********************************/
#ifdef PARALLEL
      if(id == 0){
#endif
//Gera grid
	 /* for_(i,NX) 
	     px[i] = i*dx;
	  for_(i,NY) 
	      py[i] = i*dy; 
          for_(i,NZ) 
              pz[i] = i*dz; 
	*/

	for_(i,N){
	   px[i] = i*dx;
	   py[i] = i*dy;
	   pz[i] = i*dz;		
	}	      
#ifdef PARALLEL
   }
#endif      
	iter = 0;
#endif  /*_SIM*/

    for_(i,N)
       for_(j,N)	    
	  for_(k,N)
             Tn(i,j,k) = T(i,j,k)=0.0;	  

    
/************************** Condicoes de contorno ********************************/
   xi = yi = zi = N/2 - 0.1*N;
   xf = yf = zf = N/2 + 0.1*N;

    for_(i,N)
       for_(j,N)
          for_(k,N){
             Tn(0,j,k)=T(0,j,k)=10;
	     Tn(N-1,j,k) = T(N-1,j,k) = -10; 	    
        //     if( (i>=xi) && (i<=xf) && (j>=yi) && (j<=yf) && (k>=zi) && (k<=zf))
          //   Tn(i,j,k) = T(i,j,k)= 300;
    	  }

					 
/****************************** CALCULOS *****************************************/
#ifdef _SIM
for(; t<=tf; t+=dt){    
#endif
#ifdef PARALLEL
  start=MPI_Wtime();
#endif



 for_(i,N)
    for_(j,N)
       for_(k,N)
           B(i,j,k)=alpha*Tn(i,j,k);
 
  
    for (;;) {
	maxdelta=0.0;	       
        for (i=x1; i<=x2; i++)
            for (j=1; j <= N-2; j++)
               for (k=1; k<=N-2; k++) {
//		 if( (i>=xi) && (i<=xf) && (j>=yi) && (j<=yf) && (k>=zi) && (k<=zf))
//			 continue;
		       
	         T(i,j,k) = ( Aw*Tn(i-1,j,k) + Ae*Tn(i+1,j,k) +
  			      As*Tn(i,j-1,k) + An*Tn(i,j+1,k) + 
			      Ab*Tn(i,j,k-1) + At*Tn(i,j,k+1) + B(i,j,k))*w/Ap +
	            	      (1.0 - w)*T(i,j,k);
		  delta = T(i,j,k) - Tn(i,j,k);
		  if(delta > maxdelta)
			   maxdelta = fabs(delta);
		  Tn(i,j,k) = T(i,j,k);
                 } 
    
/**************************** PASSAGEM DE MENSAGENS ******************************/ 
#ifdef PARALLEL
       if (id != procs-1){
          MPI_Send(&grid[x2*size_msg], size_msg, MPI_DOUBLE, id+1, 100, MPI_COMM_WORLD);
          MPI_Recv(&grid[(x2+1)*size_msg], size_msg, MPI_DOUBLE, id+1, 200, MPI_COMM_WORLD,&status);
		}	

        if (id != 0) {
          MPI_Recv(&grid[(x1-1)*size_msg], size_msg, MPI_DOUBLE, id-1, 100, MPI_COMM_WORLD, &status);
	  MPI_Send(&grid[x1*size_msg], size_msg, MPI_DOUBLE, id-1, 200, MPI_COMM_WORLD);
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
    
// 	t=t + dt;
	
#ifdef PARALLEL
  if (id != 0)
        MPI_Send(&grid[x1*size_msg], ((x2-x1)+2)*size_msg, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);     

  else {for(i=1; i< procs; i++)
        MPI_Recv(&grid[i*(x2+1)*size_msg], ((x2-x1)+2)*size_msg, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status);
   }
#endif

#ifdef PARALLEL
   end = MPI_Wtime();
   time = end -time;	
   if(id == 0){
#endif

#ifdef PARALLEL
   printf("%5.3f%13d%11.6f%8.3lf\n",w, iter, maximum, time);
#else 
   t = t1-t0;
   printf("%5.3f%13d%11.6f%8.3lf\n",w, iter, maximum, t);
#endif
   
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
	printf("Numero de processadores %d\n", procs);
#endif
        printf("dim %d %d %d\nomega    iteracoes    epsilon    tempo	 ", NX, NY, NZ);
	printf("\n%5.3f%13d%11.6f%8ds\n",best_w, best_iter, b_delta, best_time);
	fprintf(f_w,"%f",best_w);		 
	fclose(arq);
	fclose(f_w);
#endif /*_STAT*/      	       

/********************* GERA ARQUIVO EM FORMATO VTK **********************/
#ifdef _SIM

	printf("Gerando arquivo vtk ...");
	fprint_header_vtk(data, "teste vtk", "ASCII", "RECTILINEAR_GRID");
	fprintf(data,"DIMENSIONS %d %d %d\n", NX, NY, NZ);
	//fprintf(data, "POINTS %d float\n", N*N*N);
	
	fprintf(data, "\nX_COORDINATES %d float\n", NX);
	fprint_vector(data, px, N);
	fprintf(data, "\nY_COORDINATES %d float\n", NY);
	fprint_vector(data, py, N);
	fprintf(data, "\nZ_COORDINATES %d float\n", NZ);
	fprint_vector(data, pz, N);
	
	fprintf(data, "\n\nPOINT_DATA %d\n", ntot);
	fprintf(data, "SCALARS temperatura float 1\n");
	fprintf(data, "LOOKUP_TABLE default\n");
//	fprintf(data, "\nFIELD escalar 1\n");
//      fprintf(data, "temperatura 1 %d 1\n",ntot);
	fprint_mat(data,grid,N,N,N);
        printf("    [ OK ]\n"); 
	fclose(data);
	if(ctrl)
	fclose(r_w);
       // printf("\nt=%f\n",t);
	if(t<tf)
	{  
		strcpy(filename,file);
        	str[0] = count;
		strcat(filename,str);
		strcat(filename,format);
		strcpy(_dir,Dir);
		strcat(_dir,"/simulation/");
  		strcat(_dir,filename);
   	        count++;
	        data = open_file(_dir,"w");
		iter=0;
	}
     
   }
#endif
//fclose(init);
#ifdef PARALLEL
	}
	MPI_Finalize();
#endif

return 0; 
}
		    
 /************************************* FUNCOES  ****************************************/

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
	for(i=0; i<m; i++){ 
	   for(j=0; j<n; j++){
	      for(k=0; k<p; k++)
	         fprintf(fp,"%8.4f ",mat[i*n*p + j*p + k]);
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
				fprintf(fpoint,"%f %f %f \n", px[i], py[j], pz[k]);
			
}
void fprint_vector(FILE *fpoint, double *p, unsigned n)
{
        unsigned i;
        for(i=0; i<n; i++)
           fprintf(fpoint,"%f ", p[i]);

}
