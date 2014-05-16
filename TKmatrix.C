#include <stdlib.h>
#include <stdio.h>
#include "TKmatrix.h"


void TKmultMatrix(double **idcm,double **u,int ndata,int ndata2,int npol,double **uout)
{
#ifdef ACCEL_MULTMATRIX
   if(useT2accel){
     //      printf("Using ACCELERATED multmatrix (M.Keith 2012)\n");
      accel_multMatrix(idcm[0],u[0],ndata,ndata2,npol,uout[0]);
   } else {
#endif
     //      printf("Using SLOW multmatrix\n");
      int i,j,k;
      for (j=0;j<npol;j++)
      {
         for (i=0;i<ndata;i++)
         {
            uout[i][j]=0.0;
         }
      }
      for (k=0;k<ndata2;k++) {
         for (i=0;i<ndata;i++){
            for (j=0;j<npol;j++){
               uout[i][j]+=idcm[i][k]*u[k][j];
            }
         }
      }

#ifdef ACCEL_MULTMATRIX
   }
#endif
}


void TKmultMatrix_sq(double **idcm,double **u,int ndata,int npol,double **uout)
{
   TKmultMatrix(idcm,u,ndata,ndata,npol,uout);
}
void TKmultMatrixVec(double **idcm,double *b,int ndata,int ndata2,double *bout)
{
#ifdef ACCEL_MULTMATRIX
   if (useT2accel){
     //      printf("using accelerated multmatrix (m.keith 2012)\n");
      accel_multMatrix(idcm[0],b,ndata,ndata2,1,bout);
   }else{
#endif
      int i,j;
      for (i=0;i<ndata;i++)
      {
         bout[i] = 0.0;
      }

      for (j=0;j<ndata2;j++)
      {
         for (i=0;i<ndata;i++)
            bout[i]+=idcm[i][j]*b[j];
      }

#ifdef ACCEL_MULTMATRIX
   }
#endif
}


void TKmultMatrixVec_sq(double **idcm,double *b,int ndata,double *bout)
{
   TKmultMatrixVec(idcm,b,ndata,ndata,bout);
}



long double** malloc_2dLL(int rows,int cols){
   int i;
   long double* memory;
   long double** m;
   //   printf("Allocate %d x %d long double array (%.3f kb)\n",rows,cols, (double)(rows*cols*sizeof(long double)/1024.0));

   m=(long double**) calloc(rows,sizeof(long double*));
   memory=(long double*)malloc(sizeof(long double)*rows*cols);

   if(memory==NULL || m==NULL){
	  printf("Could not allocate %d x %d long double array (%.3f kb)\n",rows,cols, (double)(rows*cols*sizeof(double)/1024.0));
	  printf("Cannot allocate enough memory for array\n");
	  exit(1);
   }
   //   printf("Allocated mem=0x%016x m=0x%016x\n",memory,m);
   for(i=0;i<rows;i++){
	  m[i]=memory+cols*i;
   }

   //   printf("Accessible memory 0x%016x -> 0x%016x\n",memory,memory+rows*cols);

   return m;

}

void free_2dLL(long double** m){
   free(m[0]);
   free(m);
}


/**
 * Allocate uinv in a "BLAS/LAPACK" compatile way
 */
double **malloc_uinv(int n){
   return malloc_blas(n,n);
}

int get_blas_rows(double** uinv){
   int* dim=(int*)(uinv[0]-2);
   return *dim;
}

int get_blas_cols(double** uinv){
   int* dim=(int*)(uinv[0]-1);
   return *dim;
}
/**
 * Allocate uinv in a "BLAS/LAPACK" compatile way
 * store the dimensions of the array in two secret elements before
 * the main memory allocation. Useful for checks.
 * WARNING: assumes that sizeof(int) <= sizeof(double)
 */
double **malloc_blas(int rows,int cols){
   double *memory;
   int *dimN;
   int *dimM;
   double** uinv;
   int i;
   if (sizeof(int) > sizeof(double)){
      printf("Error, somehow you have a system with sizeof(int) > sizeof(double) %d %d\n",sizeof(int),sizeof(double));
      exit(1);
   }
   //   printf("Allocate %d x %d double array (%.3f kb)\n",rows,cols, (double)(rows*cols*sizeof(double)/1024.0));
   memory=(double*)calloc((rows*cols+2),sizeof(double)); // we add 2 to store dimensions
   uinv=(double**)malloc(sizeof(double*)*rows);
   if(memory==NULL || uinv==NULL){
   printf("Could not allocate %d x %d double array (%.3f kb)\n",rows,cols, (double)(rows*cols*sizeof(double)/1024.0));
      printf("Cannot allocate enough memory for array\n");
      exit(1);
   }
   //   printf("Allocated mem=0x%016x uinv=0x%016x\n",memory,uinv);
   memory+=2; // the first two bytes are for the dimensions.
   dimN=(int*)(memory-2);
   dimM=(int*)(memory-1);

   *dimN = rows;
   *dimM = cols;

   //   printf("Secret memory rows=0x%016x cols=0x%016x\n",dimN,dimM);

   //   printf("Accessible memory 0x%016x -> 0x%016x\n",memory,memory+rows*cols);
   for(i=0;i<rows;i++){
      uinv[i]=memory+cols*i;
   }
   return uinv;
}

void free_blas(double** m){
  //   if(debugFlag){
  //	  printf("free 0x%016x\n",m[0]-2);
  //	  printf("m was %d x %d\n",get_blas_rows(m),get_blas_cols(m));
	  //   }
   fflush(stdout);
   free(m[0]-2);
   //   printf("free 0x%016x\n",m);
   fflush(stdout);
   free(m);
   //   printf("leaving free_blas\n");
      fflush(stdout);
}


void free_uinv(double** uinv){
   free_blas(uinv);
}



