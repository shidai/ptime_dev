#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include "ptimeLib.h"
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2.
 *
 *    TEMPO2 is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    TEMPO2 is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    You should have received a copy of the GNU General Public License
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2,
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "TKsvd.h"
#include "TKmatrix.h"
#include "T2toolkit.h"
#include "TKfit.h"

void TKremovePoly_f(float *px,float *py,int n,int m)
{
   int i,j;
   double x[n],y[n];
   double p[m];
   double chisq;
   double v[m];

   for (i=0;i<n;i++)
   {
	  x[i] = (float)px[i];
	  y[i] = (float)py[i];
   }
   TKleastSquares_svd_noErr(x,y,n, p, m, TKfitPoly);  

   for (i=0;i<n;i++)
   {
	  TKfitPoly(x[i],v,m);
	  for (j=0;j<m;j++)
		 py[i] -= v[j]*p[j];
   }
}

void TKremovePoly_d(double *x,double *y,int n,int m)
{
   int i,j;
   double p[m];
   double chisq;
   double v[m];

   printf("Remove polynomial n=%d m=%d\n",n,m);
   TKleastSquares_svd_noErr(x,y,n, p, m, TKfitPoly);  

   for (i=0;i<n;i++)
   {
	  TKfitPoly(x[i],v,m);
	  for (j=0;j<m;j++)
		 y[i] -= v[j]*p[j];
   }
}

void TKfitPoly(double x,double *v,int m)
{
   int i;
   double t=1;
   for (i=0;i<m;i++)
   {
	  v[i] = t;
	  t*=x;
   }
}

/* Least squares fitting routines */


// returns the chisq.
double TKleastSquares(double* b, double* white_b,
	  double** designMatrix, double** white_designMatrix,
	  int n,int nf, double tol, char rescale_errors,
	  double* outP, double* e, double** cvm){

   double chisq = 0;
   int i,j,k;

   // quad precision arrays for fitting.
   long double svd_V[n];
   long double **svd_M=malloc_2dLL(n,nf);
   long double **v=malloc_2dLL(nf,nf);
   long double **u=malloc_2dLL(n,nf);
   long double w[nf],wt[nf];
   long double p[nf];

   // other variables
   long double sum,wmax,sum_w;
   bool computeErrors = (e!=NULL);
   bool computeCVM = (cvm!=NULL);
   bool computeParam = (outP!=NULL && b!=NULL);
   bool needToFreeCVM=false;
   int writeResiduals = 0;

   if(computeErrors && ! computeCVM){
	  // we can't easily compute the errors without the CVM matrix
	  // so we have to create one.
	  cvm=malloc_uinv(nf);
	  computeCVM=true;
	  needToFreeCVM=true;
   }

   if(writeResiduals==1 && white_b!=NULL && b!=NULL){
	  FILE* wFile=fopen("prefit.res","w");
	  if (!wFile){
	    printf("Unable to write out whitened residuals: cannot open file prefit.res\n");
	  }
	  else
	    {
	      for (i=0;i<n;i++)
		fprintf(wFile,"%d %lg %lg\n",i,b[i],white_b[i]);
	      fclose(wFile);
	    }

   }
   if(writeResiduals==1){
	  printf("Writing out design matrix\n");
	  FILE * wFile=fopen("design.matrix","w");
	  if (!wFile){
	    printf("Unable to write out design matrix: cannot open file design.matrix\n");
	  }
	  else
	    {
	      for (i=0;i<n;i++) {
		for (j=0;j<nf;j++){
		  fprintf(wFile,"%d %d %lg %lg\n",i,j,designMatrix[i][j],white_designMatrix[i][j]);
		}
		fprintf(wFile,"\n");
	      }
	      fclose(wFile);
	    }
   }

   // Now go to long double precision!
   for (i=0;i<n;i++){
     for (j=0;j<nf;j++) svd_M[i][j] = white_designMatrix[i][j];
   }

   /* Now carry out the singular value decomposition */
   // note that this modifies svd_M
   TKsingularValueDecomposition_lsq(svd_M,n,nf,v,w,u);

   wmax = TKfindMax_Ld(w,nf);
   long double sensible_wmax=pow(2,sizeof(long double)*8-17);
   if (wmax > sensible_wmax){
	  printf("Warning: wmax very large. Precision issues likely to break fit\nwmax=%Lf\ngood=%Lf",wmax,sensible_wmax);
   }

   for (i=0;i<nf;i++)
   {
	  if (w[i] < tol*wmax) w[i]=0.0;
   }
   /* Back substitution */


   /* Now form the covariance matrix */
   if(computeCVM){
	  for (i=0;i<nf;i++)
	  {
		 if (w[i]!=0) wt[i] = 1.0/w[i]/w[i];
		 else wt[i] = 0.0;     
	  }
	  for (i=0;i<nf;i++)
	  {
		 for (j=0;j<=i;j++)
		 {
			sum=0.0;
			for (k=0;k<nf;k++)
			   sum+=v[i][k]*v[j][k]*wt[k];
			cvm[i][j] = cvm[j][i] = (double)sum;
		 }
	  } 

	  //	  if(debugFlag==1) 
	    {
		 FILE *fout;
		 fout = fopen("cvm.matrix","w");
		 if (!fout){
		   printf("Unable to open file cvm.matrix for writing\n");
		 }
		 else{
		   for (i=0;i<nf;i++)
		     {
		       for (j=0;j<=i;j++)
			 {
			   fprintf(fout,"%+.8f ",cvm[i][j]/sqrt(cvm[i][i]*cvm[j][j]));
			 }
		       fprintf(fout,"\n");
		     }
		   fclose(fout);
		 }
	  }
	  if(computeErrors){
	  for (i=0;i<nf;i++){e[i]=sqrt(cvm[i][i]);}
	  }

   }


   if (computeParam){

     for (i=0;i<n;i++){
       svd_V[i] = white_b[i];
       //       printf("svd_V = %g\n",(double)svd_V[i]);
     }

     TKbacksubstitution_svd(v, w, svd_M, svd_V, p, n, nf);
     
     for (k=0;k<nf;k++)outP[k]=(double)(p[k]);
     // compute chisq
     chisq = 0.0;
     for (j=0;j<n;j++)
       {
	 sum = 0.0;
	 for (k=0;k<nf;k++)
	   {
	     sum+=p[k]*white_designMatrix[j][k];
	   }
	 chisq += pow((white_b[j]-sum),2);
       }
     
     if(computeErrors && rescale_errors){
       for (j=0;j<nf;j++)
	 e[j] *= sqrt(chisq/(n-nf));
     }
     
     if (writeResiduals){
       FILE* wFile=fopen("postfit.res","w");
       if (!wFile){
	 printf("Unable to open file postfit.res for writing\n");
       }
       else
	 {
	   for (i=0;i<n;i++)
	     {
	       sum=0;
	       sum_w=0;
	       for (j=0;j<nf;j++){
		 sum += designMatrix[i][j]*p[j];
		 sum_w += white_designMatrix[i][j]*p[j];
	       }
	       fprintf(wFile,"%d %lg %lg\n",i,(double)(b[i]-sum),(double)(white_b[i]-sum_w));
	     }
	   fclose(wFile);
	 }
     }
   } // computeParam
   // this funny method of freeing is because of the BLAS style matricies. M.Keith 2012
   free_2dLL(v);     // free-TKleastSquares_svd_psr_dcm-v**
   free_2dLL(u);     // free-TKleastSquares_svd_psr_dcm-u**
   free_2dLL(svd_M);

   if(needToFreeCVM){
	  // we created CVM, so free it
	  free_uinv(cvm);
   }
   return chisq;
}


// legacy method.




void TKleastSquares_svd_noErr(double *x,double *y,int n,double *p,int nf, void (*fitFuncs)(double, double [], int))
{
   double chisq=0;
   //   TKleastSquares_svd(x,y,NULL,n,p,NULL,nf,NULL,&chisq,fitFuncs,0);
}

void TKfit_getDesignMatrix(double *x,double *y,int n,int nf,void (*fitFuncs)(double, double [], int,int), double **uinv,double ***OUT_designMatrix,double ***OUT_white_designMatrix,double** OUT_b, double** OUT_wb){

   //double precision arrays for matrix algebra.
   double **designMatrix, **white_designMatrix;
   double basisFunc[nf];
   double *b,*white_b;
   int    i,j,k;
   int nrows=get_blas_rows(uinv);
   int ncols=get_blas_cols(uinv);
   if (ncols!=n){
	  printf("n=%d ncols=%d\n",n,ncols);
	  printf("uinv error. Either you did not use malloc_uinv() to create uinv or np!=ncols\n");
	  exit(1);
   }

   if (nrows!=n && nrows != 1){
	  printf("n=%d nrows=%d\n",n,nrows);
	  printf("uinv error. Either you did not use malloc_uinv() to create uinv or np!=nrows\n");
	  exit(1);
   }


   // double arrays
   white_designMatrix=malloc_blas(n,nf);
   designMatrix=malloc_blas(n,nf);
   b=(double*)malloc(sizeof(double)*n);
   white_b=(double*)malloc(sizeof(double)*n);

   /* This routine has been developed from Section 15 in Numerical Recipes */

   /* Determine the design matrix - eq 15.4.4 
	* and the vector 'b' - eq 15.4.5 
	*/
   for (i=0;i<n;i++)
   {
	  // fitFuncs is not threadsafe!
	  fitFuncs(x[i],basisFunc,nf,i);
	  for (j=0;j<nf;j++) designMatrix[i][j] = basisFunc[j];
	  b[i] = y[i];
   }
   // Take into account the data covariance matrix
   printf("Got to this bit\n");
   if(nrows==1){
	  // we have only diagonal elements
	  for (i=0;i<n;i++){
		 white_b[i]=b[i]*uinv[0][i];
		 for (j=0;j<nf;j++){
			white_designMatrix[i][j] = designMatrix[i][j]*uinv[0][i];
		 }
	  }
   } else {
	  TKmultMatrix_sq(uinv,designMatrix,n,nf,white_designMatrix);  
	  TKmultMatrixVec_sq(uinv,b,n,white_b);
   }

   /*   for (i=0;i<n;i++)
     {
       for (j=0;j<nf;j++)
	 printf("%d %d %g\n",i,j,designMatrix[i][j]);
	 }*/

   *OUT_designMatrix=designMatrix;
   *OUT_white_designMatrix=white_designMatrix;
   *OUT_b=b;
   *OUT_wb=white_b;
}



// Non-pulsar fit. No cholesky yet though...
void TKleastSquares_svd(double *x,double *y,double *sig,int *vali,int *valj,int *valk,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], tmplStruct *tmpl, int, double,int,int,int),tmplStruct *tmpl, int weight,double phi)
{
   double **designMatrix, **white_designMatrix;
   double basisFunc[nf];
   double *b,*white_b;
   int    i,j,k;

   // double arrays
   white_designMatrix=malloc_blas(n,nf);
   designMatrix=malloc_blas(n,nf);
   b=(double*)malloc(sizeof(double)*n);
   white_b=(double*)malloc(sizeof(double)*n);

   //      printf("Non pulsar least-squares fit. n=%d nf=%d\n",n,nf);
   /* Determine the design matrix - eq 15.4.4 
	* and the vector 'b' - eq 15.4.5 
	*/
   for (i=0;i<n;i++)
   {
	  // fitFuncs is not threadsafe!
	  fitFuncs(x[i],basisFunc,tmpl,nf,phi,vali[i],valj[i],valk[i]);
	  for (j=0;j<nf;j++) designMatrix[i][j] = basisFunc[j];
	  b[i] = y[i];
   }
   //   printf("Complete design matrix\n");
   // deal with the weights if we are doing a weighted fit.
   if(weight==1 && sig!=NULL){
     //	  printf("Divide by errors\n");
	  for (i=0;i<n;i++){
		 white_b[i]=b[i]/sig[i];
		 for (j=0;j<nf;j++) white_designMatrix[i][j] = designMatrix[i][j]/sig[i];
	  }
   } else {
	  for (i=0;i<n;i++){
		 white_b[i]=b[i];
		 for (j=0;j<nf;j++) white_designMatrix[i][j] = designMatrix[i][j];
	  }

   }

   // go ahead and do the fit!

   *chisq = TKleastSquares(b,white_b,designMatrix,white_designMatrix,
		 n,nf,1e-10,1,
		 p,e,cvm);

   free_blas(designMatrix); // free-TKleastSquares_svd_psr_dcm-designMatrix**
   free_blas(white_designMatrix);  // free-TKleastSquares_svd_psr_dcm-white_designMatrix**
   free(b);
   free(white_b);


}

void TKleastSquares_svd_passN(double *x,double *y,double *sig2,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int,int),int weight)
{
   printf("This method no longer implemented.\n");
   exit(-1);
}

long double TKfindMax_Ld(long double *x,int n)
{
   long double ret;
   int i;

   ret = x[0];
   for (i=0;i<n;i++)
   {
	  if (x[i] > ret) ret = x[i];
   }
   return ret;
}

