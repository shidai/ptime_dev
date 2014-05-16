//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

#ifndef __TKfit_h
#define __TKfit_h
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


#include "TKmatrix.h"
#include "ptimeLib.h"

void TKfit_getDesignMatrix(double *x,double *y,int n,int nf,void (*fitFuncs)(double, double [], int,int), double **uinv,double ***OUT_designMatrix,double ***OUT_white_designMatrix,double** OUT_b, double** OUT_wb);

double TKleastSquares(double* b, double* white_b,
	  double** designMatrix, double** white_designMatrix,
	  int n,int nf, double tol, char rescale_errors,
	  double* outP, double* e, double** CVM);

long double TKfindMax_Ld(long double *x,int n);
void TKremovePoly_f(float *px,float *py,int n,int m);
void TKremovePoly_d(double *px,double *py,int n,int m);
void TKfitPoly(double x,double *v,int m);


// legacy methods
void TKleastSquares_svd(double *x,double *y,double *sig,int *vali,int *valj,int *valk,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], tmplStruct *tmpl,int,double,int,int,int),tmplStruct *tmpl, int weight,double phi);
void TKleastSquares_svd_noErr(double *x,double *y,int n,double *p,int nf, void (*fitFuncs)(double, double [], int));     

#endif
