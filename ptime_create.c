// To do:
// 7. Text output showing what you've done
// 8. Allow to delete some components
// 9. Global fitting in frequency
// 11. Properly fill up the date/name etc.
// 12. Zoom in y-direction on "difference" and on "profile"
// 13. Move fitting to different library
// 15. Allow postscript output of plot and summary of template
// 17. Have some method for labelling components
// 18. When have a good model for one channel then auto fit to all other channels
// 19. Properly get the frequecny low and high set (done .. but needs checking)
// 22. Deal correctly if file has been dedispersed
// 23. Currently must use 'pam -S' to convert to Stokes parameters
// 24. Be able to deal with coherence parameters instead of Stokes
// 25. Properly free memory

// Program to produce Stokes and Frequency-dependent analytic templates
// This program can load a PSRFITS file that has been completely time scrunched (nsub = 1)
// It can process multiple polarisation states and multiple frequency channels
// G. Hobbs (V1.0 Dec 25th 2013)
// S. Dai (V1.01 May 14th 2014)
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "ptime.h"
#include "readPfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#define MAX_COMPONENTS 64 // Maximum number of components per polarisation per frequency channel
#define MAX_VONMISES 16 // Maximum number of Von Mises functions per component

int global_nComp;

typedef struct ptime_pol {
  int nbin;
  float *val;
  double devi;
} ptime_pol;

typedef struct ptime_chan {
  int npol;
  ptime_pol *pol;
} ptime_chan;

typedef struct ptime_observation {
  int nchan;
  ptime_chan *chan; 
} ptime_observation;

struct data {
  size_t n;
  int nComp; // number of components
  int *nVm; // number of Von Mises function for each component
  double *x;
  double *y;
  double *sigma;
};

struct datag {
  size_t n;
  int nComp; // number of components
  int *nVm; // number of Von Mises function for each component
  double *x;
  double *y;
  double *sigma;
	double *par;
};

double nonlinearFunc( double x, const double *par);
void allocateObsMemory(ptime_observation *obs,pheader *phead);
void deallocateMemory(ptime_observation *obs);
void readData(ptime_observation *obs,pheader *phead,fitsfile *fp);
void plot(ptime_observation *obs,pheader *phead,tmplStruct *tmpl);
void findMinMax(float *y,int n,float *min,float *max);
void removeBaseline(ptime_observation *obs,pheader *phead,int baselineType,float baselineFrac);
void startNewTemplate(tmplStruct *tmpl,pheader *phead);
void doFit(float *fx,float *fy,int nbin,tmplStruct *tmpl,int chan,int pol);

int Mises_f (const gsl_vector *x, void *data, gsl_vector *f);
int Mises_df (const gsl_vector *x, void *data, gsl_matrix *J);
int Mises_fdf (const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J);
void doFit_err(float *fx,float *fy, double devi, int nbin,tmplStruct *tmpl,int chan,int pol);
void print_state (size_t iter, gsl_multifit_fdfsolver * s);
void doFit_gerr(float *fx,float *fy, double devi, int nbin,tmplStruct *tmpl,int chan,int pol);
void doFit_Gerr(float *fx,float *fy, double devi, int nbin,tmplStruct *tmpl,int chan,int pol);
void updateTemplateHeader(tmplStruct *tmpl,pheader *phead);

int main(int argc,char *argv[])
{
  int baselineType = 1;
  float baselineFrac = 0.1;
  char fname[128];
  tmplStruct tmpl;
  int i;
  pheader *header;
  fitsfile *fp;
  ptime_observation *obs;
  char units[128] = "Unset";
  int setUnits=0;

  header = (pheader *)malloc(sizeof(pheader));
  obs = (ptime_observation *)malloc(sizeof(ptime_observation));
  
  initialiseTemplate(&tmpl);
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	sprintf(fname,argv[++i]);
      else if (strcmp(argv[i],"-u")==0)
	{strcpy(units,argv[++i]); setUnits = 1;}
    }
  fp = openFitsFile(fname);
	printf ("Open fits file done.\n");
  loadPrimaryHeader(fp,header);
  allocateMemoryTemplateDefault(&tmpl,header->nchan,header->npol,MAX_COMPONENTS,MAX_VONMISES);
  startNewTemplate(&tmpl,header);
  if (setUnits==1) strcpy(tmpl.units,units);
  strcpy(tmpl.profileFile,fname);
  allocateObsMemory(obs,header);
  readData(obs,header,fp);
  closeFitsFile(fp);

  // Must remove a baseline
  removeBaseline(obs,header,baselineType,baselineFrac);

  // Plot the data
  plot(obs,header,&tmpl);

  // Should deallocate the template memory
  
  // Should deallocate the observation memory
  free(header);
  deallocateMemory(obs);
}

void startNewTemplate(tmplStruct *tmpl,pheader *phead)
{
  int npol = phead->npol;
  int nchan = phead->nchan;
  int i,j;
  double f0 = phead->freq - phead->bw/2.0;
  double f1,f2;
  double fa,fb;

  strcpy(tmpl->source,phead->source);
  //printf ("%s\n", tmpl->source);

  tmpl->nchan = nchan;
  // Should set up user/date etc.
  for (i=0;i<nchan;i++)
    {
      tmpl->channel[i].nstokes = npol;
      fa = f0 + i*phead->bw/phead->nchan - (phead->bw/phead->nchan)/2.0;
      fb = f0 + i*phead->bw/phead->nchan + (phead->bw/phead->nchan)/2.0;
      if (fa < fb) {f1=fa; f2=fb;}
      else {f1=fb; f2=fa;}
      tmpl->channel[i].freqLow = f1;
      tmpl->channel[i].freqHigh = f2;
      // Should set up frequency range here
      for (j=0;j<npol;j++)
	{
	  tmpl->channel[i].pol[j].stokes = j;
	  tmpl->channel[i].pol[j].nComp = 0;
	}
    }
}

void updateTemplateHeader(tmplStruct *tmpl,pheader *phead)
{
  int npol = phead->npol;
  int nchan = phead->nchan;
  int i,j;
  double f0 = phead->freq - phead->bw/2.0;
  double f1,f2;
  double fa,fb;

  strcpy(tmpl->source,phead->source);

  tmpl->nchan = nchan;
  // Should set up user/date etc.
  for (i=0;i<nchan;i++)
    {
      tmpl->channel[i].nstokes = npol;
      fa = f0 + i*phead->bw/phead->nchan - (phead->bw/phead->nchan)/2.0;
      fb = f0 + i*phead->bw/phead->nchan + (phead->bw/phead->nchan)/2.0;
      if (fa < fb) {f1=fa; f2=fb;}
      else {f1=fb; f2=fa;}
      tmpl->channel[i].freqLow = f1;
      tmpl->channel[i].freqHigh = f2;
      // Should set up frequency range here
      for (j=0;j<npol;j++)
	{
	  tmpl->channel[i].pol[j].stokes = j;
	}
    }
}

void removeBaseline(ptime_observation *obs,pheader *phead,int baselineType,float baselineFrac)
{
  int nchan,npol,nbin;
  int i,j,k,k0,k1,it;
  double bl,bl_best;
  int setbl=0;
  int best_k0=0;

  double devi;

  nchan = phead->nchan;
  nbin = phead->nbin;
  npol = phead->npol;


  for (i=0;i<nchan;i++)
  {
    for (j=0;j<npol;j++)
    {
      if (j==0) 
      {
	setbl=0;
	for (k0=0;k0<nbin;k0++)
	{
	  for (it = 0; it < baselineFrac*nbin; it++)
	  {
	    bl=0;
	    k1 = k0+it;
	
	    if (k1 > nbin) k1 = k1-nbin;
	
	    if (k1 > k0) 
	    {
	      for (k=k0;k<k1;k++)
		bl += obs->chan[i].pol[j].val[k];
	    }
	    else if (k1 < k0) 
	    {
	      for (k=k0;k<nbin;k++)
		bl += obs->chan[i].pol[j].val[k];
	      for (k=0;k<k1;k++)
		bl += obs->chan[i].pol[j].val[k];
	    }
	  }
	  if (setbl==0) {bl_best = bl; setbl=1;}
	  else
	  {
	    if (bl_best > bl) {bl_best = bl; best_k0=k0;}
	  }
	}
	bl_best = bl_best/(baselineFrac*nbin);

	devi = 0.0;
	for (it = 0; it < baselineFrac*nbin; it++)
	{
	  k1 = best_k0 + it;
	  if (k1 > nbin) k1 = k1-nbin;
		    
	  if (k1 > best_k0) 
	  {
	    for (k = best_k0; k < k1; k++)
	      devi += (obs->chan[i].pol[j].val[k]-bl_best)*(obs->chan[i].pol[j].val[k]-bl_best);
	  }
	  else if (k1 < best_k0) 
	  {
	    for (k = best_k0; k < nbin; k++)
	      devi += (obs->chan[i].pol[j].val[k]-bl_best)*(obs->chan[i].pol[j].val[k]-bl_best);
	    for (k = 0; k < k1; k++)
	      devi += (obs->chan[i].pol[j].val[k]-bl_best)*(obs->chan[i].pol[j].val[k]-bl_best);
	  }
	}
	
	devi = sqrt(devi/(baselineFrac*nbin));
      } 
      else 
      {
	bl = 0;
	for (it = 0; it < baselineFrac*nbin; it++)
	{
	  bl=0;
	  k1 = best_k0+it;
	  if (k1 > nbin) k1 = k1-nbin;
	
	  if (k1 > best_k0) 
	  {
	    for (k=best_k0;k<k1;k++)
	      bl += obs->chan[i].pol[j].val[k];
	  }
	
	  else if (k1 < best_k0) 
	  {
	    for (k=k0;k<nbin;k++)
	      bl += obs->chan[i].pol[j].val[k];
	    for (k=0;k<k1;k++)
	      bl += obs->chan[i].pol[j].val[k];
	  }
	}
	bl_best = bl;
	bl_best = bl_best/(baselineFrac*nbin);

	devi = 0.0;

	for (it = 0; it < baselineFrac*nbin; it++)
	{
	  k1 = best_k0 + it;
	  if (k1 > nbin) k1 = k1-nbin;
	  if (k1 > best_k0) 
	  {
	    for (k = best_k0; k < k1; k++)
	      devi += (obs->chan[i].pol[j].val[k]-bl_best)*(obs->chan[i].pol[j].val[k]-bl_best);
	  }
	  else if (k1 < best_k0) 
	  {
	    for (k = best_k0; k < nbin; k++)
	      devi += (obs->chan[i].pol[j].val[k]-bl_best)*(obs->chan[i].pol[j].val[k]-bl_best);
	    for (k = 0; k < k1; k++)
	      devi += (obs->chan[i].pol[j].val[k]-bl_best)*(obs->chan[i].pol[j].val[k]-bl_best);
	  }
	}
	devi = sqrt(devi/(baselineFrac*nbin));
      }

	
      for (k = 0; k < nbin; k++)
      {
	(obs->chan[i].pol[j].val[k]) -= bl_best;
      }
      (obs->chan[i].pol[j].devi) = devi;
    }
  }
}

void plot(ptime_observation *obs,pheader *phead,tmplStruct *tmpl)
{
  int showComp=-1;
  int nchan,npol,nbin;
  float fx[phead->nbin*2],fy[phead->nbin*2];
  float ty[phead->nbin*2],dy[phead->nbin*2];
  float cy[phead->nbin*2];
  int i,j,k;
  int chan,pol,nComp;
  float minx,maxx,miny,maxy;
  float miny0,maxy0;
  float miny2,maxy2;
  float mx,my;
  char key;
  char title[128];
  int zoom=0;
  int saved=1;
  char ylabel[128];
 
  nchan = phead->nchan;
  nbin = phead->nbin;
  npol = phead->npol;

  chan = 0;
  pol = 0;

  cpgbeg(0,"/xs",1,1);
  cpgsch(1.2);
  cpgscf(2);

  do {
    sprintf(title,"Channel number %d, polarisation %d",chan,pol);
    // 
    nComp = tmpl->channel[chan].pol[pol].nComp;
    printf("Number of components = %d\n",nComp);
    for (i=0;i<nbin;i++)
      {
	fx[i] = (float)i/(float)nbin;
	fy[i] = obs->chan[chan].pol[pol].val[i];
	fx[i+nbin] = (float)(i+nbin)/(float)nbin;
	fy[i+nbin] = fy[i];
	// Evaluate template
	ty[i] = (float)evaluateTemplateChannel(tmpl,fx[i],chan,pol,0);
	ty[i+nbin] = ty[i];
	dy[i] = fy[i] - ty[i];
	dy[i+nbin] = dy[i];
      }
    findMinMax(fy,nbin,&miny0,&maxy0);
    findMinMax(dy,nbin,&miny2,&maxy2);
    miny = miny0-(maxy0-miny0)*0.1;
    maxy = maxy0+(maxy0-miny0)*0.1;
    if (zoom==0)
      {
	minx = 0;
	maxx = 2;
      }
    cpgsvp(0.1,0.9,0.5,0.9);
    cpgeras();
    cpgswin(minx,maxx,miny,maxy);
    cpgbox("BCTS",0.0,0,"BCTSN",0.0,0);
    
    //  cpgenv(minx,maxx,miny,maxy,0,1);
    sprintf(ylabel,"Signal strength (%s)",tmpl->units);
    cpglab("",ylabel,title);
    cpgbin(nbin*2,fx,fy,1);
    cpgsci(2); cpgbin(nbin*2,fx,ty,1); cpgsci(1);

    if (showComp==1)
      {
	for (k=0;k<tmpl->channel[chan].pol[pol].nComp;k++)
	  {
	    for (i=0;i<nbin;i++)
	      {
		cy[i] = (float)evaluateTemplateComponent(tmpl,fx[i],chan,pol,k,0);
		cy[i+nbin] = cy[i];
	      }
	    cpgsls(4); cpgsci(3+k);
	    cpgbin(nbin*2,fx,cy,1);
	    cpgsci(1);

	    cpgsls(1);
	  }
      }

    // Now draw the difference between the profile and the model
    cpgsvp(0.1,0.9,0.15,0.5);
    cpgswin(minx,maxx,miny2,maxy2);
    cpgbox("BCTSN",0.0,0,"BCTSN",0.0,0);
    cpglab("Phase","prof-tmpl","");
    cpgbin(nbin*2,fx,dy,1);
    cpgcurs(&mx,&my,&key);
    if (key=='c'){
      int ochan = chan;
      printf("Current channel = %d (from 0 to %d). Please enter new channel number: ",chan,nchan-1);
      scanf("%d",&chan);
      if (tmpl->channel[chan].pol[pol].nComp == 0)
	{
	  // Copy current components
	  tmpl->channel[chan].pol[pol].nComp = tmpl->channel[ochan].pol[pol].nComp;
	  for (i=0;i<tmpl->channel[chan].pol[pol].nComp;i++)
	    {
	      for (j=0;j<tmpl->channel[chan].pol[pol].comp[i].nVm;j++)
	      {
		tmpl->channel[chan].pol[pol].comp[i].vonMises[j].concentration = tmpl->channel[ochan].pol[pol].comp[i].vonMises[j].concentration;
		tmpl->channel[chan].pol[pol].comp[i].vonMises[j].height = tmpl->channel[ochan].pol[pol].comp[i].vonMises[j].height;
		tmpl->channel[chan].pol[pol].comp[i].vonMises[j].centroid = tmpl->channel[ochan].pol[pol].comp[i].vonMises[j].centroid;
	      }
	    }
	}
    }
    else if (key=='l'){ // Load
      char fname[128];
      printf("Please enter filename to load: ");
      scanf("%s",fname);
      printf("Trying to load\n");
      readTemplate(fname,tmpl);
			updateTemplateHeader(tmpl,phead);
      printf("Completed load\n");
    }
    else if (key=='p'){
      int opol = pol;
      printf("Current polarisation = %d (from 0 to %d). Please enter new channel number: ",pol,npol-1);
      scanf("%d",&pol);
    }
    else if (key=='r')
      {
	char yesno[128];
	printf("Clear the template components for this channel/stokes? (y/n) ");
	scanf("%s",yesno);
	if (yesno[0] == 'y' || yesno[0] == 'Y'){
	  tmpl->channel[chan].pol[pol].nComp = 0;
	}
      }
    else if (key=='R')
      {
	char yesno[128];
	printf("Clear the template components for all channels with this Stokes parameter? (y/n) ");
	scanf("%s",yesno);
	if (yesno[0] == 'y' || yesno[0] == 'Y'){
	  for (i=0;i<tmpl->nchan;i++){
	    tmpl->channel[i].pol[pol].nComp = 0;
	  }
	}
      }
    else if (key=='t')
      showComp*=-1;
    else if (key=='u')
      zoom=0;
    else if (key=='f')
      {doFit(fx,fy,nbin,tmpl,chan,pol);saved=0;}
    else if (key=='F')
      {doFit_err(fx,fy,obs->chan[chan].pol[pol].devi,nbin,tmpl,chan,pol);saved=0;}
    else if (key=='g')
      {doFit_gerr(fx,fy,obs->chan[chan].pol[pol].devi,nbin,tmpl,chan,pol);saved=0;}
    else if (key=='G')
      {doFit_Gerr(fx,fy,obs->chan[chan].pol[pol].devi,nbin,tmpl,chan,pol);saved=0;}
    else if (key=='z') {
      float mx2,my2;
      cpgband(4,0,mx,my,&mx2,&my2,&key);
      if (mx2 < mx) minx = mx2;
      else minx = mx;
      if (mx2 > mx) maxx = mx2;
      else maxx = mx;
      zoom=1;
    }
    /*
    else if (key=='A') // Define a component
      {
	float mx2,my2,mx3,my3,mx4,my4;
	cpgband(4,0,mx,my,&mx3,&my3,&key);	
	if (key!='q'){
	  cpgband(5,0,mx,my,&mx4,&my4,&key);	
	  if (key!='q') {
	    if (tmpl->channel[chan].pol[pol].nComp == MAX_COMPONENTS-1)
	      {
		printf("Too many components\n");
	      }
	    else
	      {
		int nComp = tmpl->channel[chan].pol[pol].nComp;
		printf("Adding a component\n");
		tmpl->channel[chan].pol[pol].comp[nComp].concentration = 0.25/pow(fabs(mx3-mx),2);
		tmpl->channel[chan].pol[pol].comp[nComp].height = my4 ;
		tmpl->channel[chan].pol[pol].comp[nComp].centroid = mx ;
		(tmpl->channel[chan].pol[pol].nComp)++;
		saved=0;
	      }
	  }
	}
	key='a';
      }
      */
    else if (key=='s')
      {
	char fname[128];
	saved=1;
	printf("Please enter filename: ");
	scanf("%s",fname);
	saveTemplate(fname,tmpl);
      }
  } while (key!='q');
  if (saved==0){
    char yesno[128];
    printf("Do you wish to save the template? (y/n) ");
    scanf("%s",yesno);
    if (yesno[0] == 'y' || yesno[0] == 'Y')
      {
	char fname[128];
	printf("Please enter filename: ");
	scanf("%s",fname);
	saveTemplate(fname,tmpl);
      }
  }

  printf("Goodbye\n");
  cpgend();
}

void findMinMax(float *y,int n,float *min,float *max)
{
  int i;
  *min = *max = y[0];
  for (i=0;i<n;i++)
    {
      if (*min > y[i]) *min = y[i];
      if (*max < y[i]) *max = y[i];
    }
}

void readData(ptime_observation *obs,pheader *phead,fitsfile *fp)
{
  int status=0;
  int i,j,k,l;
  int initflag=0;
  int nval=0;
  int colnum;
  int nchan,nbin,npol;
  float ty[phead->nbin];
  float datScl[phead->nchan*phead->npol];
  float datOffs[phead->nchan*phead->npol];
  nchan = phead->nchan;
  nbin = phead->nbin;
  npol = phead->npol;

  //
  fits_movnam_hdu(fp,BINARY_TBL,(char *)"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,(char *)"DAT_SCL",&colnum,&status);
  if (status) {
    printf("Unable to find DAT_SCL in the subint table in FITS file\n");
    exit(1);
  }
  fits_read_col_flt(fp,colnum,1,1,phead->nchan*phead->npol,nval,datScl,&initflag,&status);


  fits_get_colnum(fp,CASEINSEN,(char *)"DAT_OFFS",&colnum,&status);
  if (status) {
    printf("Unable to find DAT_OFFS in the subint table in FITS file\n");
    exit(1);
  }
  fits_read_col_flt(fp,colnum,1,1,phead->nchan*phead->npol,nval,datOffs,&initflag,&status);

  fits_get_colnum(fp,CASEINSEN,(char *)"DATA",&colnum,&status);  
  if (status) {
    printf("Unable to find data in the subint table in FITS file\n");
    exit(1);
  }

  for (j=0;j<npol;j++)
    {
      for (i=0;i<nchan;i++)
	{
	  fits_read_col_flt(fp,colnum,1,j*(nchan*nbin)+i*nbin+1,nbin,nval,ty,&initflag,&status);
	  for (k=0;k<nbin;k++)
	    {
	      obs->chan[i].pol[j].val[k] = (ty[k]+datOffs[j*nchan+i])*datScl[j*nchan+i];
	    }
	}
    }
}

void allocateObsMemory(ptime_observation *obs,pheader *phead)
{
  int nbin,nchan,npol;
  int i,j;

  nchan = phead->nchan;
  npol  = phead->npol;
  nbin  = phead->nbin;

  obs->chan = (ptime_chan *)malloc(sizeof(ptime_chan)*nchan);
  obs->nchan = nchan;
  for (i=0;i<nchan;i++)
    {
      obs->chan[i].pol = (ptime_pol *)malloc(sizeof(ptime_pol)*npol);
      obs->chan[i].npol = npol;
      for (j=0;j<npol;j++)
	{
	  obs->chan[i].pol[j].val = (float *)malloc(sizeof(float)*nbin);
	  obs->chan[i].pol[j].nbin = nbin;
	}
    }
}

void deallocateMemory(ptime_observation *obs)
{
  int i,j;
  for (i=0;i<obs->nchan;i++)
    {
      for (j=0;j<obs->chan[i].npol;j++)
	{
	  
	  free(obs->chan[i].pol[j].val);
	}
      free(obs->chan[i].pol);
    }
  free(obs->chan);
  free(obs);
}


// Routines for the fitting
/* Collection of control (input) parameters. */
typedef struct {
    double ftol;      /* relative error desired in the sum of squares. */
    double xtol;      /* relative error between last two approximations. */
    double gtol;      /* orthogonality desired between fvec and its derivs. */
    double epsilon;   /* step used to calculate the jacobian. */
    double stepbound; /* initial bound to steps in the outer loop. */
    int maxcall;      /* maximum number of iterations. */
    int scale_diag;   /* UNDOCUMENTED, TESTWISE automatical diag rescaling? */
    int printflags;   /* OR'ed to produce more noise */
} lm_control_struct;

/* Collection of status (output) parameters. */
typedef struct {
    double fnorm;     /* norm of the residue vector fvec. */
    int nfev;	      /* actual number of iterations. */
    int info;	      /* status of minimization. */
} lm_status_struct;

/* Recommended control parameter settings. */
extern const lm_control_struct lm_control_double;
extern const lm_control_struct lm_control_float;


/* Standard monitoring routine. */
void lm_printout_std( int n_par, const double *par, int m_dat,
                      const void *data, const double *fvec,
                      int printflags, int iflag, int iter, int nfev);

/* Refined calculation of Eucledian norm, typically used in printout routine. */
double lm_enorm( int, const double * );

/* The actual minimization. */
void lmmin( int n_par, double *par, int m_dat, const void *data, 
            void (*evaluate) (const double *par, int m_dat, const void *data,
                              double *fvec, int *info),
            const lm_control_struct *control, lm_status_struct *status,
            void (*printout) (int n_par, const double *par, int m_dat,
                              const void *data, const double *fvec,
                              int printflags, int iflag, int iter, int nfev) );


/** Legacy low-level interface. **/

/* Alternative to lm_minimize, allowing full control, and read-out
   of auxiliary arrays. For usage, see implementation of lm_minimize. */
void lm_lmdif( int m, int n, double *x, double *fvec, double ftol,
	       double xtol, double gtol, int maxfev, double epsfcn,
	       double *diag, int mode, double factor, int *info, int *nfev,
	       double *fjac, int *ipvt, double *qtf, double *wa1,
	       double *wa2, double *wa3, double *wa4,
               void (*evaluate) (const double *par, int m_dat, const void *data,
                                 double *fvec, int *info),
               void (*printout) (int n_par, const double *par, int m_dat,
                                 const void *data, const double *fvec,
                                 int printflags, int iflag, int iter, int nfev),
               int printflags, const void *data );

extern const char *lm_infmsg[];
extern const char *lm_shortmsg[];

typedef struct {
    const double *t;
    const double *y;
    double (*f) (double t, const double *par);
} lmcurve_data_struct;

void lmcurve_fit( int n_par, double *par, int m_dat,
                  const double *t, const double *y,
                  double (*f)( double t, const double *par),
                  const lm_control_struct *control, lm_status_struct *status );
#define LM_MACHEP     DBL_EPSILON   /* resolution of arithmetic */
#define LM_DWARF      DBL_MIN       /* smallest nonzero number */
#define LM_SQRT_DWARF sqrt(DBL_MIN) /* square should not underflow */
#define LM_SQRT_GIANT sqrt(DBL_MAX) /* square should not overflow */
//#define LM_USERTOL    30*LM_MACHEP  /* users are recommended to require this */
#define LM_USERTOL 1e-9 // IMPORTANT TO SET TO SOMETHING SENSIBLE <<<<<

const lm_control_struct lm_control_double = {
    LM_USERTOL, LM_USERTOL, LM_USERTOL, LM_USERTOL, 100., 100, 1, 0 };
const lm_control_struct lm_control_float = {
    1.e-7, 1.e-7, 1.e-7, 1.e-7, 100., 100, 0, 0 };

void lm_lmpar( int n, double *r, int ldr, int *ipvt, double *diag,
	       double *qtb, double delta, double *par, double *x,
	       double *sdiag, double *aux, double *xdi );
void lm_qrfac( int m, int n, double *a, int pivot, int *ipvt,
	       double *rdiag, double *acnorm, double *wa );
void lm_qrsolv( int n, double *r, int ldr, int *ipvt, double *diag,
	        double *qtb, double *x, double *sdiag, double *wa );

#define MIN(a,b) (((a)<=(b)) ? (a) : (b))
#define MAX(a,b) (((a)>=(b)) ? (a) : (b))
#define SQR(x)   (x)*(x)



// Do the non-linear fit of the components
void doFit(float *fx,float *fy,int nbin,tmplStruct *tmpl,int chan,int pol)
{
  lm_status_struct status;
  lm_control_struct control = lm_control_double;
  control.printflags = 0; // monitor status (+1) and parameters (+2)
  double pval[tmpl->channel[chan].pol[pol].allVm*3];
  double datX[nbin],datY[nbin];
  int i,j,nFit;

  nFit = tmpl->channel[chan].pol[pol].allVm*3;

  // Fit for all parameters of all components
  int ivm=0;
  int icomp=0;
  for (i=0;i<tmpl->channel[chan].pol[pol].allVm;i++)
    {
      if (ivm != tmpl->channel[chan].pol[pol].comp[icomp].nVm-1)
      {
				pval[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height;
				pval[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration;
				pval[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid;
				ivm++;
      }
      else
      {
				pval[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height;
				pval[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration;
				pval[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid;
				icomp++;
				ivm=0;
      }
	  }
  // Fill up the fitting array
  for (i=0;i<nbin;i++)
    {
      datX[i] = (double)fx[i];
      datY[i] = (double)fy[i];
    }
  global_nComp = tmpl->channel[chan].pol[pol].allVm;
  lmcurve_fit(nFit,pval,nbin,datX,datY,nonlinearFunc,&control,&status);
  // Update components from the fit
  for (i=0;i<nFit/3;i++)
    printf("Fit results: %g %g %g\n",pval[3*i],pval[3*i+1],pval[3*i+2]);

  ivm=0;
  icomp=0;
  for (i=0;i<tmpl->channel[chan].pol[pol].allVm;i++)
    {
      if (ivm != tmpl->channel[chan].pol[pol].comp[icomp].nVm-1)
      {
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height = pval[3*i];
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = pval[3*i+1];
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = pval[3*i+2];
				ivm++;
      }
      else
      {
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height = pval[3*i];
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = pval[3*i+1];
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = pval[3*i+2];
				icomp++;
				ivm=0;
      }
    }
}

double nonlinearFunc( double x, const double *par)
{
  int k;
  double result=0.0;

  for (k=0;k<global_nComp;k++)
    result += fabs(par[3*k+0])*exp(par[3*k+1]*(cos((x-(par[3*k+2]))*2*M_PI)-1));
    //result += par[3*k+0]*exp(par[3*k+1]*(cos((x-(par[3*k+2]))*2*M_PI)-1));
  return result;
}

// -------------
// LM functions
// -------------
void lmcurve_evaluate( const double *par, int m_dat, const void *data,
                       double *fvec, int *info)
{
    int i;
    for ( i = 0; i < m_dat; i++ )
	fvec[i] =
            ((lmcurve_data_struct*)data)->y[i] -
            ((lmcurve_data_struct*)data)->f(
                ((lmcurve_data_struct*)data)->t[i], par);
    // *info = *info; /* to prevent a 'unused variable' warning */
}


void lmcurve_fit( int n_par, double *par, int m_dat, 
                  const double *t, const double *y,
                  double (*f)( double t, const double *par),
                  const lm_control_struct *control, lm_status_struct *status )
{
    lmcurve_data_struct data = { t, y, f };

    lmmin( n_par, par, m_dat, (const void*) &data,
           lmcurve_evaluate, control, status, lm_printout_std );
}



void lm_printout_std( int n_par, const double *par, int m_dat,
                      const void *data, const double *fvec,
                      int printflags, int iflag, int iter, int nfev)
/*
 *       data  : for soft control of printout behaviour, add control
 *                 variables to the data struct
 *       iflag : 0 (init) 1 (outer loop) 2(inner loop) -1(terminated)
 *       iter  : outer loop counter
 *       nfev  : number of calls to *evaluate
 */
{
    if( !printflags )
        return;

    int i;

    if( printflags & 1 ){
        /* location of printout call within lmdif */
        if (iflag == 2) {
            printf("trying step in gradient direction  ");
        } else if (iflag == 1) {
            printf("determining gradient (iteration %2d)", iter);
        } else if (iflag == 0) {
            printf("starting minimization              ");
        } else if (iflag == -1) {
            printf("terminated after %3d evaluations   ", nfev);
        }
    }

    if( printflags & 2 ){
        printf("  par: ");
        for (i = 0; i < n_par; ++i)
            printf(" %18.11g", par[i]);
        printf(" => norm: %18.11g", lm_enorm(m_dat, fvec));
    }

    if( printflags & 3 )
        printf( "\n" );

    if ( (printflags & 8) || ((printflags & 4) && iflag == -1) ) {
	printf("  residuals:\n");
	for (i = 0; i < m_dat; ++i)
	    printf("    fvec[%2d]=%12g\n", i, fvec[i] );
    }
}


/*****************************************************************************/
/*  lm_minimize (intermediate-level interface)                               */
/*****************************************************************************/

void lmmin( int n_par, double *par, int m_dat, const void *data,
            void (*evaluate) (const double *par, int m_dat, const void *data,
                              double *fvec, int *info),
            const lm_control_struct *control, lm_status_struct *status,
            void (*printout) (int n_par, const double *par, int m_dat,
                              const void *data, const double *fvec,
                              int printflags, int iflag, int iter, int nfev) )
{

/*** allocate work space. ***/

    double *fvec, *diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
    int *ipvt;

    int n = n_par;
    int m = m_dat;

    if ( (fvec = (double *) malloc(m * sizeof(double))) == NULL ||
	 (diag = (double *) malloc(n * sizeof(double))) == NULL ||
	 (qtf  = (double *) malloc(n * sizeof(double))) == NULL ||
	 (fjac = (double *) malloc(n*m*sizeof(double))) == NULL ||
	 (wa1  = (double *) malloc(n * sizeof(double))) == NULL ||
	 (wa2  = (double *) malloc(n * sizeof(double))) == NULL ||
	 (wa3  = (double *) malloc(n * sizeof(double))) == NULL ||
	 (wa4  = (double *) malloc(m * sizeof(double))) == NULL ||
	 (ipvt = (int *)    malloc(n * sizeof(int)   )) == NULL    ) {
	status->info = 9;
	return;
    }

    int j;
    if( ! control->scale_diag )
        for( j=0; j<n_par; ++j )
            diag[j] = 1;

/*** perform fit. ***/

    status->info = 0;

    /* this goes through the modified legacy interface: */
    lm_lmdif( m, n, par, fvec, control->ftol, control->xtol, control->gtol,
              control->maxcall * (n + 1), control->epsilon, diag,
              ( control->scale_diag ? 1 : 2 ),
              control->stepbound, &(status->info),
              &(status->nfev), fjac, ipvt, qtf, wa1, wa2, wa3, wa4,
              evaluate, printout, control->printflags, data );

    if ( printout )
        (*printout)( n, par, m, data, fvec,
                     control->printflags, -1, 0, status->nfev );
    status->fnorm = lm_enorm(m, fvec);
    if ( status->info < 0 )
	status->info = 11;

/*** clean up. ***/

    free(fvec);
    free(diag);
    free(qtf);
    free(fjac);
    free(wa1);
    free(wa2);
    free(wa3);
    free(wa4);
    free(ipvt);
} /*** lm_minimize. ***/

/*****************************************************************************/
/*  lm_enorm (Euclidean norm)                                                */
/*****************************************************************************/

double lm_enorm(int n, const double *x)
{
/*     Given an n-vector x, this function calculates the
 *     euclidean norm of x.
 *
 *     The euclidean norm is computed by accumulating the sum of
 *     squares in three different sums. The sums of squares for the
 *     small and large components are scaled so that no overflows
 *     occur. Non-destructive underflows are permitted. Underflows
 *     and overflows do not occur in the computation of the unscaled
 *     sum of squares for the intermediate components.
 *     The definitions of small, intermediate and large components
 *     depend on two constants, LM_SQRT_DWARF and LM_SQRT_GIANT. The main
 *     restrictions on these constants are that LM_SQRT_DWARF**2 not
 *     underflow and LM_SQRT_GIANT**2 not overflow.
 *
 *     Parameters
 *
 *	n is a positive integer input variable.
 *
 *	x is an input array of length n.
 */
    int i;
    double agiant, s1, s2, s3, xabs, x1max, x3max, temp;

    s1 = 0;
    s2 = 0;
    s3 = 0;
    x1max = 0;
    x3max = 0;
    agiant = LM_SQRT_GIANT / n;

    /** sum squares. **/

    for (i = 0; i < n; i++) {
	xabs = fabs(x[i]);
	if (xabs > LM_SQRT_DWARF) {
            if ( xabs < agiant ) {
                s2 += xabs * xabs;
            } else if ( xabs > x1max ) {
		temp = x1max / xabs;
		s1 = 1 + s1 * SQR(temp);
		x1max = xabs;
	    } else {
		temp = xabs / x1max;
		s1 += SQR(temp);
	    }
	} else if ( xabs > x3max ) {
	    temp = x3max / xabs;
	    s3 = 1 + s3 * SQR(temp);
	    x3max = xabs;
	} else if (xabs != 0.) {
            temp = xabs / x3max;
            s3 += SQR(temp);
	}
    }

    /** calculation of norm. **/

    if (s1 != 0)
	return x1max * sqrt(s1 + (s2 / x1max) / x1max);
    else if (s2 != 0)
        if (s2 >= x3max)
            return sqrt(s2 * (1 + (x3max / s2) * (x3max * s3)));
        else
            return sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
    else
        return x3max * sqrt(s3);

} /*** lm_enorm. ***/

void lm_lmdif( int m, int n, double *x, double *fvec, double ftol,
	       double xtol, double gtol, int maxfev, double epsfcn,
	       double *diag, int mode, double factor, int *info, int *nfev,
	       double *fjac, int *ipvt, double *qtf, double *wa1,
	       double *wa2, double *wa3, double *wa4,
               void (*evaluate) (const double *par, int m_dat, const void *data,
                                 double *fvec, int *info),
               void (*printout) (int n_par, const double *par, int m_dat,
                                 const void *data, const double *fvec,
                                 int printflags, int iflag, int iter, int nfev),
	       int printflags, const void *data )
{
/*
 *   The purpose of lmdif is to minimize the sum of the squares of
 *   m nonlinear functions in n variables by a modification of
 *   the levenberg-marquardt algorithm. The user must provide a
 *   subroutine evaluate which calculates the functions. The jacobian
 *   is then calculated by a forward-difference approximation.
 *
 *   The multi-parameter interface lm_lmdif is for users who want
 *   full control and flexibility. Most users will be better off using
 *   the simpler interface lmmin provided above.
 *
 *   Parameters:
 *
 *	m is a positive integer input variable set to the number
 *	  of functions.
 *
 *	n is a positive integer input variable set to the number
 *	  of variables; n must not exceed m.
 *
 *	x is an array of length n. On input x must contain an initial
 *        estimate of the solution vector. On OUTPUT x contains the
 *        final estimate of the solution vector.
 *
 *	fvec is an OUTPUT array of length m which contains
 *	  the functions evaluated at the output x.
 *
 *	ftol is a nonnegative input variable. Termination occurs when
 *        both the actual and predicted relative reductions in the sum
 *        of squares are at most ftol. Therefore, ftol measures the
 *        relative error desired in the sum of squares.
 *
 *	xtol is a nonnegative input variable. Termination occurs when
 *        the relative error between two consecutive iterates is at
 *        most xtol. Therefore, xtol measures the relative error desired
 *        in the approximate solution.
 *
 *	gtol is a nonnegative input variable. Termination occurs when
 *        the cosine of the angle between fvec and any column of the
 *        jacobian is at most gtol in absolute value. Therefore, gtol
 *        measures the orthogonality desired between the function vector
 *        and the columns of the jacobian.
 *
 *	maxfev is a positive integer input variable. Termination
 *	  occurs when the number of calls to lm_fcn is at least
 *	  maxfev by the end of an iteration.
 *
 *	epsfcn is an input variable used in choosing a step length for
 *        the forward-difference approximation. The relative errors in
 *        the functions are assumed to be of the order of epsfcn.
 *
 *	diag is an array of length n. If mode = 1 (see below), diag is
 *        internally set. If mode = 2, diag must contain positive entries
 *        that serve as multiplicative scale factors for the variables.
 *
 *	mode is an integer input variable. If mode = 1, the
 *	  variables will be scaled internally. If mode = 2,
 *	  the scaling is specified by the input diag.
 *
 *	factor is a positive input variable used in determining the
 *	  initial step bound. This bound is set to the product of
 *	  factor and the euclidean norm of diag*x if nonzero, or else
 *	  to factor itself. In most cases factor should lie in the
 *	  interval (0.1,100.0). Generally, the value 100.0 is recommended.
 *
 *	info is an integer OUTPUT variable that indicates the termination
 *        status of lm_lmdif as follows:
 *
 *        info < 0  termination requested by user-supplied routine *evaluate;
 *
 *	  info = 0  fnorm almost vanishing;
 *
 *	  info = 1  both actual and predicted relative reductions
 *		    in the sum of squares are at most ftol;
 *
 *	  info = 2  relative error between two consecutive iterates
 *		    is at most xtol;
 *
 *	  info = 3  conditions for info = 1 and info = 2 both hold;
 *
 *	  info = 4  the cosine of the angle between fvec and any
 *		    column of the jacobian is at most gtol in
 *		    absolute value;
 *
 *	  info = 5  number of calls to lm_fcn has reached or
 *		    exceeded maxfev;
 *
 *	  info = 6  ftol is too small: no further reduction in
 *		    the sum of squares is possible;
 *
 *	  info = 7  xtol is too small: no further improvement in
 *		    the approximate solution x is possible;
 *
 *	  info = 8  gtol is too small: fvec is orthogonal to the
 *		    columns of the jacobian to machine precision;
 *
 *	  info =10  improper input parameters;
 *
 *	nfev is an OUTPUT variable set to the number of calls to the
 *        user-supplied routine *evaluate.
 *
 *	fjac is an OUTPUT m by n array. The upper n by n submatrix
 *	  of fjac contains an upper triangular matrix r with
 *	  diagonal elements of nonincreasing magnitude such that
 *
 *		pT*(jacT*jac)*p = rT*r
 *
 *              (NOTE: T stands for matrix transposition),
 *
 *	  where p is a permutation matrix and jac is the final
 *	  calculated jacobian. Column j of p is column ipvt(j)
 *	  (see below) of the identity matrix. The lower trapezoidal
 *	  part of fjac contains information generated during
 *	  the computation of r.
 *
 *	ipvt is an integer OUTPUT array of length n. It defines a
 *        permutation matrix p such that jac*p = q*r, where jac is
 *        the final calculated jacobian, q is orthogonal (not stored),
 *        and r is upper triangular with diagonal elements of
 *        nonincreasing magnitude. Column j of p is column ipvt(j)
 *        of the identity matrix.
 *
 *	qtf is an OUTPUT array of length n which contains
 *	  the first n elements of the vector (q transpose)*fvec.
 *
 *	wa1, wa2, and wa3 are work arrays of length n.
 *
 *	wa4 is a work array of length m, used among others to hold
 *        residuals from evaluate.
 *
 *      evaluate points to the subroutine which calculates the
 *        m nonlinear functions. Implementations should be written as follows:
 *
 *        void evaluate( double* par, int m_dat, void *data,
 *                       double* fvec, int *info )
 *        {
 *           // for ( i=0; i<m_dat; ++i )
 *           //     calculate fvec[i] for given parameters par;
 *           // to stop the minimization, 
 *           //     set *info to a negative integer.
 *        }
 *
 *      printout points to the subroutine which informs about fit progress.
 *        Call with printout=0 if no printout is desired.
 *        Call with printout=lm_printout_std to use the default implementation.
 *
 *      printflags is passed to printout.
 *
 *      data is an input pointer to an arbitrary structure that is passed to
 *        evaluate. Typically, it contains experimental data to be fitted.
 *
 */
    int i, iter, j;
    double actred, delta, dirder, eps, fnorm, fnorm1, gnorm, par, pnorm,
	prered, ratio, step, sum, temp, temp1, temp2, temp3, xnorm;
    static double p1 = 0.1;
    static double p0001 = 1.0e-4;

    *nfev = 0;			/* function evaluation counter */
    iter = 0;			/* outer loop counter */
    par = 0;			/* levenberg-marquardt parameter */
    delta = 0;	 /* to prevent a warning (initialization within if-clause) */
    xnorm = 0;	 /* ditto */
    temp = MAX(epsfcn, LM_MACHEP);
    eps = sqrt(temp); /* for calculating the Jacobian by forward differences */

/*** lmdif: check input parameters for errors. ***/

    if ((n <= 0) || (m < n) || (ftol < 0.)
	|| (xtol < 0.) || (gtol < 0.) || (maxfev <= 0) || (factor <= 0.)) {
	*info = 10;		// invalid parameter
	return;
    }
    if (mode == 2) {		/* scaling by diag[] */
	for (j = 0; j < n; j++) {	/* check for nonpositive elements */
	    if (diag[j] <= 0.0) {
		*info = 10;	// invalid parameter
		return;
	    }
	}
    }
#ifdef LMFIT_DEBUG_MESSAGES
    printf("lmdif\n");
#endif

/*** lmdif: evaluate function at starting point and calculate norm. ***/

    *info = 0;
    (*evaluate) (x, m, data, fvec, info);
    ++(*nfev);
    if( printout )
        (*printout) (n, x, m, data, fvec, printflags, 0, 0, *nfev);
    if (*info < 0)
	return;
    fnorm = lm_enorm(m, fvec);
    if( fnorm <= LM_DWARF ){
        *info = 0;
        return;
    }

/*** lmdif: the outer loop. ***/

    do {
#ifdef LMFIT_DEBUG_MESSAGES
	printf("lmdif/ outer loop iter=%d nfev=%d fnorm=%.10e\n",
	       iter, *nfev, fnorm);
#endif

/*** outer: calculate the Jacobian. ***/

	for (j = 0; j < n; j++) {
	    temp = x[j];
	    step = MAX(eps*eps, eps * fabs(temp));
	    x[j] = temp + step; /* replace temporarily */
	    *info = 0;
	    (*evaluate) (x, m, data, wa4, info);
            ++(*nfev);
            if( printout )
                (*printout) (n, x, m, data, wa4, printflags, 1, iter, *nfev);
	    if (*info < 0)
		return;	/* user requested break */
	    for (i = 0; i < m; i++)
		fjac[j*m+i] = (wa4[i] - fvec[i]) / step;
	    x[j] = temp; /* restore */
	}
#ifdef LMFIT_DEBUG_MATRIX
	/* print the entire matrix */
	for (i = 0; i < m; i++) {
	    for (j = 0; j < n; j++)
		printf("%.5e ", fjac[j*m+i]);
	    printf("\n");
	}
#endif

/*** outer: compute the qr factorization of the Jacobian. ***/

	lm_qrfac(m, n, fjac, 1, ipvt, wa1, wa2, wa3);
        /* return values are ipvt, wa1=rdiag, wa2=acnorm */

	if (!iter) { 
            /* first iteration only */
	    if (mode != 2) {
                /* diag := norms of the columns of the initial Jacobian */
		for (j = 0; j < n; j++) {
		    diag[j] = wa2[j];
		    if (wa2[j] == 0.)
			diag[j] = 1.;
		}
	    }
            /* use diag to scale x, then calculate the norm */
	    for (j = 0; j < n; j++)
		wa3[j] = diag[j] * x[j];
	    xnorm = lm_enorm(n, wa3);
            /* initialize the step bound delta. */
	    delta = factor * xnorm;
	    if (delta == 0.)
		delta = factor;
	} else {
            if (mode != 2) {
                for (j = 0; j < n; j++)
                    diag[j] = MAX( diag[j], wa2[j] );
            }
        }

/*** outer: form (q transpose)*fvec and store first n components in qtf. ***/

	for (i = 0; i < m; i++)
	    wa4[i] = fvec[i];

	for (j = 0; j < n; j++) {
	    temp3 = fjac[j*m+j];
	    if (temp3 != 0.) {
		sum = 0;
		for (i = j; i < m; i++)
		    sum += fjac[j*m+i] * wa4[i];
		temp = -sum / temp3;
		for (i = j; i < m; i++)
		    wa4[i] += fjac[j*m+i] * temp;
	    }
	    fjac[j*m+j] = wa1[j];
	    qtf[j] = wa4[j];
	}

/*** outer: compute norm of scaled gradient and test for convergence. ***/

	gnorm = 0;
        for (j = 0; j < n; j++) {
            if (wa2[ipvt[j]] == 0)
                continue;
            sum = 0.;
            for (i = 0; i <= j; i++)
                sum += fjac[j*m+i] * qtf[i];
            gnorm = MAX( gnorm, fabs( sum / wa2[ipvt[j]] / fnorm ) );
        }

	if (gnorm <= gtol) {
	    *info = 4;
	    return;
	}

/*** the inner loop. ***/
	do {
#ifdef LMFIT_DEBUG_MESSAGES
	    printf("lmdif/ inner loop iter=%d nfev=%d\n", iter, *nfev);
#endif

/*** inner: determine the levenberg-marquardt parameter. ***/

	    lm_lmpar( n, fjac, m, ipvt, diag, qtf, delta, &par,
                      wa1, wa2, wa4, wa3 );
            /* used return values are fjac (partly), par, wa1=x, wa3=diag*x */

	    for (j = 0; j < n; j++)
		wa2[j] = x[j] - wa1[j]; /* new parameter vector ? */

	    pnorm = lm_enorm(n, wa3);

            /* at first call, adjust the initial step bound. */

	    if (*nfev <= 1+n)
		delta = MIN(delta, pnorm);

/*** inner: evaluate the function at x + p and calculate its norm. ***/

	    *info = 0;
	    (*evaluate) (wa2, m, data, wa4, info);
            ++(*nfev);
            if( printout )
                (*printout) (n, wa2, m, data, wa4, printflags, 2, iter, *nfev);
	    if (*info < 0)
		return; /* user requested break. */

	    fnorm1 = lm_enorm(m, wa4);
#ifdef LMFIT_DEBUG_MESSAGES
	    printf("lmdif/ pnorm %.10e  fnorm1 %.10e  fnorm %.10e"
		   " delta=%.10e par=%.10e\n",
		   pnorm, fnorm1, fnorm, delta, par);
#endif

/*** inner: compute the scaled actual reduction. ***/

	    if (p1 * fnorm1 < fnorm)
		actred = 1 - SQR(fnorm1 / fnorm);
	    else
		actred = -1;

/*** inner: compute the scaled predicted reduction and 
     the scaled directional derivative. ***/

	    for (j = 0; j < n; j++) {
		wa3[j] = 0;
		for (i = 0; i <= j; i++)
		    wa3[i] -= fjac[j*m+i] * wa1[ipvt[j]];
	    }
	    temp1 = lm_enorm(n, wa3) / fnorm;
	    temp2 = sqrt(par) * pnorm / fnorm;
	    prered = SQR(temp1) + 2 * SQR(temp2);
	    dirder = -(SQR(temp1) + SQR(temp2));

/*** inner: compute the ratio of the actual to the predicted reduction. ***/

	    ratio = prered != 0 ? actred / prered : 0;
#ifdef LMFIT_DEBUG_MESSAGES
	    printf("lmdif/ actred=%.10e prered=%.10e ratio=%.10e"
		   " sq(1)=%.10e sq(2)=%.10e dd=%.10e\n",
		   actred, prered, prered != 0 ? ratio : 0.,
		   SQR(temp1), SQR(temp2), dirder);
#endif

/*** inner: update the step bound. ***/

	    if (ratio <= 0.25) {
		if (actred >= 0.)
		    temp = 0.5;
		else
		    temp = 0.5 * dirder / (dirder + 0.55 * actred);
		if (p1 * fnorm1 >= fnorm || temp < p1)
		    temp = p1;
		delta = temp * MIN(delta, pnorm / p1);
		par /= temp;
	    } else if (par == 0. || ratio >= 0.75) {
		delta = pnorm / 0.5;
		par *= 0.5;
	    }

/*** inner: test for successful iteration. ***/

	    if (ratio >= p0001) {
                /* yes, success: update x, fvec, and their norms. */
		for (j = 0; j < n; j++) {
		    x[j] = wa2[j];
		    wa2[j] = diag[j] * x[j];
		}
		for (i = 0; i < m; i++)
		    fvec[i] = wa4[i];
		xnorm = lm_enorm(n, wa2);
		fnorm = fnorm1;
		iter++;
	    }
#ifdef LMFIT_DEBUG_MESSAGES
	    else {
		printf("ATTN: iteration considered unsuccessful\n");
	    }
#endif

/*** inner: test for convergence. ***/

            if( fnorm<=LM_DWARF ){
                *info = 0;
                return;
            }

	    *info = 0;
	    if (fabs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1)
		*info = 1;
	    if (delta <= xtol * xnorm)
		*info += 2;
	    if (*info != 0)
		return;

/*** inner: tests for termination and stringent tolerances. ***/

	    if (*nfev >= maxfev){
		*info = 5;
                return;
            }
	    if (fabs(actred) <= LM_MACHEP &&
		prered <= LM_MACHEP && 0.5 * ratio <= 1){
		*info = 6;
                return;
            }
	    if (delta <= LM_MACHEP * xnorm){
		*info = 7;
                return;
            }
	    if (gnorm <= LM_MACHEP){
		*info = 8;
		return;
            }

/*** inner: end of the loop. repeat if iteration unsuccessful. ***/

	} while (ratio < p0001);

/*** outer: end of the loop. ***/

    } while (1);

} /*** lm_lmdif. ***/


/*****************************************************************************/
/*  lm_qrfac (QR factorisation, from lapack)                                 */
/*****************************************************************************/

void lm_qrfac(int m, int n, double *a, int pivot, int *ipvt,
	      double *rdiag, double *acnorm, double *wa)
{
/*
 *     This subroutine uses householder transformations with column
 *     pivoting (optional) to compute a qr factorization of the
 *     m by n matrix a. That is, qrfac determines an orthogonal
 *     matrix q, a permutation matrix p, and an upper trapezoidal
 *     matrix r with diagonal elements of nonincreasing magnitude,
 *     such that a*p = q*r. The householder transformation for
 *     column k, k = 1,2,...,min(m,n), is of the form
 *
 *	    i - (1/u(k))*u*uT
 *
 *     where u has zeroes in the first k-1 positions. The form of
 *     this transformation and the method of pivoting first
 *     appeared in the corresponding linpack subroutine.
 *
 *     Parameters:
 *
 *	m is a positive integer input variable set to the number
 *	  of rows of a.
 *
 *	n is a positive integer input variable set to the number
 *	  of columns of a.
 *
 *	a is an m by n array. On input a contains the matrix for
 *	  which the qr factorization is to be computed. On OUTPUT
 *	  the strict upper trapezoidal part of a contains the strict
 *	  upper trapezoidal part of r, and the lower trapezoidal
 *	  part of a contains a factored form of q (the non-trivial
 *	  elements of the u vectors described above).
 *
 *	pivot is a logical input variable. If pivot is set true,
 *	  then column pivoting is enforced. If pivot is set false,
 *	  then no column pivoting is done.
 *
 *	ipvt is an integer OUTPUT array of length lipvt. This array
 *	  defines the permutation matrix p such that a*p = q*r.
 *	  Column j of p is column ipvt(j) of the identity matrix.
 *	  If pivot is false, ipvt is not referenced.
 *
 *	rdiag is an OUTPUT array of length n which contains the
 *	  diagonal elements of r.
 *
 *	acnorm is an OUTPUT array of length n which contains the
 *	  norms of the corresponding columns of the input matrix a.
 *	  If this information is not needed, then acnorm can coincide
 *	  with rdiag.
 *
 *	wa is a work array of length n. If pivot is false, then wa
 *	  can coincide with rdiag.
 *
 */
    int i, j, k, kmax, minmn;
    double ajnorm, sum, temp;

/*** qrfac: compute initial column norms and initialize several arrays. ***/

    for (j = 0; j < n; j++) {
	acnorm[j] = lm_enorm(m, &a[j*m]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	if (pivot)
	    ipvt[j] = j;
    }
#ifdef LMFIT_DEBUG_MESSAGES
    printf("qrfac\n");
#endif

/*** qrfac: reduce a to r with householder transformations. ***/

    minmn = MIN(m, n);
    for (j = 0; j < minmn; j++) {
	if (!pivot)
	    goto pivot_ok;

        /** bring the column of largest norm into the pivot position. **/

	kmax = j;
	for (k = j + 1; k < n; k++)
	    if (rdiag[k] > rdiag[kmax])
		kmax = k;
	if (kmax == j)
	    goto pivot_ok;

	for (i = 0; i < m; i++) {
	    temp = a[j*m+i];
	    a[j*m+i] = a[kmax*m+i];
	    a[kmax*m+i] = temp;
	}
	rdiag[kmax] = rdiag[j];
	wa[kmax] = wa[j];
	k = ipvt[j];
	ipvt[j] = ipvt[kmax];
	ipvt[kmax] = k;

      pivot_ok:
        /** compute the Householder transformation to reduce the
            j-th column of a to a multiple of the j-th unit vector. **/

	ajnorm = lm_enorm(m-j, &a[j*m+j]);
	if (ajnorm == 0.) {
	    rdiag[j] = 0;
	    continue;
	}

	if (a[j*m+j] < 0.)
	    ajnorm = -ajnorm;
	for (i = j; i < m; i++)
	    a[j*m+i] /= ajnorm;
	a[j*m+j] += 1;

        /** apply the transformation to the remaining columns
            and update the norms. **/

	for (k = j + 1; k < n; k++) {
	    sum = 0;

	    for (i = j; i < m; i++)
		sum += a[j*m+i] * a[k*m+i];

	    temp = sum / a[j + m * j];

	    for (i = j; i < m; i++)
		a[k*m+i] -= temp * a[j*m+i];

	    if (pivot && rdiag[k] != 0.) {
		temp = a[m * k + j] / rdiag[k];
		temp = MAX(0., 1 - temp * temp);
		rdiag[k] *= sqrt(temp);
		temp = rdiag[k] / wa[k];
		if ( 0.05 * SQR(temp) <= LM_MACHEP ) {
		    rdiag[k] = lm_enorm(m-j-1, &a[m*k+j+1]);
		    wa[k] = rdiag[k];
		}
	    }
	}

	rdiag[j] = -ajnorm;
    }
}

/*****************************************************************************/
/*  lm_lmpar (determine Levenberg-Marquardt parameter)                       */
/*****************************************************************************/

void lm_lmpar(int n, double *r, int ldr, int *ipvt, double *diag,
	      double *qtb, double delta, double *par, double *x,
	      double *sdiag, double *aux, double *xdi)
{
/*     Given an m by n matrix a, an n by n nonsingular diagonal
 *     matrix d, an m-vector b, and a positive number delta,
 *     the problem is to determine a value for the parameter
 *     par such that if x solves the system
 *
 *	    a*x = b  and  sqrt(par)*d*x = 0
 *
 *     in the least squares sense, and dxnorm is the euclidean
 *     norm of d*x, then either par=0 and (dxnorm-delta) < 0.1*delta,
 *     or par>0 and abs(dxnorm-delta) < 0.1*delta.
 *
 *     Using lm_qrsolv, this subroutine completes the solution of the problem
 *     if it is provided with the necessary information from the
 *     qr factorization, with column pivoting, of a. That is, if
 *     a*p = q*r, where p is a permutation matrix, q has orthogonal
 *     columns, and r is an upper triangular matrix with diagonal
 *     elements of nonincreasing magnitude, then lmpar expects
 *     the full upper triangle of r, the permutation matrix p,
 *     and the first n components of qT*b. On output
 *     lmpar also provides an upper triangular matrix s such that
 *
 *	    pT*(aT*a + par*d*d)*p = sT*s.
 *
 *     s is employed within lmpar and may be of separate interest.
 *
 *     Only a few iterations are generally needed for convergence
 *     of the algorithm. If, however, the limit of 10 iterations
 *     is reached, then the output par will contain the best
 *     value obtained so far.
 *
 *     parameters:
 *
 *	n is a positive integer input variable set to the order of r.
 *
 *	r is an n by n array. on input the full upper triangle
 *	  must contain the full upper triangle of the matrix r.
 *	  on OUTPUT the full upper triangle is unaltered, and the
 *	  strict lower triangle contains the strict upper triangle
 *	  (transposed) of the upper triangular matrix s.
 *
 *	ldr is a positive integer input variable not less than n
 *	  which specifies the leading dimension of the array r.
 *
 *	ipvt is an integer input array of length n which defines the
 *	  permutation matrix p such that a*p = q*r. column j of p
 *	  is column ipvt(j) of the identity matrix.
 *
 *	diag is an input array of length n which must contain the
 *	  diagonal elements of the matrix d.
 *
 *	qtb is an input array of length n which must contain the first
 *	  n elements of the vector (q transpose)*b.
 *
 *	delta is a positive input variable which specifies an upper
 *	  bound on the euclidean norm of d*x.
 *
 *	par is a nonnegative variable. on input par contains an
 *	  initial estimate of the levenberg-marquardt parameter.
 *	  on OUTPUT par contains the final estimate.
 *
 *	x is an OUTPUT array of length n which contains the least
 *	  squares solution of the system a*x = b, sqrt(par)*d*x = 0,
 *	  for the output par.
 *
 *	sdiag is an array of length n which contains the
 *	  diagonal elements of the upper triangular matrix s.
 *
 *	aux is a multi-purpose work array of length n.
 *
 *	xdi is a work array of length n. On OUTPUT: diag[j] * x[j].
 *
 */
    int i, iter, j, nsing;
    double dxnorm, fp, fp_old, gnorm, parc, parl, paru;
    double sum, temp;
    static double p1 = 0.1;

#ifdef LMFIT_DEBUG_MESSAGES
    printf("lmpar\n");
#endif

/*** lmpar: compute and store in x the gauss-newton direction. if the
     jacobian is rank-deficient, obtain a least squares solution. ***/

    nsing = n;
    for (j = 0; j < n; j++) {
	aux[j] = qtb[j];
	if (r[j * ldr + j] == 0 && nsing == n)
	    nsing = j;
	if (nsing < n)
	    aux[j] = 0;
    }
#ifdef LMFIT_DEBUG_MESSAGES
    printf("nsing %d ", nsing);
#endif
    for (j = nsing - 1; j >= 0; j--) {
	aux[j] = aux[j] / r[j + ldr * j];
	temp = aux[j];
	for (i = 0; i < j; i++)
	    aux[i] -= r[j * ldr + i] * temp;
    }

    for (j = 0; j < n; j++)
	x[ipvt[j]] = aux[j];

/*** lmpar: initialize the iteration counter, evaluate the function at the
     origin, and test for acceptance of the gauss-newton direction. ***/

    iter = 0;
    for (j = 0; j < n; j++)
	xdi[j] = diag[j] * x[j];
    dxnorm = lm_enorm(n, xdi);
    fp = dxnorm - delta;
    if (fp <= p1 * delta) {
#ifdef LMFIT_DEBUG_MESSAGES
	printf("lmpar/ terminate (fp<p1*delta)\n");
#endif
	*par = 0;
	return;
    }

/*** lmpar: if the jacobian is not rank deficient, the newton
     step provides a lower bound, parl, for the 0. of
     the function. otherwise set this bound to 0.. ***/

    parl = 0;
    if (nsing >= n) {
	for (j = 0; j < n; j++)
	    aux[j] = diag[ipvt[j]] * xdi[ipvt[j]] / dxnorm;

	for (j = 0; j < n; j++) {
	    sum = 0.;
	    for (i = 0; i < j; i++)
		sum += r[j * ldr + i] * aux[i];
	    aux[j] = (aux[j] - sum) / r[j + ldr * j];
	}
	temp = lm_enorm(n, aux);
	parl = fp / delta / temp / temp;
    }

/*** lmpar: calculate an upper bound, paru, for the 0. of the function. ***/

    for (j = 0; j < n; j++) {
	sum = 0;
	for (i = 0; i <= j; i++)
	    sum += r[j * ldr + i] * qtb[i];
	aux[j] = sum / diag[ipvt[j]];
    }
    gnorm = lm_enorm(n, aux);
    paru = gnorm / delta;
    if (paru == 0.)
	paru = LM_DWARF / MIN(delta, p1);

/*** lmpar: if the input par lies outside of the interval (parl,paru),
     set par to the closer endpoint. ***/

    *par = MAX(*par, parl);
    *par = MIN(*par, paru);
    if (*par == 0.)
	*par = gnorm / dxnorm;
#ifdef LMFIT_DEBUG_MESSAGES
    printf("lmpar/ parl %.4e  par %.4e  paru %.4e\n", parl, *par, paru);
#endif

/*** lmpar: iterate. ***/

    for (;; iter++) {

        /** evaluate the function at the current value of par. **/

	if (*par == 0.)
	    *par = MAX(LM_DWARF, 0.001 * paru);
	temp = sqrt(*par);
	for (j = 0; j < n; j++)
	    aux[j] = temp * diag[j];

	lm_qrsolv( n, r, ldr, ipvt, aux, qtb, x, sdiag, xdi );
        /* return values are r, x, sdiag */

	for (j = 0; j < n; j++)
	    xdi[j] = diag[j] * x[j]; /* used as output */
	dxnorm = lm_enorm(n, xdi);
	fp_old = fp;
	fp = dxnorm - delta;
        
        /** if the function is small enough, accept the current value
            of par. Also test for the exceptional cases where parl
            is zero or the number of iterations has reached 10. **/

	if (fabs(fp) <= p1 * delta
	    || (parl == 0. && fp <= fp_old && fp_old < 0.)
	    || iter == 10)
	    break; /* the only exit from the iteration. */
        
        /** compute the Newton correction. **/

	for (j = 0; j < n; j++)
	    aux[j] = diag[ipvt[j]] * xdi[ipvt[j]] / dxnorm;

	for (j = 0; j < n; j++) {
	    aux[j] = aux[j] / sdiag[j];
	    for (i = j + 1; i < n; i++)
		aux[i] -= r[j * ldr + i] * aux[j];
	}
	temp = lm_enorm(n, aux);
	parc = fp / delta / temp / temp;

        /** depending on the sign of the function, update parl or paru. **/

	if (fp > 0)
	    parl = MAX(parl, *par);
	else if (fp < 0)
	    paru = MIN(paru, *par);
	/* the case fp==0 is precluded by the break condition  */
        
        /** compute an improved estimate for par. **/
        
	*par = MAX(parl, *par + parc);
        
    }

} /*** lm_lmpar. ***/

/*****************************************************************************/
/*  lm_qrsolv (linear least-squares)                                         */
/*****************************************************************************/

void lm_qrsolv(int n, double *r, int ldr, int *ipvt, double *diag,
	       double *qtb, double *x, double *sdiag, double *wa)
{
/*
 *     Given an m by n matrix a, an n by n diagonal matrix d,
 *     and an m-vector b, the problem is to determine an x which
 *     solves the system
 *
 *	    a*x = b  and  d*x = 0
 *
 *     in the least squares sense.
 *
 *     This subroutine completes the solution of the problem
 *     if it is provided with the necessary information from the
 *     qr factorization, with column pivoting, of a. That is, if
 *     a*p = q*r, where p is a permutation matrix, q has orthogonal
 *     columns, and r is an upper triangular matrix with diagonal
 *     elements of nonincreasing magnitude, then qrsolv expects
 *     the full upper triangle of r, the permutation matrix p,
 *     and the first n components of (q transpose)*b. The system
 *     a*x = b, d*x = 0, is then equivalent to
 *
 *	    r*z = qT*b,  pT*d*p*z = 0,
 *
 *     where x = p*z. If this system does not have full rank,
 *     then a least squares solution is obtained. On output qrsolv
 *     also provides an upper triangular matrix s such that
 *
 *	    pT *(aT *a + d*d)*p = sT *s.
 *
 *     s is computed within qrsolv and may be of separate interest.
 *
 *     Parameters
 *
 *	n is a positive integer input variable set to the order of r.
 *
 *	r is an n by n array. On input the full upper triangle
 *	  must contain the full upper triangle of the matrix r.
 *	  On OUTPUT the full upper triangle is unaltered, and the
 *	  strict lower triangle contains the strict upper triangle
 *	  (transposed) of the upper triangular matrix s.
 *
 *	ldr is a positive integer input variable not less than n
 *	  which specifies the leading dimension of the array r.
 *
 *	ipvt is an integer input array of length n which defines the
 *	  permutation matrix p such that a*p = q*r. Column j of p
 *	  is column ipvt(j) of the identity matrix.
 *
 *	diag is an input array of length n which must contain the
 *	  diagonal elements of the matrix d.
 *
 *	qtb is an input array of length n which must contain the first
 *	  n elements of the vector (q transpose)*b.
 *
 *	x is an OUTPUT array of length n which contains the least
 *	  squares solution of the system a*x = b, d*x = 0.
 *
 *	sdiag is an OUTPUT array of length n which contains the
 *	  diagonal elements of the upper triangular matrix s.
 *
 *	wa is a work array of length n.
 *
 */
    int i, kk, j, k, nsing;
    double qtbpj, sum, temp;
    double _sin, _cos, _tan, _cot; /* local variables, not functions */

/*** qrsolv: copy r and (q transpose)*b to preserve input and initialize s.
     in particular, save the diagonal elements of r in x. ***/

    for (j = 0; j < n; j++) {
	for (i = j; i < n; i++)
	    r[j * ldr + i] = r[i * ldr + j];
	x[j] = r[j * ldr + j];
	wa[j] = qtb[j];
    }
#ifdef LMFIT_DEBUG_MESSAGES
    printf("qrsolv\n");
#endif

/*** qrsolv: eliminate the diagonal matrix d using a Givens rotation. ***/

    for (j = 0; j < n; j++) {

/*** qrsolv: prepare the row of d to be eliminated, locating the
     diagonal element using p from the qr factorization. ***/

	if (diag[ipvt[j]] == 0.)
	    goto L90;
	for (k = j; k < n; k++)
	    sdiag[k] = 0.;
	sdiag[j] = diag[ipvt[j]];

/*** qrsolv: the transformations to eliminate the row of d modify only 
     a single element of qT*b beyond the first n, which is initially 0. ***/

	qtbpj = 0.;
	for (k = j; k < n; k++) {

            /** determine a Givens rotation which eliminates the
                appropriate element in the current row of d. **/

	    if (sdiag[k] == 0.)
		continue;
	    kk = k + ldr * k;
	    if (fabs(r[kk]) < fabs(sdiag[k])) {
		_cot = r[kk] / sdiag[k];
		_sin = 1 / sqrt(1 + SQR(_cot));
		_cos = _sin * _cot;
	    } else {
		_tan = sdiag[k] / r[kk];
		_cos = 1 / sqrt(1 + SQR(_tan));
		_sin = _cos * _tan;
	    }

            /** compute the modified diagonal element of r and
                the modified element of ((q transpose)*b,0). **/

	    r[kk] = _cos * r[kk] + _sin * sdiag[k];
	    temp = _cos * wa[k] + _sin * qtbpj;
	    qtbpj = -_sin * wa[k] + _cos * qtbpj;
	    wa[k] = temp;

            /** accumulate the tranformation in the row of s. **/

	    for (i = k + 1; i < n; i++) {
		temp = _cos * r[k * ldr + i] + _sin * sdiag[i];
		sdiag[i] = -_sin * r[k * ldr + i] + _cos * sdiag[i];
		r[k * ldr + i] = temp;
	    }
	}

      L90:
        /** store the diagonal element of s and restore
            the corresponding diagonal element of r. **/

	sdiag[j] = r[j * ldr + j];
	r[j * ldr + j] = x[j];
    }

/*** qrsolv: solve the triangular system for z. if the system is
     singular, then obtain a least squares solution. ***/

    nsing = n;
    for (j = 0; j < n; j++) {
	if (sdiag[j] == 0. && nsing == n)
	    nsing = j;
	if (nsing < n)
	    wa[j] = 0;
    }

    for (j = nsing - 1; j >= 0; j--) {
	sum = 0;
	for (i = j + 1; i < nsing; i++)
	    sum += r[j * ldr + i] * wa[i];
	wa[j] = (wa[j] - sum) / sdiag[j];
    }

/*** qrsolv: permute the components of z back to components of x. ***/

    for (j = 0; j < n; j++)
	x[ipvt[j]] = wa[j];

} /*** lm_qrsolv. ***/

/**************************************************************/

int Mises_f (const gsl_vector *x, void *data, gsl_vector *f)
{
  size_t n = ((struct data *)data)->n;  // n is the nbin
  int nComp = ((struct data *)data)->nComp;
  int *nVm = ((struct data *)data)->nVm;
  double *t = ((struct data *)data)->x;
  double *y = ((struct data *)data)->y;
  double *sigma = ((struct data *)data)->sigma;

  double result;
  int i, k, j, index;
  for (i = 0; i < n; i++)
  {
    result = 0.0;
	  index = 0;
    for (k = 0; k < nComp; k++)
    {
      for (j = 0; j < nVm[k]; j++)
      {
				result += fabs((gsl_vector_get(x,index+3*j)))*exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
      }
			index += nVm[k]*3;
      //result += fabs(par[3*k+0])*exp(par[3*k+1]*(cos((x-(par[3*k+2]))*2*M_PI)-1));
      //result += par[3*k+0]*exp(par[3*k+1]*(cos((x-(par[3*k+2]))*2*M_PI)-1));
    }
    gsl_vector_set (f, i, (result - y[i])/sigma[i]);
  }

  return GSL_SUCCESS;
}

int Mises_df (const gsl_vector *x, void *data, gsl_matrix *J)
{
  size_t n = ((struct data *)data)->n;  // n is the nbin
	int nComp = ((struct data *)data)->nComp;
  int *nVm = ((struct data *)data)->nVm;
	double *t = ((struct data *)data)->x;
	double *y = ((struct data *)data)->y;
	double *sigma = ((struct data *)data)->sigma;

	int i, k, j, index;
	for (i = 0; i < n; i++)
	{
		index = 0;
    for (k = 0; k < nComp; k++)
    {
      for (j = 0; j < nVm[k]; j++)
      {
				double phase = t[i];
				double s = sigma[i];
				double e1;
			//double e1 = exp(gsl_vector_get(x,3*k+1)*(cos((t[i]-(gsl_vector_get(x,3*k+2)))*2*M_PI)-1));
				if (gsl_vector_get(x,index+3*j) >= 0)
				{
					e1 = exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
				}
				else
				{
					e1 = -exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
				}
				double e2 = (gsl_vector_get(x,index+3*j))*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1)*exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
				double e3 = (gsl_vector_get(x,index+3*j))*exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1))*(gsl_vector_get(x,index+3*j+1)*sin((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)*2*M_PI);
				gsl_matrix_set(J,i,index+3*j, e1/s);
				gsl_matrix_set(J,i,index+3*j+1, e2/s);
				gsl_matrix_set(J,i,index+3*j+2, e3/s);
			}
			index += nVm[k]*3;
		}
	}
	return GSL_SUCCESS;
}

int Mises_fdf (const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
{
	Mises_f (x, data, f);
	Mises_df (x, data, J);
	return GSL_SUCCESS;
}

// Do the non-linear fit of the components, and also output the error of parameters
void doFit_err(float *fx,float *fy, double devi, int nbin,tmplStruct *tmpl,int chan,int pol)
{
  double pval[tmpl->channel[chan].pol[pol].allVm*3];
  double datX[nbin],datY[nbin];
  double sigma[nbin];
  unsigned int i;
  int nFit;

  nFit = tmpl->channel[chan].pol[pol].allVm*3;

  // Fit for all parameters of all components
  int ivm=0;
  int icomp=0;
  for (i=0;i<tmpl->channel[chan].pol[pol].allVm;i++)
    {
      if (ivm != tmpl->channel[chan].pol[pol].comp[icomp].nVm-1)
      {
				pval[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height;
				pval[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration;
				pval[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid;
				ivm++;
      }
      else
      {
				pval[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height;
				pval[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration;
				pval[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid;
				icomp++;
				ivm=0;
      }
	  }
	// Fill up the fitting array
	for (i=0;i<nbin;i++)
    {
		datX[i] = (double)fx[i];
		datY[i] = (double)fy[i];
		sigma[i] = devi;
		//printf ("%lf\n",sigma[i]);
    }

	int nComp = tmpl->channel[chan].pol[pol].nComp;
	//global_nComp = tmpl->channel[chan].pol[pol].nComp;
	
	int nVm[nComp];
	for (i = 0; i < nComp; i++)
	{
		nVm[i] = tmpl->channel[chan].pol[pol].comp[i].nVm;
	}

	struct data dat = {nbin, nComp, nVm, datX, datY, sigma}; 

	// Initialize the solver
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int iter = 0;
	const size_t n = nbin;
	const size_t p = nFit;

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_function_fdf f;

	//double x_init[nFit] = { 1.0, 0.0, 0.0 };
	gsl_vector_view x = gsl_vector_view_array (pval, p); // pval is the x_init

	const gsl_rng_type * type;
	gsl_rng * r;
	gsl_rng_env_setup();

	type = gsl_rng_default;
	r = gsl_rng_alloc (type);

	f.f = &Mises_f;
	f.df = &Mises_df;
	f.fdf = &Mises_fdf;
	f.n = n;
	f.p = p;
	f.params = &dat;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	//print_state (iter, s);
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		//printf ("status = %s\n", gsl_strerror (status));
		//print_state (iter, s);
		if (status)
			break;
		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-4, 1e-4);
	}
	while (status == GSL_CONTINUE && iter < 500);
	gsl_multifit_covar (s->J, 0.0, covar);

	#define FIT(i) gsl_vector_get(s->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

	double chi = gsl_blas_dnrm2(s->f);
	double dof = n - p;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));
	//printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
	//printf ("A = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
	//printf ("lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
	//printf ("b = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
	//printf ("status = %s\n", gsl_strerror (status));

	// Update components from the fit
	int k, j;
  int index = 0;
  for (k = 0; k < nComp; k++)
  {
    for (j = 0; j < nVm[k]; j++)
    {
			printf("Fit results: COMP%d %g (%g) %g (%g) %g (%g)\n", k, FIT(index+3*j), c*ERR(index+3*j), FIT(index+3*j+1), c*ERR(index+3*j+1), FIT(index+3*j+2), c*ERR(index+3*j+2));
		}
		index += nVm[k]*3;
	}

  ivm=0;
  icomp=0;
  for (i=0;i<tmpl->channel[chan].pol[pol].allVm;i++)
    {
      if (ivm != tmpl->channel[chan].pol[pol].comp[icomp].nVm-1)
      {
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height = FIT(3*i);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = FIT(3*i+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = FIT(3*i+2);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err = c*ERR(3*i);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration_err = c*ERR(3*i+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid_err = c*ERR(3*i+2);
				ivm++;
      }
      else
      {
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height = FIT(3*i);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = FIT(3*i+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = FIT(3*i+2);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err = c*ERR(3*i);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration_err = c*ERR(3*i+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid_err = c*ERR(3*i+2);
				icomp++;
				ivm=0;
      }
    }

	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	gsl_rng_free (r);
}

int Mises_gf (const gsl_vector *x, void *datag, gsl_vector *f)
{
  size_t n = ((struct datag *)datag)->n;  // n is the nbin
  int nComp = ((struct datag *)datag)->nComp;
  int *nVm = ((struct datag *)datag)->nVm;
  double *t = ((struct datag *)datag)->x;
  double *y = ((struct datag *)datag)->y;
  double *sigma = ((struct datag *)datag)->sigma;
  double *par = ((struct datag *)datag)->par;

  double result;
  int i, k, j, index;
  for (i = 0; i < n; i++)
  {
    result = 0.0;
	  index = 0;
    for (k = 0; k < nComp; k++)
    {
      for (j = 0; j < nVm[k]; j++)
      {
				result += fabs(gsl_vector_get(x,2*k))*fabs(par[index+3*j])*exp(par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,2*k+1))*2*M_PI)-1));
      }
			index += nVm[k]*3;
      //result += fabs(par[3*k+0])*exp(par[3*k+1]*(cos((x-(par[3*k+2]))*2*M_PI)-1));
      //result += par[3*k+0]*exp(par[3*k+1]*(cos((x-(par[3*k+2]))*2*M_PI)-1));
    }
    gsl_vector_set (f, i, (result - y[i])/sigma[i]);
  }

  return GSL_SUCCESS;
}

int Mises_gdf (const gsl_vector *x, void *datag, gsl_matrix *J)
{
  size_t n = ((struct datag *)datag)->n;  // n is the nbin
	int nComp = ((struct datag *)datag)->nComp;
  int *nVm = ((struct datag *)datag)->nVm;
	double *t = ((struct datag *)datag)->x;
	double *y = ((struct datag *)datag)->y;
	double *sigma = ((struct datag *)datag)->sigma;
  double *par = ((struct datag *)datag)->par;

	int i, k, j, index;
	for (i = 0; i < n; i++)
	{
		index = 0;
		double s = sigma[i];
    for (k = 0; k < nComp; k++)
    {
			double e1 = 0.0;
			double e2 = 0.0;
      for (j = 0; j < nVm[k]; j++)
      {
			//double e1 = exp(gsl_vector_get(x,3*k+1)*(cos((t[i]-(gsl_vector_get(x,3*k+2)))*2*M_PI)-1));
				if (gsl_vector_get(x,2*k) >= 0)
				{
					//e1 = exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
					e1 += fabs(par[index+3*j])*exp(par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,2*k+1))*2*M_PI)-1));
				}
				else
				{
					//e1 = -exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
					e1 += -fabs(par[index+3*j])*exp(par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,2*k+1))*2*M_PI)-1));
				}
				//double e2 = (gsl_vector_get(x,index+3*j))*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1)*exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
				e2 += fabs(gsl_vector_get(x,2*k))*fabs(par[index+3*j])*exp(par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,2*k+1))*2*M_PI)-1))*(2*M_PI*par[index+3*j+1]*sin((t[i]-par[index+3*j+2]-gsl_vector_get(x,2*k+1))*2*M_PI));
			}
			gsl_matrix_set(J,i,2*k, e1/s);
			gsl_matrix_set(J,i,2*k+1, e2/s);
			index += nVm[k]*3;
		}
	}
	return GSL_SUCCESS;
}

int Mises_gfdf (const gsl_vector *x, void *datag, gsl_vector *f, gsl_matrix *J)
{
	Mises_gf (x, datag, f);
	Mises_gdf (x, datag, J);
	return GSL_SUCCESS;
}

// Do the non-linear fit of the components, and also output the error of parameters
void doFit_gerr(float *fx,float *fy, double devi, int nbin,tmplStruct *tmpl,int chan,int pol)
{
  double pval[tmpl->channel[chan].pol[pol].nComp*2];
  double par[tmpl->channel[chan].pol[pol].allVm*3];
  double par_err[tmpl->channel[chan].pol[pol].allVm*3];
  double datX[nbin],datY[nbin];
  double sigma[nbin];
  unsigned int i;
  int nFit;

  nFit = tmpl->channel[chan].pol[pol].nComp*2;

  // Fit for all parameters of all components
  int ivm=0;
  int icomp=0;
  for (i=0;i<tmpl->channel[chan].pol[pol].allVm;i++)
    {
      if (ivm != tmpl->channel[chan].pol[pol].comp[icomp].nVm-1)
      {
				par[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height;
				par[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration;
				par[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid;
				par_err[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err;
				par_err[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration_err;
				par_err[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid_err;
				ivm++;
      }
      else
      {
				par[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height;
				par[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration;
				par[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid;
				par_err[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err;
				par_err[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration_err;
				par_err[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid_err;
				icomp++;
				ivm=0;
      }
	  }

  for (i=0;i<tmpl->channel[chan].pol[pol].nComp;i++)
    {
			pval[2*i] = 1.0;
			pval[2*i+1] = 0.0;
		}

	// Fill up the fitting array
	for (i=0;i<nbin;i++)
    {
		datX[i] = (double)fx[i];
		datY[i] = (double)fy[i];
		sigma[i] = devi;
		//printf ("%lf\n",sigma[i]);
    }

	int nComp = tmpl->channel[chan].pol[pol].nComp;
	//global_nComp = tmpl->channel[chan].pol[pol].nComp;
	
	int nVm[nComp];
	for (i = 0; i < nComp; i++)
	{
		nVm[i] = tmpl->channel[chan].pol[pol].comp[i].nVm;
	}

	struct datag dat = {nbin, nComp, nVm, datX, datY, sigma, par}; 

	// Initialize the solver
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int iter = 0;
	const size_t n = nbin;
	const size_t p = nFit;

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_function_fdf f;

	//double x_init[nFit] = { 1.0, 0.0, 0.0 };
	gsl_vector_view x = gsl_vector_view_array (pval, p); // pval is the x_init

	const gsl_rng_type * type;
	gsl_rng * r;
	gsl_rng_env_setup();

	type = gsl_rng_default;
	r = gsl_rng_alloc (type);

	f.f = &Mises_gf;
	f.df = &Mises_gdf;
	f.fdf = &Mises_gfdf;
	f.n = n;
	f.p = p;
	f.params = &dat;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	printf("/////////////////\n");
	//print_state (iter, s);
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		//printf ("status = %s\n", gsl_strerror (status));
		//print_state (iter, s);
		if (status)
			break;
		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-4, 1e-4);
	}
	while (status == GSL_CONTINUE && iter < 500);
	gsl_multifit_covar (s->J, 0.0, covar);

	#define FIT(i) gsl_vector_get(s->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

	double chi = gsl_blas_dnrm2(s->f);
	double dof = n - p;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));
	//printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
	//printf ("A = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
	//printf ("lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
	//printf ("b = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
	//printf ("status = %s\n", gsl_strerror (status));

	// Update components from the fit
	int k, j;
  int index = 0;
  for (k = 0; k < nComp; k++)
  {
    for (j = 0; j < nVm[k]; j++)
    {
			printf("Fit results: COMP%d %g (%g) %g (%g) %g (%g)\n", k, FIT(2*k)*par[index+3*j], c*ERR(2*k), par[index+3*j+1], par_err[index+3*j+1], FIT(2*k+1)+par[index+3*j+2], c*ERR(2*k+1));
		}
		index += nVm[k]*3;
	}

  ivm=0;
  icomp=0;
  for (i=0;i<tmpl->channel[chan].pol[pol].allVm;i++)
    {
      if (ivm != tmpl->channel[chan].pol[pol].comp[icomp].nVm-1)
      {
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height = par[3*i]*FIT(2*icomp);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = par[3*i+1];
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = par[3*i+2]+FIT(2*icomp+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err = c*ERR(2*icomp);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = par_err[3*i+1];
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = c*ERR(2*icomp+1);
				ivm++;
      }
      else
      {
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height = par[3*i]*FIT(2*icomp);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = par[3*i+1];
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = par[3*i+2]+FIT(2*icomp+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err = c*ERR(2*icomp);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = par_err[3*i+1];
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = c*ERR(2*icomp+1);
				icomp++;
				ivm=0;
      }
    }

	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	gsl_rng_free (r);
}

int Mises_Gf (const gsl_vector *x, void *datag, gsl_vector *f)
{
  size_t n = ((struct datag *)datag)->n;  // n is the nbin
  int nComp = ((struct datag *)datag)->nComp;
  int *nVm = ((struct datag *)datag)->nVm;
  double *t = ((struct datag *)datag)->x;
  double *y = ((struct datag *)datag)->y;
  double *sigma = ((struct datag *)datag)->sigma;
  double *par = ((struct datag *)datag)->par;

  double result;
  int i, k, j, index;
  for (i = 0; i < n; i++)
  {
    result = 0.0;
	  index = 0;
    for (k = 0; k < nComp; k++)
    {
      for (j = 0; j < nVm[k]; j++)
      {
				result += fabs(gsl_vector_get(x,3*k))*fabs(par[index+3*j])*exp(gsl_vector_get(x,3*k+1)*par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,3*k+2))*2*M_PI)-1));
      }
			index += nVm[k]*3;
      //result += fabs(par[3*k+0])*exp(par[3*k+1]*(cos((x-(par[3*k+2]))*2*M_PI)-1));
      //result += par[3*k+0]*exp(par[3*k+1]*(cos((x-(par[3*k+2]))*2*M_PI)-1));
    }
    gsl_vector_set (f, i, (result - y[i])/sigma[i]);
  }

  return GSL_SUCCESS;
}

int Mises_Gdf (const gsl_vector *x, void *datag, gsl_matrix *J)
{
  size_t n = ((struct datag *)datag)->n;  // n is the nbin
	int nComp = ((struct datag *)datag)->nComp;
  int *nVm = ((struct datag *)datag)->nVm;
	double *t = ((struct datag *)datag)->x;
	double *y = ((struct datag *)datag)->y;
	double *sigma = ((struct datag *)datag)->sigma;
  double *par = ((struct datag *)datag)->par;

	int i, k, j, index;
	for (i = 0; i < n; i++)
	{
		index = 0;
		double s = sigma[i];
    for (k = 0; k < nComp; k++)
    {
			double e1 = 0.0;
			double e2 = 0.0;
			double e3 = 0.0;
      for (j = 0; j < nVm[k]; j++)
      {
			//double e1 = exp(gsl_vector_get(x,3*k+1)*(cos((t[i]-(gsl_vector_get(x,3*k+2)))*2*M_PI)-1));
				if (gsl_vector_get(x,3*k) >= 0)
				{
					//e1 = exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
					e1 += fabs(par[index+3*j])*exp(gsl_vector_get(x,3*k+1)*par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,3*k+2))*2*M_PI)-1));
				}
				else
				{
					//e1 = -exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
					e1 += -fabs(par[index+3*j])*exp(gsl_vector_get(x,3*k+1)*par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,3*k+2))*2*M_PI)-1));
				}
				//double e2 = (gsl_vector_get(x,index+3*j))*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1)*exp(gsl_vector_get(x,index+3*j+1)*(cos((t[i]-(gsl_vector_get(x,index+3*j+2)))*2*M_PI)-1));
				e2 += fabs(gsl_vector_get(x,3*k))*fabs(par[index+3*j])*exp(gsl_vector_get(x,3*k+1)*par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,3*k+2))*2*M_PI)-1))*(par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,3*k+1))*2*M_PI)-1));
				e3 += fabs(gsl_vector_get(x,3*k))*fabs(par[index+3*j])*exp(gsl_vector_get(x,3*k+1)*par[index+3*j+1]*(cos((t[i]-par[index+3*j+2]-gsl_vector_get(x,3*k+2))*2*M_PI)-1))*(2*M_PI*gsl_vector_get(x,3*k+1)*par[index+3*j+1]*sin((t[i]-par[index+3*j+2]-gsl_vector_get(x,3*k+2))*2*M_PI));
			}
			gsl_matrix_set(J,i,3*k, e1/s);
			gsl_matrix_set(J,i,3*k+1, e2/s);
			gsl_matrix_set(J,i,3*k+2, e3/s);
			index += nVm[k]*3;
		}
	}
	return GSL_SUCCESS;
}

int Mises_Gfdf (const gsl_vector *x, void *datag, gsl_vector *f, gsl_matrix *J)
{
	Mises_Gf (x, datag, f);
	Mises_Gdf (x, datag, J);
	return GSL_SUCCESS;
}

// Do the non-linear fit of the components, and also output the error of parameters
void doFit_Gerr(float *fx,float *fy, double devi, int nbin,tmplStruct *tmpl,int chan,int pol)
{
  double pval[tmpl->channel[chan].pol[pol].nComp*3];
  double par[tmpl->channel[chan].pol[pol].allVm*3];
  double par_err[tmpl->channel[chan].pol[pol].allVm*3];
  double datX[nbin],datY[nbin];
  double sigma[nbin];
  unsigned int i;
  int nFit;

  nFit = tmpl->channel[chan].pol[pol].nComp*3;

  // Fit for all parameters of all components
  int ivm=0;
  int icomp=0;
  for (i=0;i<tmpl->channel[chan].pol[pol].allVm;i++)
    {
      if (ivm != tmpl->channel[chan].pol[pol].comp[icomp].nVm-1)
      {
				par[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height;
				par[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration;
				par[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid;
				par_err[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err;
				par_err[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration_err;
				par_err[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid_err;
				ivm++;
      }
      else
      {
				par[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height;
				par[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration;
				par[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid;
				par_err[3*i] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err;
				par_err[3*i+1] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration_err;
				par_err[3*i+2] = tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid_err;
				icomp++;
				ivm=0;
      }
	  }

  for (i=0;i<tmpl->channel[chan].pol[pol].nComp;i++)
    {
			pval[3*i] = 1.0;
			pval[3*i+1] = 1.0;
			pval[3*i+2] = 0.0;
		}

	// Fill up the fitting array
	for (i=0;i<nbin;i++)
    {
		datX[i] = (double)fx[i];
		datY[i] = (double)fy[i];
		sigma[i] = devi;
		//printf ("%lf\n",sigma[i]);
    }

	int nComp = tmpl->channel[chan].pol[pol].nComp;
	//global_nComp = tmpl->channel[chan].pol[pol].nComp;
	
	int nVm[nComp];
	for (i = 0; i < nComp; i++)
	{
		nVm[i] = tmpl->channel[chan].pol[pol].comp[i].nVm;
	}

	struct datag dat = {nbin, nComp, nVm, datX, datY, sigma, par}; 

	// Initialize the solver
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int iter = 0;
	const size_t n = nbin;
	const size_t p = nFit;

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_function_fdf f;

	//double x_init[nFit] = { 1.0, 0.0, 0.0 };
	gsl_vector_view x = gsl_vector_view_array (pval, p); // pval is the x_init

	const gsl_rng_type * type;
	gsl_rng * r;
	gsl_rng_env_setup();

	type = gsl_rng_default;
	r = gsl_rng_alloc (type);

	f.f = &Mises_Gf;
	f.df = &Mises_Gdf;
	f.fdf = &Mises_Gfdf;
	f.n = n;
	f.p = p;
	f.params = &dat;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	printf("/////////////////\n");
	//print_state (iter, s);
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		//printf ("status = %s\n", gsl_strerror (status));
		//print_state (iter, s);
		if (status)
			break;
		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-4, 1e-4);
	}
	while (status == GSL_CONTINUE && iter < 500);
	gsl_multifit_covar (s->J, 0.0, covar);

	#define FIT(i) gsl_vector_get(s->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

	double chi = gsl_blas_dnrm2(s->f);
	double dof = n - p;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));
	//printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
	//printf ("A = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
	//printf ("lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
	//printf ("b = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
	//printf ("status = %s\n", gsl_strerror (status));

	// Update components from the fit
	int k, j;
  int index = 0;
  for (k = 0; k < nComp; k++)
  {
    for (j = 0; j < nVm[k]; j++)
    {
			printf("Fit results: COMP%d %g (%g) %g (%g) %g (%g)\n", k, FIT(3*k)*par[index+3*j], c*ERR(3*k), FIT(3*k+1)*par[index+3*j+1], c*ERR(3*k+1), FIT(3*k+2)+par[index+3*j+2], c*ERR(3*k+2));
		}
		index += nVm[k]*3;
	}

  ivm=0;
  icomp=0;
  for (i=0;i<tmpl->channel[chan].pol[pol].allVm;i++)
    {
      if (ivm != tmpl->channel[chan].pol[pol].comp[icomp].nVm-1)
      {
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height = par[3*i]*FIT(3*icomp);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = par[3*i+1]*FIT(3*icomp+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = par[3*i+2]+FIT(3*icomp+2);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err = c*ERR(3*icomp);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration_err = c*ERR(3*icomp+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid_err = c*ERR(3*icomp+2);
				ivm++;
      }
      else
      {
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height = par[3*i]*FIT(3*icomp);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration = par[3*i+1]*FIT(3*icomp+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid = par[3*i+2]+FIT(3*icomp+2);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].height_err = c*ERR(3*icomp);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].concentration_err = c*ERR(3*icomp+1);
				tmpl->channel[chan].pol[pol].comp[icomp].vonMises[ivm].centroid_err = c*ERR(3*icomp+2);
				icomp++;
				ivm=0;
      }
    }

	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	gsl_rng_free (r);
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
	printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
			"|f(x)| = %g\n",
			iter,
			gsl_vector_get (s->x, 0),
			gsl_vector_get (s->x, 1),
			gsl_vector_get (s->x, 2),
			gsl_blas_dnrm2 (s->f));
}

