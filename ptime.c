
// SHOULD TRY WITH ONLY ~4 FREQUENCY CHANNELS -- TAKES TOO LONG!

// To do:
// 1. Must deallocate memory
// 3. Must calculate period
// 4. Must calculate frequency
// 5. Must add in frequency dependence
// 6. Deal with subintegrations
// 7. Must implement Taylor algorithm
// 8. Must have sensibe output file name for ToAs.
// 9. Include filenames on plots
// 10. Should correct set frequency if de-dispersed
// 11. Should have plot of difference between template and profile (include in log)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ptimeLib.h"
#include <cpgplot.h>
#include "readPfits.h"
#include "fitsio.h"
#include "tempo2pred.h"

// Linear LSQ fitting
#include "TKsvd.h"
#include "TKmatrix.h"
#include "T2toolkit.h"
#include "TKfit.h"


// Program to determine pulse ToAs
typedef struct toaStruct {
  double gof; // Goodness of fit
  double dphi; // Phase offset
  double dphiErr; // Error in phase offset
  char fname[1024]; // Filename
} toaStruct;

void processFile(char *fname,tmplStruct *tmpl,FILE *fout,int logit);
void getTOA_alg1(ptime_observation *obs,pheader *header,tmplStruct *tmpl,toaStruct *toa,FILE *fout_log);
void findMinMax(int n,float *y,float *miny,float *maxy);
void calcArrivalTime(ptime_observation *obs, pheader *header,tmplStruct *tmpl,toaStruct *toa,double *offs_sub,double *datFreq,long double period,FILE *fout,FILE *fout_log);
void fitFunc(double x,double *p,tmplStruct *tmpl,int nParam,double phi,int ival,int jval,int kval);

int main(int argc,char *argv[])
{
  tmplStruct tmpl;
  int i;
  char tmplFile[1024];
  char fname[128];
  char outName[128] = "ptime.tim";
  FILE *fout;
  int logit=1;

  for (i=0;i<argc;i++)
    {

      if (strcmp(argv[i],"-f")==0)
	strcpy(tmplFile,argv[++i]);
    }


  initialiseTemplate(&tmpl);
  printf("Reading template\n");
  readTemplate(tmplFile,&tmpl); 
  printf("Complete reading template\n");
  fout = fopen(outName,"w");
  fprintf(fout,"FORMAT 1\n");
  fflush(fout);
  // Now go through all the observations to process
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0){
	i++;
      }  else {
	strcpy(fname,argv[i]);
	processFile(fname,&tmpl,fout,logit);
      }	
    }
  fclose(fout);
  // Should deallocate template memory

}

void processFile(char *fname,tmplStruct *tmpl,FILE *fout,int logit)
{
  pheader *header;
  ptime_observation *obs;
  toaStruct toa;
  fitsfile *fp;
  int i;
  float baselineFrac = 0.1;
  int baselineType = 1;
  int toaAlgorithm=1;
  double dPhi;
  double *offs_sub;
  double *datFreq;
  int nSub = 1;
  FILE *fout_log;
  char logFname[128];
  long double period;

  sprintf(logFname,"log.ptime.%s",fname);
  fout_log = fopen(logFname,"w");
  offs_sub = (double *)malloc(sizeof(double)*nSub);
  obs = (ptime_observation *)malloc(sizeof(ptime_observation));
  header = (pheader *)malloc(sizeof(pheader));
  fp = openFitsFile(fname);
  loadPrimaryHeader(fp,header);

  datFreq= (double *)malloc(sizeof(double)*header->nchan);
  allocateObsMemory(obs,header);
  readData(obs,header,fp);
  // Now get an array of subintegration times
  readSubintOffs(obs,offs_sub,fp);
  // Get an array of frequencies
  readDatFreq(obs,datFreq,fp,header->nchan);

  printf("offs_sub = %g\n",offs_sub[0]);
  // Calculate the period
  {
    FILE *pred_out;
    int status=0;
    long nrows,row;
    char nval[128]="UNKNOWN";
    int anynul=0;
    long double mjd0;
    long double frequency;
    char **line;
    T2Predictor pred;
    int ret;

    line = (char **)malloc(sizeof(char *));
    line[0] = (char *)malloc(sizeof(char)*1024);     


    if (!(pred_out = fopen("ptime.pred","w"))){
      printf("Unable to open file >%s<\n","ptime.pred");
    }
    fits_movnam_hdu(fp,BINARY_TBL,(char *)"T2PREDICT",1,&status);
    if (status)
      {
	printf("No predictor table in PSRFITS file\n");
	status=0;
      }
    fits_get_num_rows(fp,&nrows,&status);
    printf("NROWS = %d\n",nrows);
    for (row = 1; row <= nrows ; row++){
      fits_read_col_str(fp,1,row,1,1,nval,line,&anynul,&status);
      printf("Have read %s\n",line[0]);
      fprintf(pred_out,"%s\n",line[0]);
    }
    free(line[0]);
    free(line);

    fclose(pred_out);

    T2Predictor_Init(&pred);  // prepare the predictor                                                                                     
    if (ret=T2Predictor_Read(&pred,(char *)"ptime.pred"))
      {
	printf("Error: unable to read predictor\n");
	exit(1);
      }
    mjd0 = header->imjd + (header->smjd + header->stt_offs)/86400.0L;
    frequency = datFreq[0];  //(freq in MHz) WHAT SHOULD THIS BE SET TO!!
    period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,frequency);
    if (period==-1){
      printf("Error: period returned from predictor = -1. This cannot be correct.\n");
      exit(1);
    }
    T2Predictor_Destroy(&pred);
    printf("Period = %.15Lf %.15Lf %.15Lf\n",period,mjd0,frequency);
  }

  closeFitsFile(fp);
  // Must remove a baseline
  removeBaseline(obs,header,baselineType,baselineFrac);
  printf("Here\n");
  // Now get the shift
  if (toaAlgorithm == 1)
    getTOA_alg1(obs,header,tmpl,&toa,fout_log);
  strcpy(toa.fname,fname);



  // Calculate the arrival time
  calcArrivalTime(obs,header,tmpl,&toa,offs_sub,datFreq,period,fout,fout_log);
  free(header);
  free(obs);
  free(offs_sub);
  free(datFreq);
  fclose(fout_log);
}


// MUST PASS IN THE OFFSET, ERROR, GOF AND FLAG LINES ETc.
void calcArrivalTime(ptime_observation *obs, pheader *header,tmplStruct *tmpl,toaStruct *toa,double *offs_sub,double *datFreq,long double period,FILE *fout,FILE *fout_log)
{
  long double arrivalTime;
  //  long double period; // = 0.016451278438646; // MUST FIX
  long double error;
  int sub=0;

  // Must calculate period


  // Now calculate arrival time
  arrivalTime = header->imjd + (offs_sub[sub] + header->smjd + header->stt_offs + toa->dphi*period)/86400.0L;
  error = toa->dphiErr * period * 1e6; // Error in microseconds

  printf("%s %.3f %.15Lf %.3Lf %s\n",toa->fname,datFreq[0],arrivalTime,error,header->telescope);
  fprintf(fout,"%s %.3f %.15Lf %.4Lf %s\n",toa->fname,datFreq[0],arrivalTime,error,header->telescope);
  fprintf(fout_log,"SUB: %d\n",sub);
  fprintf(fout_log,"IMJD: %d\n",header->imjd);
  fprintf(fout_log,"SMJD: %.5f\n",header->smjd);
  fprintf(fout_log,"STT_OFFS: %.5f\n",header->stt_offs);
  fprintf(fout_log,"period: %.15Lf\n",period);
  fprintf(fout_log,"OFF_SUB[sub]: %.5f\n",offs_sub[sub]);

  fprintf(fout_log,"%s %.3f %.15Lf %.4Lf %s\n",toa->fname,datFreq[0],arrivalTime,error,header->telescope);
  fflush(fout);
}

void getTOA_alg1(ptime_observation *obs,pheader *header,tmplStruct *tmpl,toaStruct *toa,FILE *fout_log)
{
  int i,j,k;
  int nbin = header->nbin;
  int nchan = header->nchan;
  int npol = header->npol;
  double chisq;
  double diffVals[nbin];
  double phiRot = 0;
  double step1 = 1.0/nbin/2.0;
  double step = step1;
  double baseline = 0;
  double scale = 1;
  double tmplEval;
  int chan = 0;
  int pol = 0;
  double phi,bestPhi;
  double phi0,phi1,bestChisq;
  int it;
  double chisqVals[nbin*2];
  double phiVals[nbin*2];
  int nChisqVals,ibest;
  int iterateAgain;
  int maxIterations = 10;
  int plotOutput=0;
  double error;
  double *fitX;
  double *fitY;
  double *fitE;
  int    *fitI,*fitJ,*fitK;
  int nfit = npol+1; // Up to 4 baselines per profile + 1 scaling factor
  double outputParams_v[nfit];
  double outputParams_e[nfit];
  double results_v[npol*nchan+nchan];
  double results_e[npol*nchan+nchan];
  int weight = 1;
  double bestParameters[npol*nchan+nchan];
  double chisqTot;

  // Covariance matrix
  double **cvm;
  cvm = (double **)malloc(sizeof(double*)*nfit);
  for (i=0;i<nfit;i++)
    cvm[i] = (double *)malloc(sizeof(double)*nfit);
  if (!(fitX = (double *)malloc(sizeof(double)*npol*nchan*nbin))){
    printf("Unable to allocate enough memory for fitX\n");
    exit(1);
  }
  if (!(fitY = (double *)malloc(sizeof(double)*npol*nchan*nbin))){
    printf("Unable to allocate enough memory for fitY\n");
    exit(1);
  }
  if (!(fitE = (double *)malloc(sizeof(double)*npol*nchan*nbin))){
    printf("Unable to allocate enough memory for fitE\n");
    exit(1);
  }
  if (!(fitI = (int *)malloc(sizeof(int)*npol*nchan*nbin))){
    printf("Unable to allocate enough memory for fitI\n");
    exit(1);
  }
  if (!(fitJ = (int *)malloc(sizeof(int)*npol*nchan*nbin))){
    printf("Unable to allocate enough memory for fitJ\n");
    exit(1);
  }
  if (!(fitK = (int *)malloc(sizeof(int)*npol*nchan*nbin))){
    printf("Unable to allocate enough memory for fitK\n");
    exit(1);
  }

  printf("Baseline sdev = %g, mean = %g\n",obs->chan[0].pol[0].sdev,obs->chan[0].pol[0].baselineVal);
  printf("On the next line\n");
  printf("npol = %d \n",npol);
  printf("Doing fit\n");
  it = 0;
  do {
    printf("Iteration %d\n",it+1);
    if (it == 0) {
      phi0 = -0.5;
      phi1 = 0.5;
    } else {
      phi0 = bestPhi - step;
      phi1 = bestPhi + step;
      step/=(double)10.0;
    }
    nChisqVals = 0;
    for (phiRot = phi0;phiRot < phi1;phiRot += step)
      {
	printf("Complete: %.1f percent\n",(phiRot-phi0)/(phi1-phi0)*100);
	//	phiRot = 0.0;
	// Least squares fit for baseline and amplitude at given phiRot
	//	printf("Doing fit\n");
	chisqTot=0;
	for (j=0;j<nchan;j++){
	  for (i=0;i<npol;i++){
	    for (k=0;k<nbin;k++){
	      //	printf("Setting %d %d %d %d %d\n",i,j,k,npol*nchan*nbin,i*(nchan*nbin)+j*nbin+k);
	      fitI[i*nbin+k] = i;
	      fitJ[i*nbin+k] = j;
	      fitK[i*nbin+k] = k;
	      fitX[i*nbin+k] = (double)k/(double)nbin; 
	      //	printf("This far\n");
	      //	printf("Searching for %g\n",obs->chan[j].pol[i].val[k]);
	      fitY[i*nbin+k] = obs->chan[j].pol[i].val[k];
	      fitE[i*nbin+k] = obs->chan[j].pol[i].sdev; 
	    }
	  }

	  TKleastSquares_svd(fitX,fitY,fitE,fitI,fitJ,fitK,npol*nbin,outputParams_v,outputParams_e,nfit,cvm, &chisq, fitFunc, tmpl, weight,phiRot);
	  chisqTot += chisq;
	  for (i=0;i<npol+1;i++)
	    {
	      results_v[(npol+1)*j + i] = outputParams_v[i];
	      results_e[(npol+1)*j + i] = outputParams_e[i];
	    }
	  //	  for (i=0;i<npol*nbin;i++)
	  //	    {
	  //	      printf("FitVals = %g %g %g\n",fitX[i],fitY[i],fitE[i]);
	  //	    }
	  //	  for (i=0;i<nbin;i++)
	  //	    printf("Best %g %g\n",fitX[i],outputParams_v[4]*evaluateTemplateChannel(tmpl,fitX[i],j,0,phiRot)+outputParams_v[0]);
	  //	  for (i=0;i<nbin;i++)
	  //	    printf("Best %g %g\n",fitX[i],outputParams_v[4]*evaluateTemplateChannel(tmpl,fitX[i],j,1,phiRot)+outputParams_v[1]);
	  //	  for (i=0;i<nbin;i++)
	  //	    printf("Best %g %g\n",fitX[i],outputParams_v[4]*evaluateTemplateChannel(tmpl,fitX[i],j,2,phiRot)+outputParams_v[2]);
	  //	  for (i=0;i<nbin;i++)
	  //	    printf("Best %g %g\n",fitX[i],outputParams_v[4]*evaluateTemplateChannel(tmpl,fitX[i],j,3,phiRot)+outputParams_v[3]);

	  //	for (i=0;i<nfit;i++){
	  //	  printf("%d %g %g\n",i,outputParams_v[i],outputParams_e[i]);


	}

	//	printf("Done fit\n");
	//	baseline = outputParams_v[0];
	//	scale = outputParams_v[1];
	//	for (i=0;i<nfit;i++){
	//	  printf("%d %g %g\n",i,outputParams_v[i],outputParams_e[i]);
	//	}
	//	exit(1);	
	chisqVals[nChisqVals] = chisqTot;
	phiVals[nChisqVals] = phiRot;
	if (nChisqVals==0){
	  bestPhi = phiRot;
	  bestChisq = chisqTot;
	  for (i=0;i<nchan*npol+nchan;i++){
	    bestParameters[i] = outputParams_v[i];
	  }
	  ibest = nChisqVals;
	} else {
	  if (bestChisq > chisqTot){
	    bestChisq = chisqTot;
	    bestPhi = phiRot;
	    ibest = nChisqVals;
	    for (i=0;i<nchan*npol+nchan;i++){
	      bestParameters[i] = outputParams_v[i];
	    }
	  }
	}
	//	printf("Chisq = %g\n",chisq);
	//	exit(1);
	nChisqVals++;
      }
    printf("nvals =%d\n",nChisqVals);
    // Should check if we need to iterate again - do check based on how chisq is changing - i.e, must get a good measure of the chisq increasing by 1
    iterateAgain = 1;
    for (i=ibest+1;i<nChisqVals;i++)
      {
	if (chisqVals[i] < bestChisq + 1){
	  iterateAgain = 0;
	  break;
	}
      }
    it++;
  } while (iterateAgain == 1 && it < maxIterations);
  printf("ibest = %d, nChisqVals = %d, bestPhi = %g\n",ibest,nChisqVals,bestPhi);
  //  exit(1);
  {
    int foundStart = 0;
    double start,end;

    // Should think how to improve this method
    for (i=0;i<nChisqVals;i++)
      {
	fprintf(fout_log,"chisqVals %g %g\n",phiVals[i],chisqVals[i]);
	if (foundStart==0 && chisqVals[i] <= bestChisq + 1) 
	  {
	    foundStart = 1;
	    start = phiVals[i]; 
	  }
	else if (foundStart == 1 && chisqVals[i] >= bestChisq + 1)
	  {
	    end = phiVals[i];
	    break;
	  }
      }
    error = (end-start)/2.0;
  }
    
  printf("Number of iterations = %d\n",it);
  printf("bestPhi = %g\n",bestPhi);
  printf("bestChisq = %g\n",bestChisq);
  for (i=0;i<nchan*npol+nchan;i++){
    printf("Best parameter: %d %g\n",i,bestParameters[i]);
  }

  fprintf(fout_log,"Number of iterations = %d\n",it);
  fprintf(fout_log,"bestPhi = %g\n",bestPhi);
  fprintf(fout_log,"bestChisq = %g\n",bestChisq);
  for (i=0;i<nchan*npol+nchan;i++){
    fprintf(fout_log,"Best parameter: %d %g\n",i,bestParameters[i]);
  }

  if (plotOutput == 1){
    float fx[nbin*2],fy[nbin*2],ft[nbin*2],dy[nbin*2];
    float miny,maxy;
    float miny2,maxy2;

    for (i=0;i<nbin;i++)
      {
	fx[i] = (float)i/(float)nbin;
	fx[i+nbin] = fx[i]+1;
	//	tmplEval = bestScale*evaluateTemplateChannel(tmpl,fx[i],chan,pol,bestPhi)+bestBaseline;
	
	fy[i] = obs->chan[chan].pol[pol].val[i];
	fy[i+nbin] = fy[i];
	ft[i] = tmplEval;
	ft[i+nbin] = ft[i];
	dy[i] = fy[i]-ft[i];
	dy[i+nbin] = dy[i];
      }
    findMinMax(nbin,fy,&miny,&maxy);
    findMinMax(nbin,dy,&miny2,&maxy2);

    cpgbeg(0,"/xs",1,1);
    cpgsvp(0.1,0.9,0.5,0.9);
    cpgeras();
    cpgswin(0,2,miny,maxy);
    cpgbox("BCTS",0.0,0,"BCTSN",0.0,0);


    cpgbin(nbin*2,fx,fy,1);
    cpgsci(2); cpgline(nbin*2,fx,ft); cpgsci(1);

    cpgsvp(0.1,0.9,0.15,0.5);
    cpgswin(0,2,miny2,maxy2);
    cpgbox("BCTSN",0.0,0,"BCTSN",0.0,0);
    cpglab("Phase","prof-tmpl","");
    cpgbin(nbin*2,fx,dy,1);



    cpgend();
  }
  toa->dphi = bestPhi;
  toa->dphiErr = error;

  for (i=0;i<nfit;i++)
    free(cvm[i]);
  free(cvm);
  free(fitX); free(fitY); free(fitE); free(fitI); free(fitJ); free(fitK);
  return;
}

void findMinMax(int n,float *y,float *miny,float *maxy)
{
  int i;
  *miny = *maxy = y[0];
  for (i=0;i<n;i++)
    {
      if (y[i] > *maxy) *maxy = y[i];
      if (y[i] < *miny) *miny = y[i];
    }
}

void fitFunc(double x,double *p,tmplStruct *tmpl,int nParam,double phi,int ival,int jval,int kval)
{
  int chan;
  int pol;
  double phiRot = phi;
  int i,j;

  //  printf("ival = %d %d %d %d\n",ival,jval,kval,nParam);
  chan = jval;
  pol = ival;
  //  printf("In here with %g %g\n",x,phi);
  for (i=0;i<tmpl->channel[0].nstokes;i++)
    {
      if (ival == i)
	p[i] = 1;
      else
	p[i] = 0;
      
    }

  p[tmpl->channel[0].nstokes] = evaluateTemplateChannel(tmpl,x,chan,pol,phiRot);


}
