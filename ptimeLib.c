#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ptimeLib.h"


// Initialises the template, but does not allocate memory for the profiles
// This routine should be called at the start of the program
void initialiseTemplate(tmplStruct *tmpl)
{
  strcpy(tmpl->dte,"UNSET");
  strcpy(tmpl->user,"UNSET");
  tmpl->templateVersion = 0;
  strcpy(tmpl->source,"UNSET");
  strcpy(tmpl->profileFile,"UNSET");
  strcpy(tmpl->units,"UNSET");
  tmpl->dedispersed = 0;
  tmpl->nchan = 0;
  tmpl->channelMemoryAllocated=0;
  tmpl->nChannelAllocated=0;
}

// Reads a template from disk
void readTemplate(char *file,tmplStruct *tmpl)
{
  FILE *fin;
  int i;
  char line[4096];
  char firstword[4096];
  char dummy[4096];
  int nchan=-1;
  int nstokes=1;
  int chan,stokes,comp;
  char stokesStr[1024];
  double f1,f2;

  // Read primary header
  if (!(fin = fopen(file,"r"))){
    printf("Unable to open file: >%s<\n",file);
    exit(1);
  }
  while (!feof(fin))
    {
      if (fgets(line,4096,fin) != NULL){
	if (line[0] == '#') // Comment line
	  {
	    // Do nothing
	  }
	else {
	  sscanf(line,"%s",firstword);
	  if (strcasecmp(firstword,"TEMPLATE_VERSION:")==0)
	    sscanf(line,"%s %f",dummy,&(tmpl->templateVersion));
	  else if (strcasecmp(firstword,"SOURCE:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->source));
	  else if (strcasecmp(firstword,"PROFILE_FILE:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->profileFile));
	  else if (strcasecmp(firstword,"DATE:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->dte));
	  else if (strcasecmp(firstword,"UNITS:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->units));
	  else if (strcasecmp(firstword,"ID:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->user));
	  else if (strcasecmp(firstword,"DM_CORRECTION:")==0)
	    sscanf(line,"%s %lf",dummy,&(tmpl->dedispersed));
	  else if (strcasecmp(firstword,"NCHAN:")==0)
	    sscanf(line,"%s %d",dummy,&nchan);
	  else if (strcasecmp(firstword,"STOKES:")==0)
	    {
	      sscanf(line,"%s %s",dummy,stokesStr);
	      if (strcmp(stokesStr,"I")==0)
		nstokes=1;
	      else if (strcmp(stokesStr,"Q")==0 || strcmp(stokesStr,"U")==0 || strcmp(stokesStr,"V")==0)
		nstokes=4;
	    }
	}
      }
    }
  fclose(fin);
  // Do some checks
  if (nchan < 0){
    printf("Have not defined any channels. Unable to continue\n");
    exit(1);
  }
  // Allocate memory for these channels
  tmpl->nchan = nchan;
  if (tmpl->channelMemoryAllocated == 0){
    if (!(tmpl->channel = (channelStruct *)malloc(sizeof(channelStruct)*nchan))){
      printf("ERROR in allocated memory for channels\n");
      exit(1);
    }
    tmpl->channelMemoryAllocated = 1;
    tmpl->nChannelAllocated = nchan;
    for (i=0;i<nchan;i++)
      tmpl->channel[i].polMemoryAllocated = 0;
  }

  chan = -1;
  stokes = -1;
  comp=-1;
  // Now read the data
  if (!(fin = fopen(file,"r"))){
    printf("Unable to open file: >%s<\n",file);
    exit(1);
  }
  while (!feof(fin))
    {
      if (fgets(line,4096,fin)!=NULL){
	if (line[0] == '#') // Comment line
	  {
	    // Do nothing
	  }
	else {
	  sscanf(line,"%s",firstword);
	  if (strcasecmp(firstword,"STOKES:")==0){
	    sscanf(line,"%s %s",dummy,stokesStr);
	    if (strcmp(stokesStr,"I")==0)
	      stokes=0;
	    else if (strcmp(stokesStr,"Q")==0)
	      stokes = 1;
	    else if (strcmp(stokesStr,"U")==0)
	      stokes = 2;
	    else if (strcmp(stokesStr,"V")==0)
	      stokes = 3;
	  } 
	  else if (strcasecmp(firstword,"FREQUENCY_RANGE:")==0)
	    {
	      if (stokes==0)
		chan++;
	      sscanf(line,"%s %lf %lf",dummy,&f1,&f2);
	      if (f1 < f2){
		tmpl->channel[chan].freqLow = f1;
		tmpl->channel[chan].freqHigh = f2;
	      } else {
		tmpl->channel[chan].freqLow = f2;
		tmpl->channel[chan].freqHigh = f1;
	      }
	      tmpl->channel[chan].nstokes = nstokes;
	      // Allocate memory
	      if (tmpl->channel[chan].polMemoryAllocated==0){
		if (!(tmpl->channel[chan].pol = (polStruct *)malloc(sizeof(polStruct)*nstokes))){
		  printf("Error in allocated memory for Stokes\n");
		  exit(1);
		}
		tmpl->channel[chan].polMemoryAllocated = 1;
		for (i=0;i<nstokes;i++)
		  tmpl->channel[chan].pol[i].compMemoryAllocated = 0;
	      }
	      //	    }
	      tmpl->channel[chan].nPolAllocated = nstokes;
	    }
	  else if (strcasecmp(firstword,"NCOMP:")==0)
	    {
	      int ncomp;
	      sscanf(line,"%s %d",dummy,&ncomp);
	      tmpl->channel[chan].pol[stokes].nComp = ncomp;
	      tmpl->channel[chan].pol[stokes].stokes = stokes;
	      if (tmpl->channel[chan].pol[stokes].compMemoryAllocated==0){
		if (!(tmpl->channel[chan].pol[stokes].comp = (component *)malloc(sizeof(component)*ncomp))){
		  printf("Error in allocated memory for components\n");
		  exit(1);
		}
		tmpl->channel[chan].pol[stokes].compMemoryAllocated = 1;
	      }
	      tmpl->channel[chan].pol[stokes].nCompAllocated = ncomp;
	      comp=0;
	    }
	  else // Look for component
	    {
	      char substr[4096];
	      strcpy(substr,firstword);
	      substr[4]='\0';
	      if (strcasecmp(substr,"COMP")==0)
		{
		  sscanf(line,"%s %lf %lf %lf",dummy,
			 &(tmpl->channel[chan].pol[stokes].comp[comp].height),
			 &(tmpl->channel[chan].pol[stokes].comp[comp].concentration),
			 &(tmpl->channel[chan].pol[stokes].comp[comp].centroid));
		  comp++;
		}
	    }
	}
      }
    }
  fclose(fin);
}


// Evaluate a single template component
double evaluateTemplateComponent(tmplStruct *tmpl,double phi,int chan,int stokes,int comp,double phiRot)
{
  double result=0;
  result = tmpl->channel[chan].pol[stokes].comp[comp].height *
    exp(tmpl->channel[chan].pol[stokes].comp[comp].concentration*
	(cos((phi - tmpl->channel[chan].pol[stokes].comp[comp].centroid - phiRot)*2*M_PI)-1));
  return result;
}

// Evaluate a given frequency channel and polarisation
double evaluateTemplateChannel(tmplStruct *tmpl,double phi,int chan,int stokes,double phiRot)
{
  double result=0;
  int k;
  for (k=0;k<tmpl->channel[chan].pol[stokes].nComp;k++)
    result += evaluateTemplateComponent(tmpl,phi,chan,stokes,k,phiRot);
  return result;
}

// Allocate specific amount of memory
void allocateMemoryTemplateDefault(tmplStruct *tmpl,int nchan,int npol,int ncomp)
{
  int i,j;
  tmpl->channel = (channelStruct *)malloc(sizeof(channelStruct)*nchan);
  tmpl->channelMemoryAllocated = 1;
  tmpl->nChannelAllocated = nchan;

  for (i=0;i<nchan;i++)
    {
      tmpl->channel[i].pol = (polStruct *)malloc(sizeof(polStruct)*npol);
      tmpl->channel[i].polMemoryAllocated = 1;
      tmpl->channel[i].nPolAllocated = npol;

      for (j=0;j<npol;j++)
	{
	  tmpl->channel[i].pol[j].comp = (component *)malloc(sizeof(component)*ncomp);
	  tmpl->channel[i].pol[j].compMemoryAllocated = 1;
	  tmpl->channel[i].pol[j].nCompAllocated = ncomp;
	}
    }
}


void saveTemplate(char *fname,tmplStruct *tmpl)
{
  FILE *fout;
  char version[128] = "1.0";
  int i,j,k;
  double centre;

  if (!(fout = fopen(fname,"w")))
    {
      printf("Unable to open file %s\n",fname);
      return;
    }
  fprintf(fout,"TEMPLATE_VERSION: %s\n",version);
  fprintf(fout,"PROFILE_FILE: %s\n",tmpl->profileFile);
  fprintf(fout,"SOURCE: %s\n",tmpl->source);
  fprintf(fout,"DATE: %s\n",tmpl->dte);
  fprintf(fout,"ID: %s\n",tmpl->user);
  fprintf(fout,"UNITS: %s\n",tmpl->units);
  fprintf(fout,"NCHAN: %d\n",tmpl->nchan);
  for (i=0;i<tmpl->nchan;i++)
    {
      fprintf(fout,"#\n");
      for (j=0;j<tmpl->channel[i].nstokes;j++)
	{
	  fprintf(fout,"#\n");
	  if (tmpl->channel[i].pol[j].stokes == 0)
	    fprintf(fout,"STOKES: I\n");
	  else if (tmpl->channel[i].pol[j].stokes == 1)
	    fprintf(fout,"STOKES: Q\n");
	  else if (tmpl->channel[i].pol[j].stokes == 2)
	    fprintf(fout,"STOKES: U\n");
	  else if (tmpl->channel[i].pol[j].stokes == 3)
	    fprintf(fout,"STOKES: V\n");
	  fprintf(fout,"FREQUENCY_RANGE: %f %f\n",tmpl->channel[i].freqLow,tmpl->channel[i].freqHigh);
	  fprintf(fout,"NCOMP: %d\n",tmpl->channel[i].pol[j].nComp);
	  for (k=0;k<tmpl->channel[i].pol[j].nComp;k++){
	    centre = tmpl->channel[i].pol[j].comp[k].centroid;
	    if (centre > 1) centre-=1;
	    if (centre < 0) centre+=1;
	    fprintf(fout,"COMP%d: %g %g %g\n",k+1,tmpl->channel[i].pol[j].comp[k].height,
		    tmpl->channel[i].pol[j].comp[k].concentration,
		    centre);
	  }
	}
    }
  fclose(fout);
}

