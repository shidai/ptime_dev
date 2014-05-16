// To do:
// 1. Must deallocate memory

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ptime.h"
#include <cpgplot.h>

// Program to modify templates

void writePaasModelFile(tmplStruct *tmpl,char *file);

int main(int argc,char *argv[])
{
  tmplStruct tmpl;
  int i;
  int nbin = 1024;
  char fname[128];
  char paasModelFile[128];
  int paasModel=0;

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-paasModel")==0)
	{
	  paasModel=1;
	  strcpy(paasModelFile,argv[++i]);
	}
    }
  
  initialiseTemplate(&tmpl);
  printf("Reading template\n");
  readTemplate(fname,&tmpl); 
  printf("Complete reading template\n");
  
  if (paasModel==1){
    writePaasModelFile(&tmpl,paasModelFile);
  }
  
  // Must deallocate the template memory
}

void writePaasModelFile(tmplStruct *tmpl,char *file)
{
  FILE *fout;
  if (!(fout = fopen(file,"w")))
    {
      printf("Unable to open output file >%s<\n",file);
    }
  else
    {
      int i;
      for (i=0;i<tmpl->channel[0].pol[0].nComp;i++)
	fprintf(fout,"%.5f %.5f %.5f\n",
		tmpl->channel[0].pol[0].comp[i].centroid,
		tmpl->channel[0].pol[0].comp[i].concentration,
		tmpl->channel[0].pol[0].comp[i].height);
      printf("Completed writing to: %s\n",file);
      fclose(fout);
    }
}
