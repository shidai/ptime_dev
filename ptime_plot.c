// To do:
// 1. Plot invariant interval
// 2. Plot percentage polarised
// 3. Plot percentage linearly polarised
// 4. Plot percentage circularly polarised
// 7. Overlay other profiles
// 8. Plot multiple frequency channels
// 9. Make summary plot with all templates
// 10. Have titles on plot
// 11. Toggle displaying legend (and position of the legend)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ptime.h"
#include <cpgplot.h>

// Program to plot templates

void plot(tmplStruct *tmpl,int nbin);
void findMinMax(int n,float *y,float *miny,float *maxy);
void findMinMax4(int n,float *y1,float *y2,float *y3,float *y4,float *miny,float *maxy);

int main(int argc,char *argv[])
{
  tmplStruct tmpl;
  int i;
  int nbin = 1024;
  char fname[128];

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-nbin")==0)
	sscanf(argv[++i],"%d",&nbin);
    }

  initialiseTemplate(&tmpl);
  printf("Reading template\n");
  readTemplate(fname,&tmpl); 
  printf("Complete reading template\n");
  plot(&tmpl,nbin);
  // Must deallocate the template memory
}

void plot(tmplStruct *tmpl,int nbin)
{
  float *fx;
  float *pol1,*pol2,*pol3,*pol4,*lin,*tpol,*pa,*pa_x;
  int pa_n=0;
  float *fyc;
  int n=0,n0=0;
  float minx,maxx,miny,maxy;
  double phiRot=0.5;
  char yaxis[128];
  float mx,my;
  char key;
  double phi;
  int plotComponents=0;
  int i,j;
  int channel=0;
  int pol=0;
  int nchan = tmpl->nchan;
  int npol = tmpl->channel[0].nstokes; 
  int plot=1;
  int zoom=0;
  double smallVal = 0.1;

  if (nchan == 1 && npol == 1)
    {
      plot = 1;
      plotComponents=1;
    }
  else
    {
      plot = 3;
    }
  fx = (float *)malloc(sizeof(float)*nbin);
  pol1 = (float *)malloc(sizeof(float)*nbin);
  pol2 = (float *)malloc(sizeof(float)*nbin);
  pol3 = (float *)malloc(sizeof(float)*nbin);
  pol4 = (float *)malloc(sizeof(float)*nbin);
  lin = (float *)malloc(sizeof(float)*nbin);
  tpol = (float *)malloc(sizeof(float)*nbin);
  pa = (float *)malloc(sizeof(float)*nbin);
  pa_x = (float *)malloc(sizeof(float)*nbin);

  fyc = (float *)malloc(sizeof(float)*nbin);



  sprintf(yaxis,"Amplitude (%s)",tmpl->units);

  cpgbeg(0,"/xs",1,1);
  cpgsch(1.4);
  cpgscf(2);
  cpgask(0);
    


  do {
    n=0;
    pa_n=0;
    for (phi=0;phi<1;phi+=1.0/(double)nbin)
      {
	fx[n] = (float)phi;
	pol1[n] = (float)evaluateTemplateChannel(tmpl,phi,channel,0,phiRot);
	if (npol==4)
	  {
	    pol2[n] = (float)evaluateTemplateChannel(tmpl,phi,channel,1,phiRot);
	    pol3[n] = (float)evaluateTemplateChannel(tmpl,phi,channel,2,phiRot);
	    pol4[n] = (float)evaluateTemplateChannel(tmpl,phi,channel,3,phiRot);
	    lin[n] = sqrt(pol2[n]*pol2[n] + pol3[n]*pol3[n]);
	    tpol[n] = sqrt(pol2[n]*pol2[n] + pol3[n]*pol3[n] + pol4[n]*pol4[n]);
	    if (fabs(pol3[n]) > smallVal && fabs(pol2[n]) > smallVal){
	      pa[pa_n] = atan2(pol2[n],pol3[n])*180.0/M_PI;
	      pa_x[pa_n] = fx[n];
	      pa_n++;
	    }
	  }
	n++;
      }
    if (zoom==0){
      minx = 0;
      maxx = 1;  


      if (plot==1 || plot == 4)
	{
	  float tmaxy,tminy;
	  findMinMax(n,pol1,&tminy,&tmaxy);
	  maxy = tmaxy + 0.1*(tmaxy-tminy);
	  miny = tminy - 0.1*(tmaxy-tminy);
	}
      else {
	float tmaxy,tminy;
	findMinMax4(n,pol1,pol2,pol3,pol4,&tminy,&tmaxy);
	  maxy = tmaxy + 0.1*(tmaxy-tminy);
	  miny = tminy - 0.1*(tmaxy-tminy);
      }
    }
    if (plot==1)
      {
	cpgenv(minx,maxx,miny,maxy,0,1);
	cpglab("Phase",yaxis,"");
	cpgline(n,fx,pol1);
	
	if (plotComponents==1)
	  {
	    for (i=0;i<tmpl->channel[0].pol[0].nComp;i++){
	      n0=0;
	      for (phi=0;phi<1;phi+=1.0/(double)nbin)
		{
		  fyc[n0] = (float)evaluateTemplateComponent(tmpl,phi,channel,pol,i,phiRot);
		  n0++;
		}
	      cpgsls(4);
	      cpgsci(i+2);
	      cpgline(n,fx,fyc);
	      cpgsls(1); cpgsci(1);
	    }
	  }
      } 
    else if (plot == 2) 
      {
	cpgenv(minx,maxx,miny,maxy,0,1);
	cpglab("Phase",yaxis,"");

	cpgsci(1); cpgline(n,fx,pol1);
	cpgsci(2); cpgline(n,fx,pol2);
	cpgsci(3); cpgline(n,fx,pol3);
	cpgsci(4); cpgline(n,fx,pol4);
	cpgsci(1);
      }
    else if (plot == 3) 
      {
	//	cpgenv(minx,maxx,miny,maxy,0,1);
	//	cpglab("Phase",yaxis,"");
	cpgeras();

	cpgsvp(0.1,0.9,0.55,0.9);
	cpgswin(minx,maxx,-90,90);
	cpgbox("BCTS",0.0,0,"BCTSN",0.0,0);
	cpglab("","PA","");
	cpgpt(pa_n,pa_x,pa,1);

	cpgsvp(0.1,0.9,0.15,0.55);
	cpgswin(minx,maxx,miny,maxy);
	cpgbox("BCTNS",0.0,0,"BCTSN",0.0,0);
	cpglab("Phase","Template","");
	cpgsci(1); cpgline(n,fx,pol1);
	cpgsci(2); cpgline(n,fx,lin);
	cpgsci(4); cpgline(n,fx,pol4);
	cpgsci(1);
      }
    else if (plot == 4) 
      {
	cpgenv(minx,maxx,miny,maxy,0,1);
	cpglab("Phase",yaxis,"");

	cpgsci(1); cpgline(n,fx,pol1);
	cpgsci(2); cpgline(n,fx,tpol);
	cpgsci(1);
      }
    else if (plot == 5) 
      {
	float tmaxy,tminy;
	cpgenv(minx,maxx,0,nchan,0,1);
	//	cpgbox("G",0,0,"G",0,0);
	cpglab("Phase","Frequency channel","");
	for (i=0;i<nchan;i++)
	  {
	    n0=0;
	    for (phi=0;phi<1;phi+=1.0/(double)nbin)
	      {
		fyc[n0] = (float)evaluateTemplateChannel(tmpl,phi,i,pol,phiRot);
		n0++;
	      }
	    findMinMax(n0,fyc,&tminy,&tmaxy);
	    for (j=0;j<n0;j++)
	      {
		fyc[j]=(fyc[j]+tminy)/(tmaxy-tminy) + i;
	      }
	    cpgsci(1); cpgline(n0,fx,fyc);
	  }

	cpgsci(1);
      }

    else if (plot == 6) 
      {

	float tfx[8192];
	double p1,p2,p3,p4;
	cpgenv(minx,maxx,-90,90,0,1);
	//	cpgbox("G",0,0,"G",0,0);
	cpglab("Phase","PA","");
	for (i=0;i<nchan;i++)
	  {
	    n0=0;
	    n=0;
	    for (phi=0;phi<1;phi+=1.0/(double)nbin)
	      {
		p1 = (float)evaluateTemplateChannel(tmpl,phi,i,0,phiRot);
		p2 = (float)evaluateTemplateChannel(tmpl,phi,i,1,phiRot);
		p3 = (float)evaluateTemplateChannel(tmpl,phi,i,2,phiRot);
		p4 = (float)evaluateTemplateChannel(tmpl,phi,i,3,phiRot);
		if (fabs(p3) > smallVal && fabs(p2) > smallVal){
		  fyc[n0] = atan2(p2,p3)*180.0/M_PI;
		  tfx[n0] = fx[n];
		  n0++;
		}
		n++;
	      }
	    cpgsci(i+1);
	    cpgpt(n0,tfx,fyc,1);
	    cpgsci(1);
	  }

	cpgsci(1);
      }
	


    cpgcurs(&mx,&my,&key);
    if (key=='r'){
      printf("Current rotation is %g (phase). Please enter new rotation ",phiRot);
      scanf("%lf",&phiRot);
    }
    else if (key=='s'){
      printf("The small value used for determining PA is set to: %g.  What value would you like? ",smallVal);
      scanf("%lf",&smallVal);
    }
    else if (key=='1') plot = 1;
    else if (key=='2') plot = 2;
    else if (key=='3') plot = 3;
    else if (key=='4') plot = 4;
    else if (key=='5') plot = 5;
    else if (key=='6') plot = 6;
    else if (key=='z')
      {
	float mx2,my2;
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	if (mx2 < mx) minx = mx2;
	else minx = mx;
	if (mx2 > mx) maxx = mx2;
	else maxx = mx;

	if (my2 < my) miny = my2;
	else miny = my;
	if (my2 > my) maxy = my2;
	else maxy = my;
	zoom=1;
      }
    else if (key=='u')
      zoom=0;
  } while (key!='q');
    cpgend();


  // Must deallocate the memory
  free(fx);
  free(pol1);
  free(pol2);
  free(pol3);
  free(pol4);
  free(lin);
  free(tpol);
  free(pa);
  free(pa_x);
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

void findMinMax4(int n,float *y1,float *y2,float *y3,float *y4,float *miny,float *maxy)
{
  int i;
  float temp[4*n];

  for (i=0;i<n;i++)
    {
      temp[i] = y1[i];
      temp[n+i] = y2[i];
      temp[2*n+i] = y3[i];
      temp[3*n+i] = y4[i];
    }
  findMinMax(4*n,temp,miny,maxy);
}
