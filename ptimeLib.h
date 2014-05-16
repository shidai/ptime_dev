#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef __ptimeLib_h
#define __ptimeLib_h

typedef struct component {
  double concentration; // Concentration of each component
  double height;        // Height of each component
  double centroid;      // Centroid of each component
} component;

typedef struct polStruct {
  int stokes;            // 1 = I, 2 = Q, 3 = U, 4 = V
  int nComp;             // Number of components for each channel for each Stokes
  component *comp;
  int compMemoryAllocated;
  int nCompAllocated;
} polStruct;

typedef struct channelStruct {
  int nstokes;             // 1 = I, 4 = I,Q,U,V
  double freqLow;
  double freqHigh;
  polStruct *pol; 
  int polMemoryAllocated;
  int nPolAllocated;
} channelStruct;

typedef struct tmplStruct {
  // Common to all templates
  char dte[1024];      // Date template made
  char user[1024];   // Person who made the template
  float templateVersion; // Version of template header
  char source[1024]; // Source name
  char profileFile[1024]; // Profile file name
  char units[1024];   // Unit definition
  double dedispersed; // = 0 by default

  int nchan; // Number of channels
  channelStruct *channel; // Channels
  int channelMemoryAllocated;
  int nChannelAllocated;
} tmplStruct;



void initialiseTemplate(tmplStruct *tmpl);
void readTemplate(char *file,tmplStruct *tmpl);
double evaluateTemplateComponent(tmplStruct *tmpl,double phi,int chan,int stokes,int comp,double phiRot);
double evaluateTemplateChannel(tmplStruct *tmpl,double phi,int chan,int stokes,double phiRot);
void allocateMemoryTemplateDefault(tmplStruct *tmpl,int nchan,int npol,int ncomp);
void saveTemplate(char *fname,tmplStruct *tmpl);

#endif
