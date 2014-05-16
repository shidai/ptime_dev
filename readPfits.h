// Header information for pfits software

#include "fitsio.h"

typedef struct pheader {
  int nhead;
  //
  // All header information
  char **keyname;
  char **val;
  char **comment;
  //
  // Most useful header information
  int nchan; // Number of channels
  int nbits; // Number of bits
  int nsamp; // Number of samples
  int nsub;  // Number of subintegrations
  int nsblk; // Number of samples per subintegration
  int nbin;  // Number of bins
  int npol;  // Number of polarisations
  float tsamp; // Sample time
  float freq;  // Centre frequency for observation
  float bw;    // Observation bandwidth
  float zeroOff; // Zero offsets
  int   imjd;  // Integer start time (day)
  float smjd;  // Fraction start time (sec)
  float stt_offs; // Start time offset (seconds)
  float chanbw;   // Channel bandwidth (MHz)
  float dm;       // Pulsar's dispersion measure (cm^-3 pc)
  float period;   // Pulsar's period (s)   
  char obsMode[128]; // PSR, CAL, SEARCH
  char source[128];
} pheader;

void closeFitsFile(fitsfile *fp);
fitsfile * openFitsFile(char *fname);
void loadPrimaryHeader(fitsfile *fp,pheader *phead);
