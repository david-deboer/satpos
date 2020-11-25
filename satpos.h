#ifndef _MODELS_H
#define _MODELS_H

/* If MAIN is defined, this header will generate definitions.  Otherwise,
   only declarations will be generated. */

#ifdef MAIN   /* generate DEFINITIONS */
	#define CLASS
	#define EQUALS(N) = N
struct satTime {
	int Year;
	int Month;
	int Day;
	int Hour;
	int Minute;
	double Second;
	double time;
};
struct observer {
	char Name[12];
	char Code[8];
	double lat;
	double lng;
	double alt;
	double h;
	double range;
	double Horizon;
	int timeZone;
	int daylightSavings;
	struct satTime tstart;
	struct satTime tstop;
	struct satTime tnow;
	double tstep;
	double Xlam;
	double Ylam;
	double Zlam;
};
#else       /* generate DECLARATIONS */
	#define CLASS extern
	#define EQUALS(N)
#endif

enum {NOPRINT, YESPRINT, PRINTLABEL};

#include "astroliblite.h"
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SQ(X)  ( (X)*(X) )
#undef VERBOSE

int readObserver(struct observer *obs);
int readLocation(struct observer *obs);
int invxyz(double x, double y, double z, double *lng, double *lat, double *h);
int otherTerms(struct observer obs, double jd, double *ro, struct observer *subsat,
			   double *Az, double *El, double *RA, double *Dec);
int writeFootprint(double lngs, double lats, double rsat);
void addTime(struct satTime *base_t, struct satTime diff_t);
void resetTime(struct satTime *t);

#endif
