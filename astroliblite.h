#ifndef _ASTROLIB_H
#define _ASTROLIB_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifndef PI
  #define PI (3.14159654)
#endif
#ifdef TORAD
  #undef TORAD
#endif
#ifdef TWOPI
  #undef TWOPI
#endif
#ifdef FOURPI
  #undef FOURPI
#endif
#ifdef REARTH
  #undef REARTH
#endif
#define TORAD (0.017453292)
#define TWOPI (6.28318530718)
#define FOURPI (12.56637061)
#define REARTH (6378.135)  /* In kilometers */
#define REARTH_POLAR (6356.750052)
#define CENTURY(X) ( ( (X)>57 )?(1900.0):(2000.0) )
#define SQ(X)  ( (X)*(X) )
#define FABS(X) ( sqrt( SQ(X) ) )

struct GeoPosition {
      double Longitude;   /* In degrees */
      double Latitude;    /* In degrees */
      double Altitude;    /* above sea-level in meters */
      double HrsToUTC;    /* local+HrsToUTC = UTC : this one can change with DST etc */
	  double StdHrsToUTC; /* local standard time + StdHrsToUTC + DST_flag = UTC */
      char   Name[25];
	  char   longName[100];
      double LocalHorizon;/* elevation of local horizon in Deg  */
};

double Time2LST(double lng, double tz, double time, int day, int mon, int year);
double Date2Julian(int year, int mon, int day, double HR);
double parallacticAngle(double HA, double dec, double lat);
int Eq2Hor(double xt,double yt,double *pt,double *qt,double eslat);
double Parallax(double r, double gHA, double gDec, double *lHA, double *lDec, struct GeoPosition *Pos);
double AntiParallax(double rp, double lHA, double lDec, double *gHA, double *gDec, struct GeoPosition *Pos);
double AngularDistEq(double RA1, double Dec1, double RA2, double Dec2);
double subazel(double lngs,double Lats,double *Az,double *El,double *range, double rsat, struct GeoPosition *pos);
void   __IncrementDateTime(int *Year, int *Month, int *Day, double *Time, double Step);
int    __LeapYear(int y);
void   __IncrementDay(int *mn, int *dy, int *yr);
void   __DecrementDay(int *mn, int *dy, int *yr);
double __Modulo(double x, double y);
void   ConvertToDDMMSS(double X, char *DDMMSS, int LeadingSign);
int    MonthNumber(char *mon);
double fmod2p(double x);
double dint(double x);
void LatLng2ECI(double LST, struct GeoPosition *pos, double *eci);
void ECI2subsat(double LST, double x, double y, double z, double *lat, double *lng);
int getline(FILE *fp, char *s, int lim);
double Esqrt(double x);

#endif
