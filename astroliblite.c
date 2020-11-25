/*-------------------------Astro Library---------------------------*/
/* Astronomy library program (with several utilities)              */
/* See AstroLib.h for listing of functions                         */
/* ...                                                             */
/* 04APR1999 David R. DeBoer                                       */
/* 26NOV1999 DRD fixed ConvertToDDMMSS for dropping seconds        */
/* 28FEB2001 DRD made this AstroLibLite to try and debug           */
/*-----------------------------------------------------------------*/

#include "astroliblite.h"

/*---------------------------EQ2HOR------------------------------------*/
int Eq2Hor(double xt,double yt,double *pt,double *qt,double eslat)
{
/* This subroutine converts Az/El to Hr Ang/dec and vice versa.
   converts xt,yt -> pt,qt where
            xt = Az or HA, yt = El or Dec
            pt = HA or Az, qt = Dec or El
   All values are in degrees.
   eslat is the latitude of the earth station.
   Taken from Peter Duffett-Smith, "Astronomy on your Personal
   Computer". Extensive error-checking is done.*/

      double x, y, phi, sphi, cphi, sy, cy, sx, cx,
           sq, q, cq, a, cp, cptemp, p;

      phi = eslat*PI/180.0;
      x   = xt *PI/180.0;
      y   = yt *PI/180.0;

      sphi = sin(phi);
      cphi = cos(phi);
      sy = sin(y);
      cy = cos(y);
      sx = sin(x);
      cx = cos(x);
      sq = (sy*sphi) + (cy*cx*cphi);
      if (fabs(sq)>1.0) sq = fabs(sq)/sq;
      if(sq==1.0)
            q  =  PI/2.0;
      else if(sq==-1.0)
            q  = -PI/2.0;
      else
            q  = atan2( sq, (Esqrt(1.0-sq*sq)+1.0E-20) );
      cq = cos(q);
      a  = cphi*cq;
      if (a<1.0E-20)
            a = 1.0E-20;
      cp = (sy-(sphi*sq))/a;
      if (fabs(cp)>1.0) cp = fabs(cp)/cp;
      if (cp==1.0)
            cptemp = 0.0;
      else if (cp==-1.0)
            cptemp = -PI;
      else
            cptemp  = atan2( cp, (Esqrt(1.0-cp*cp)+1.0E-20) );
      p  =  PI/2.0 - cptemp;
      if (sx>0.0)
            p = TWOPI - p;

      *pt = p*180.0/PI;
      *qt = q*180.0/PI;

      return(1);
}



/*------------------------------------------------------------------------------------*/
double Time2LST(double lng, double tz, double time, int day, int mon, int year)
{
/*
; NAME:
;       Time2LST
; PURPOSE:
;       To convert from Local Civil Time to Local Mean Sidereal Time.
; INPUTS:
;       lng   The longitude in degrees of the place for which the local 
;               sidereal time is desired, scalar.   The Greenwich mean sidereal
;               time (GMST) can be found by setting Lng = 0.
;       tz    The time zone of the site in hours.  Use this to easily account 
;               for Daylight Savings time (e.g. 4=EDT, 5 = EST/CDT), scalar
;       time  The time of day of the specified date in decimal hours. 
;       day   The day of the month (1-31)
;       mon   The month, in numerical format (1-12)
;       year  The year (e.g. 1987)
;
; OUTPUTS:
;       Lst   The Local Sideral Time for the date/time specified in hours.
;
; RESTRICTIONS:
;       If specified, the date should be in numerical form.  The year should
;       appear as yyyy.
;
; PROCEDURE:
;       The Julian date of the day and time is question is used to determine
;       the number of days to have passed since 0 Jan 2000.  This is used
;       in conjunction with the GST of that date to extrapolate to the current
;       GST; this is then used to get the LST.    See Astronomical Algorithms
;       by Jean Meeus, p. 83 for the constants used.
;
; EXAMPLE:
;       Find the Greenwich mean sidereal time (GMST) on 1987 April 10, 19h21m UT
;
;       For GMST, we set lng=0, and for UT we set Tz = 0
;
;       IDL> Time2LST, lst, 0, 0,ten(19,21), 10, 4, 1987
;
;               ==> lst =  8.5825249 hours  (= 8h 34m 57.0896s)
;
; PROCEDURES USED:
;       Date2Julian - Convert from year, month, day, hour to julian date
;
; MODIFICATION HISTORY:
;       Adapted from the FORTRAN program GETSD by Michael R. Greason, STX, 
;               27 October 1988.
;       Use IAU 1984 constants Wayne Landsman, HSTX, April 1995, results 
;               differ by about 0.1 seconds  
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
*/
      double jd, c[4], jd2000, t0, t, theta, lst, gst;

/* compute the Julian date */
      if (tz != 0.0)   /* Account for time-zone if not UTC */
            __IncrementDateTime(&year,&mon,&day,&time,tz);
      jd=Date2Julian(year, mon, day, time);

/*                            Useful constants, see Meeus, p.84 */
      c[0] = 280.46061837;
      c[1] = 360.98564736629;
      c[2] = 0.000387933;
      c[3] = 38710000.0;
      jd2000 = 2451545.0;
      t0 = jd - jd2000;
      t = t0/36525.0;

/*                            Compute GST in seconds */
      theta = c[0] + (c[1] * t0) + t*t*(c[2] - t/ c[3] );

/*                            Compute GST in hours */
      lst = (theta+lng)/15.0;
      if (lst<0.0) gst = 24.0 + __Modulo(lst,24.0);
      lst=__Modulo(lst,24.0);

      return lst;
}


double Date2Julian(int year, int mon, int day, double HR)
{
/*
; NAME:
;       Date2Julian
; PURPOSE:
;       Converts Gregorian dates to Julian days   
;
; INPUTS:
;       year  = Year (integer)  
;       month = Month (integer 1-12)
;       day   = Day  (integer 1-31) 
;       HR    = Hours and fractions of hours of universal time (U.T.)
;               
; OUTPUT:
;       JULIAN = Julian date
;
; EXAMPLE:
;       To find the Julian Date at 1978 January 1, 0h (U.T.)
;
;       IDL> Date2Julian, 1978, 1, 1, 0., JULIAN
;
;       will give JULIAN = 2443509.5
; NOTES:
;       () JULDATE is an alternate procedure to perform the same function
;
; REVISON HISTORY:
;       Converted to IDL from Don Yeomans Comet Ephemeris Generator,
;       B. Pfarr, STX, 6/15/88
;       Converted to IDL V5.0   W. Landsman   September 1997
*/

      int L;
      long jd;
      double julian;

      L = (mon-14)/12;         /*In leap years, -1 for Jan, Feb, else 0*/

      jd = day - 32075L + 1461L*(year+4800L+L)/4L + 
             367L*(mon - 2L-L*12L)/12L - 3L*((year+4900L+L)/100L)/4L;
      julian = (double) jd + (HR/24.0) - 0.5;

      return julian;
}


/*--------------------------------------------------------------------------*/
/* Returns 1 if y is leap year in Gregorian calendar, 0 otherwise.  To
signify BC year use a negative number, e.g., -1 for 1 BC. */
int __LeapYear(int y)
{
      /* For leap year determination purposes, 1 BC = 0, 2 BC = -1, etc.
      Therefore, we must make an adjustment to y if it's BC. */
      if (y < 0)
            ++y;
      return !(y%4) && y%100  ||  !(y%400);
}

/*--------------------------------------------------------------------------*/
void __IncrementDateTime(int *Year, int *Month, int *Day, double *Time, double Step)
{
      *Time += Step;

      if (*Time>=24.0)
      {
            *Time-=24.0;
            __IncrementDay(Month,Day,Year);
      }
      else if (*Time<0.0)
      {
            *Time+=24.0;
            __DecrementDay(Month,Day,Year);
      }
}

/*--------------------------------------------------------------------------*/
void __IncrementDay(int *mn, int *dy, int *yr)
{
      int leap;

      ++ (*dy);
      if (*dy < 29) return;

      leap = __LeapYear(*yr);   /* From util.cpp */

      if ((*mn==1 || *mn==3 || *mn==5 || *mn==7 || *mn==8 || *mn==10) && *dy==32)
      {
            ++(*mn);
            *dy=1;
      }
      else if ((*mn==4 || *mn==6 || *mn==9 || *mn==11) && *dy==31)
      {
            ++(*mn);
            *dy=1;
      }
      else if (*mn==2 && *dy==29+leap)
      {
            *mn=3;
            *dy=1;
      }
      else if (*mn==12 && *dy==32)
      {
            *mn=1;
            *dy=1;
            ++(*yr);
      }

      return;
}

/*--------------------------------------------------------------------------*/
void __DecrementDay(int *mn, int *dy, int *yr)
{
      int leap;

      leap = __LeapYear(*yr);

      --(*dy);
      if (*dy > 1) return;

      if ( *mn==2 || *mn==4 || *mn==6 || *mn==8 || *mn==9 || *mn==11 )
      {
            --(*mn);
            *dy=31;
      }
      else if ((*mn==5 || *mn==7 || *mn==10 || *mn==12) && *dy==31)
      {
            --(*mn);
            *dy=30;
      }
      else if (*mn==1)
      {
            *mn=12;
            *dy=31;
            --(*yr);
      }
      else if (*mn==3)
      {
            *mn=2;
            *dy=28+leap;
      }

      return;
}

/*-----------------------------------------------------------------------*/
double __Modulo(double x, double y)
{
      double ddiv, m;

      ddiv = x/y;

      ddiv = floor(ddiv);

      m = x - ddiv*y;

      return m;
}

/* Convert decimal to Degrees, Minutes and Seconds */
/* If LeadingSign is true include it, else don't   */

void ConvertToDDMMSS(double X, char *DDMMSS, int LeadingSign)
{
      int iD, iM, iS, digits;
      double absX, absM, absS, sign, DD, MM, SS;

	  digits = fabs(LeadingSign);  //partially switching over to the "digits" concept...
	  if (digits < 2) digits = 2;
	  if (LeadingSign < 0)
		  LeadingSign = 0;

      absX = fabs(X);
      sign = absX/X;
      DD   = floor(absX);
      iD   = (int) DD;
      absM = (absX - DD)*60.0;
      MM   = floor(absM);
      iM   = (int) MM;
      absS = (absM - MM)*60.0;
      SS   = absS;
      iS   = (int) SS;

      if ( (absX - DD + MM/60.0 + SS/3600.0) > 0.000278 )
            ++iS;
      if (iS == 60) {iS=0; ++iM;}
      if (iM == 60) {iM=0; ++iD;}

      if (LeadingSign)
      {
            if (sign < 0.0)
                  sprintf(DDMMSS,"-%03d %02d %02.2lf",iD,iM,SS);
            else
                  sprintf(DDMMSS,"+%03d %02d %02.2lf",iD,iM,SS);
      }
      else
	  {
		if (digits==2)
			sprintf(DDMMSS,"%c%02d %02d %02.2lf",(sign<0.0)?'-':' ',iD,iM,SS);
		else if (digits==3)
			sprintf(DDMMSS,"%c%03d %02d %02.2lf",(sign<0.0)?'-':' ',iD,iM,SS);
		else
			sprintf(DDMMSS,"%c%d %02d %02.2lf",(sign<0.0)?'-':' ',iD,iM,SS);
	  }
}

/*---------------------------SUBAZEL----------------------------------*/
double subazel(double lngs, double Lats, double *Az,double *El,double rsat,struct GeoPosition *pos)
{
      double cg,sg,cosEl,Late,lnge,alpha, r, d, sa,rd,lsiv, mle, mls, alt;

/* Convert subsatellite point to Az/El
   See Pratt and Bostian p. 24 */

      /*  Calculate El */
      alt   = pos->Altitude/1000.0;
      Late  = pos->Latitude;
      Lats *= PI/180.0;
      lnge  = -1.0*pos->Longitude;
      lngs *= -1.0*PI/180.0;
      cg = cos(Late)*cos(Lats)*cos(lngs-lnge)+sin(Late)*sin(Lats);
      if(cg>1.0)
            cg = 1.0;
      else if(cg<-1.0)
            cg = -1.0;
      sg    = Esqrt( 1.0 - cg*cg );
      r     = (REARTH+alt)/rsat;
      d     = rsat* Esqrt(1.0 + r*r - 2.0*r*cg);
      rd    = (REARTH+alt)/d;

      if( (cg>=rd) && (cg<=1.0) )         /*  0 < EL < 90 */
      {
            cosEl = rsat*sg/d;
            if(cosEl>1.0)
                  *El = 0.0;
            else if(cosEl<-1.0)
                  *El = -PI;
            else
                  *El = acos( cosEl );
      }
      else                                /*  -90 < EL < 0 */
      {
            cosEl = -rsat*sg/d;
            if(cosEl>1.0)
                  *El = 0.0;
            else if(cosEl<-1.0)
                  *El = -PI;
            else
                  *El = acos( cosEl ) - PI;
      }

      *El    = 180.0*(*El)/PI;

      /*  Calculate Az */
      if(sg==0.0)
            sa = 1.0;
      else
            sa = sin(fabs(lnge-lngs))/sg;

      if(fabs(sa)>1.0)
            alpha = PI/2.0;
      else
            alpha = fabs( asin(sa) );

      if(     (lngs>lnge) && (Lats<=Late) )         /* s/c SW of g/s */
            *Az = 180.0*(1.0 + alpha/PI);
      else if( (lngs<=lnge) && (Lats<=Late) )       /* s/c SE of g/s */
            *Az = 180.0*(1.0 - alpha/PI);
      else if( (lngs>lnge) && (Lats>Late) )         /* s/c NW of g/s */
            *Az = 180.0*(2.0 - alpha/PI);
      else if( (lngs<=lnge) && (Lats>Late) )        /* s/c NE of g/s */
            *Az = 180.0*alpha/PI;
      else
            printf("Whoops\n");

      /* Kluge to allow for date-line straddling */
      if(lngs != 0.0)
            lsiv = lnge/lngs;
      else
            lsiv = 1.0;
      mle  = 180.0*fabs(lnge)/PI;
      mls  = 180.0*fabs(lngs)/PI;
      if ( (lsiv<=0.0) && (mle>90.0) && (mls>90.0) )
            *Az = 360.0-*Az;

      return(d);
}

/*---------------------------Parallax----------------------------------------------*/
/* Convert from geocentric HA/Dec to "located" HA/Dec                              */
/* Adapted from Peter Duffett-Smith "Practical Astronomy with your Calculator"     */
/*   Args:  r    -- Distance to object in kilometers from earth's center (double)  */
/*          gHA  -- Geocentric Hour Angle in degrees (double)                      */
/*          gDec -- Geocentric Declination in degrees (double)                     */
/*          lHA  -- Located Hour Angle in degrees (double *)                       */
/*          lDec -- Located Declination in degrees (double *)                      */
/*          pos  -- Geographic location structure (struct GeoPosition *)           */
/*   Return: Angular separation of geocentric to located HA/Dec in degrees (double)*/
/*                                                                                 */
/* Written 3/25/99  David R. DeBoer                                                */
/*---------------------------------------------------------------------------------*/

double Parallax(double r, double gHA, double gDec, double *lHA, double *lDec, struct GeoPosition *Pos)
{
      double hp, u, rsp, rcp, D, dist, lat;

      gHA  *= PI/180.0;
      gDec *= PI/180.0;
      lat   = Pos->Latitude;
      r    /=REARTH;

      hp = ( (Pos->Altitude) /1000.0)/REARTH;
      u  = atan(0.996647*tan(lat));

      rsp = 0.996647*sin(u) + hp*sin(lat);
      rcp = cos(u) + hp*cos(lat);

      D = atan2( (rcp*sin(gHA)), (r*cos(gDec) - rcp*cos(gHA)) );
      *lHA = gHA + D;
      *lDec= atan( cos(*lHA)*(r*sin(gDec) - rsp)/(r*cos(gDec)*cos(gHA) - rcp) );

      *lHA  *= 180.0/PI;
      *lDec *= 180.0/PI;
      gHA   *= 180.0/PI;
      gDec  *= 180.0/PI;
      dist = AngularDistEq(gHA,gDec,*lHA,*lDec);
      return dist;
}

/*---------------------------AntiParallax------------------------------------------*/
/* Convert from "located" HA/Dec to geocentric HA/Dec                              */
/* See Smart, Textbook on Spherical Astronomy and Duffett-Smith                    */
/*     pp 199-206                                   p 36                           */
/*   Args:  rp   -- Distance to object in kilometers from surface (double)         */
/*          lHA  -- Located Hour Angle in degrees (double *)                       */
/*          lDec -- Located Declination in degrees (double *)                      */
/*          gHA  -- Geocentric Hour Angle in degrees (double)                      */
/*          gDec -- Geocentric Declination in degrees (double)                     */
/*          pos  -- Geographic location structure (struct GeoPosition *)           */
/*   Return: Angular separation of geocentric to located HA/Dec in degrees (double)*/
/*                                                                                 */
/* Written 2/27/01  David R. DeBoer                                                */
/*---------------------------------------------------------------------------------*/

double AntiParallax(double rp, double lHA, double lDec, double *gHA, double *gDec, struct GeoPosition *Pos)
{
      double hp, u, rsp, rcp, D, dist, lat;

      lHA  *= PI/180.0;
      lDec *= PI/180.0;
      lat   = Pos->Latitude;
      rp   /=REARTH;

      hp = ( (Pos->Altitude) /1000.0)/REARTH;
      u  = atan(0.996647*tan(lat));

      rsp = 0.996647*sin(u) + hp*sin(lat);
      rcp = cos(u) + hp*cos(lat);

      D = atan2( (rcp*sin(lHA)), (rp*cos(lDec) + rcp*cos(lHA)) );
      *gHA = lHA - D;
      *gDec= atan( cos(*gHA)*(rp*sin(lDec) + rsp)/(rp*cos(lDec)*cos(lHA) + rcp) );

      *gHA  *= 180.0/PI;
      *gDec *= 180.0/PI;
      lHA   *= 180.0/PI;
      lDec  *= 180.0/PI;
      dist = AngularDistEq(lHA,lDec,*gHA,*gDec);
      return dist;
}

double AngularDistEq(double RA1, double Dec1, double RA2, double Dec2)
{
      double cd, d;

      RA1 *= PI/180.0;
      Dec1*= PI/180.0;
      RA2 *= PI/180.0;
      Dec2*= PI/180.0;

      cd = sin(Dec1)*sin(Dec2) + cos(Dec1)*cos(Dec2)*cos(RA1-RA2);

      d = acos(cd);

      return 180.0*d/PI;
}

int MonthNumber(char *mon)
{
      char m[4];

      m[0]=tolower(mon[0]);
      m[1]=tolower(mon[1]);
      m[2]=tolower(mon[2]);

      if (!strncmp(m,"jan",3)) return 1;
      if (!strncmp(m,"feb",3)) return 2;
      if (!strncmp(m,"mar",3)) return 3;
      if (!strncmp(m,"apr",3)) return 4;
      if (!strncmp(m,"may",3)) return 5;
      if (!strncmp(m,"jun",3)) return 6;
      if (!strncmp(m,"jul",3)) return 7;
      if (!strncmp(m,"aug",3)) return 8;
      if (!strncmp(m,"sep",3)) return 9;
      if (!strncmp(m,"oct",3)) return 10;
      if (!strncmp(m,"nov",3)) return 11;
      if (!strncmp(m,"dec",3)) return 12;
      return 0;
}

double fmod2p(double x)
/* Reduces x to range 0 - 2pi */
{
      x /= TWOPI;
      x = (x - dint(x)) * TWOPI;
      if (x < 0.)
            return x + TWOPI;
      else
            return x;
}

double dint(double x)
{
      if (x >= 0.)
            return floor(x);
      else
            return ceil(x);
}

/* Calculate the ECI coordinates given the geodetic coordinates and LST */
/* Taken from www.celestrak.com, Orbital Elements III                   */
void LatLng2ECI(double LST, struct GeoPosition *pos, double *eci)
{
      double f, C, S, r;

      f = 1.0/298.26;

      LST = (LST*15.0)*PI/180.0; /* Convert to radians */

      r = REARTH + pos->Altitude/1000.0;

      C = r/Esqrt( 1.0 + f*(f-2.0)*SQ(sin(pos->Latitude)) );

      S = SQ(1.0-f)*C;

      eci[0] = C*cos(pos->Latitude)*cos(LST);
      eci[1] = C*cos(pos->Latitude)*sin(LST);
      eci[2] = S*sin(pos->Latitude);

      return;
}

/*Return length of line*/
int getlineDDB(FILE *fp, char *s, int lim)
{
      char c;
      int i;

      for (i=0;i<lim-2 && (c=getc(fp))!=EOF && c!='\n';++i)
          s[i] = c;

      s[i] = '\0';
      return (i);
}

double Esqrt(double x)
{
	if (x<0.0)
   {
   	printf("SQRT Error (Esqrt):  0.0 returned.\n");
      return 0.0;
   }
   else
   	return sqrt(x);
}
