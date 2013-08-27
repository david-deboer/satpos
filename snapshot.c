/*****************************main.c***************************/
//copy of satpos, but do a snapshot of satellites
#define MAIN
#include "satpos.h"
#undef MAIN

#include "sgp4ext.h"
#include "sgp4unit.h"
#include "sgp4io.h"

/**************************************************************/
/* All internal time is UTC.  Some I/O in local */

#undef VERBOSE
int main(int argc, char *argv[])
{
	int argSC=0, argTLE=0, useSC, i, j, k, n, valid;
	long mcnt, nup;
	int transit, printOutTransit, satsnap, finished;
	char TLEFile[25], SCName[40], SCSearch[10], line[50];
	double RA, HA, Dec, Az, El, dEl, dAz, oldel, oldaz, olddel;
	double DOY, olddoy, jday1, jnow, timeconv[3], lngconv[3];
	double SCVector[15];
	double w, cw, cDcH, cDsH, sD, period;
	FILE *fp, *fp2;
	struct observer obs, subsat, now;
	/////////////////This is cut-and-pasted in from testcpp.cpp///////////////
	char str[2];
	double ro[3];
	double vo[3];
	char typeinput, typerun, opsmode;
	gravconsttype  whichconst;
	int whichcon;
	FILE *outfile[4];
	// ----------------------------  locals  -------------------------------
	double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;
	double sec,  jd, rad, tsince, startmfe, stopmfe, deltamin, jdstart, jdstop;
	double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;
	int  year, mon, day, hr, min;
	char longstr1[130];
	typedef char str3[4];
	str3 monstr[13];
	char longstr2[130];
	elsetrec satrec;
	
	printOutTransit = PRINTLABEL;
	rad = 180.0 / pi;
	// ------------------------  implementation   --------------------------
	strcpy(monstr[1], "Jan");
	strcpy(monstr[2], "Feb");
	strcpy(monstr[3], "Mar");
	strcpy(monstr[4], "Apr");
	strcpy(monstr[5], "May");
	strcpy(monstr[6], "Jun");
	strcpy(monstr[7], "Jul");
	strcpy(monstr[8], "Aug");
	strcpy(monstr[9], "Sep");
	strcpy(monstr[10], "Oct");
	strcpy(monstr[11], "Nov");
	strcpy(monstr[12], "Dec");
	
	printf("%s\n",SGP4Version );
    printf("stk.v.4.3 \n"); // must use 4.3...
    printf("CoordinateSystem:  TEME-to-PEF\n\n");
    
    system("ls -1 tle/*.tle > tlefiles.dat");
	readObserver(&obs);
	jday( obs.tstart.Year,1,0,0,0,0.0, jday1 );
	jday( obs.tnow.Year,obs.tnow.Month,obs.tnow.Day,obs.tnow.Hour,obs.tnow.Minute,obs.tnow.Second, jnow );
	opsmode = 'a';
	typeinput = 'e';
	typerun = 'c';  // use 'c' (see testcpp.cpp) to just read in tle's.  start/stop gets set below (so overwrites the 'c' settings in sgp4io)
	whichcon = 72;
	if (whichcon == 721) whichconst = wgs72old;
	if (whichcon == 72) whichconst = wgs72;
	if (whichcon == 84) whichconst = wgs84;
	getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

	
	// ---------------- setup files for operation ------------------
	outfile[0] = fopen("snapshotXYZ.out", "w");
    outfile[1] = fopen("snapshotLLH.out", "w");
    outfile[2] = fopen("snapshotAERD.out","w");
    fp = fopen("tlefiles.dat","r");
    satsnap = 1;
    while (satsnap)
    {
        getline(fp,TLEFile,55);
        if (strlen(TLEFile)<4)
        {
            satsnap = 0;
        }
        else
        {
            fp2 = fopen(TLEFile,"r");
            finished = 0;
            while(!finished)
            {
                getline(fp2,SCName,50);
                for (i=strlen(SCName)-1; i>0; --i)
                {
                    if (isalnum(SCName[i]) or SCName[i]==')')
                    {
                        SCName[i+1] = '\0';
                        break;
                    }
                }
                getline(fp2,longstr1,100);
                n=getline(fp2,longstr2,100);
                if (n < 10)
                {
                    finished = 1;
                }
                else
                {
                    twoline2rv( longstr1, longstr2, typerun, typeinput, opsmode, whichconst,
                               startmfe, stopmfe, deltamin, satrec );
                    period = 2.0*pi/satrec.no;
                    if (period < 1400.0 || period > 1500.0)
                        continue;
                    // This is pulled out of sgp4io.cpp to change the start/stop times
                    jday( obs.tstart.Year,obs.tstart.Month,obs.tstart.Day,obs.tstart.Hour,obs.tstart.Minute,obs.tstart.Second, jdstart );
                    jday( obs.tstart.Year,obs.tstop.Month,obs.tstop.Day,obs.tstop.Hour,obs.tstop.Minute,obs.tstop.Second, jdstop );
                    startmfe = (jdstart - satrec.jdsatepoch) * 1440.0;
                    stopmfe  = (jdstop - satrec.jdsatepoch) * 1440.0;
                    deltamin = obs.tstep;
	
                    //fprintf(outfile, "%ld xx\n", satrec.satnum);
                    // call the propagator to get the initial state vector value
                    sgp4 (whichconst, satrec,  0.0, ro,  vo);

                    jd = satrec.jdsatepoch;
                    invjday( jd, year,mon,day,hr,min, sec );
                    printf("ScenarioEpoch:  %3i %3s%5i%3i:%2i:%12.9f \n",day,monstr[mon],year,hr,min,sec );

                    // First do it for 'now'
                    tsince = (jnow - satrec.jdsatepoch)*1440.0;
                    sgp4 (whichconst, satrec,  tsince, ro,  vo);
                    otherTerms(obs,jnow,ro,&now,&Az,&El,&RA,&Dec);
                    printf("===> %s <===\n",SCName);
                    printf("\t(x,y,z,r):  %lf, %lf, %lf, %lf\n",ro[0],ro[1],ro[2],now.alt);
                    printf("\t(sub_lng, sub_lat, h, range):  %lf, %lf, %lf %lf\n",now.lng,now.lat,now.h,now.range);
                    printf("\t(Az, El, RA, Dec):  %lf, %lf, %lf, %lf\n",Az,El,RA/15.0,Dec);
                    fprintf(outfile[0],"%lf\t%lf\t%lf\t%lf\n",ro[0],ro[1],ro[2],now.alt);
                    fprintf(outfile[1],"%lf\t%lf\t%lf\t%lf\n",now.lng,now.lat,now.h,now.range);
                    fprintf(outfile[2],"%lf\t%lf\t%lf\t%lf\n",Az,El,RA/15.0,Dec);
                } // if satrec.error == 0
            }
        }
    }
	
	return 1;
}

/* Get observer and date data */
int readObserver(struct observer *obs)
{
	int commentLine, nconv, timeZone;
	char line[300], code[8], *local_s;
	time_t time_now;
	tm *tm_now;
	struct satTime dt;
	FILE *fp;
	
	///////////  Get system time
	time_now = time(NULL);
	local_s = ctime(&time_now);
	//printf("%s",local_s);
	tm_now = localtime(&time_now);
	timeZone = tm_now->tm_hour;
	printf("System:\n");
	printf("\tLocal:  %02d:%02d:%02d   %02d/%02d/%d\n",tm_now->tm_hour, tm_now->tm_min, tm_now->tm_sec, tm_now->tm_mon+1, tm_now->tm_mday, tm_now->tm_year+1900);
	tm_now = gmtime(&time_now);
	timeZone-= tm_now->tm_hour;
	if (timeZone > 12)
		timeZone = timeZone - 24;
	else if (timeZone < -12)
		timeZone = timeZone + 24;
	printf("\tUTC:    %02d:%02d:%02d   %02d/%02d/%d\n",tm_now->tm_hour, tm_now->tm_min, tm_now->tm_sec, tm_now->tm_mon+1, tm_now->tm_mday, tm_now->tm_year+1900);
	printf("\tTimeZone:  %d\n",timeZone);
	
	///////////  Read input file
	fp=fopen("satpos.cfg","r");
	if (fp==NULL)
	{
		printf("satpos.cfg not found.\n");
		exit(1);
	}
	//obs_to_use daylight_savings_time
	commentLine = 1;
	while (commentLine)
	{
		getline(fp,line,299);
		if (line[0] != '#')
			commentLine = 0;
	}
	sscanf(line,"%s %d",obs->Code,&(obs->daylightSavings));
	nconv = readLocation(obs);
	if (!nconv)
	{
		printf("Location not set\n.");
		exit(1);
	}
	
	// current Year Month Day hour min second (UTC)
	commentLine = 1;
	while (commentLine)
	{
		getline(fp,line,299);
		if (line[0] != '#')
			commentLine = 0;
	}
	nconv = sscanf(line,"%d %d %d %d %d %lf",&(obs->tnow.Year),&(obs->tnow.Month),&(obs->tnow.Day),&(obs->tnow.Hour),&(obs->tnow.Minute),&(obs->tnow.Second));
	if (nconv==6)
	{
		resetTime(&dt);
		dt.Hour = -1*(obs->timeZone + obs->daylightSavings);
		addTime(&(obs->tnow),dt);  /// Convert to UTC
	}
	else 
	{
		obs->timeZone = timeZone;
		obs->tnow.Second = (float) tm_now->tm_sec;
		obs->tnow.Minute = tm_now->tm_min;
		obs->tnow.Hour = tm_now->tm_hour;
		obs->tnow.Day = tm_now->tm_mday;
		obs->tnow.Month = tm_now->tm_mon+1;
		obs->tnow.Year = tm_now->tm_year+1900;
	}
	obs->tnow.time = obs->tnow.Hour/24.0 + obs->tnow.Minute/24.0/60.0 + obs->tnow.Second/24.0/60.0/60.0;

	// start Year Month Day hour min second
	commentLine = 1;
	while (commentLine)
	{
		getline(fp,line,299);
		if (line[0] != '#')
			commentLine = 0;
	}
	nconv = sscanf(line,"%d %d %d %d %d %lf",&(obs->tstart.Year),&(obs->tstart.Month),&(obs->tstart.Day),&(obs->tstart.Hour),&(obs->tstart.Minute),&(obs->tstart.Second));
	if (nconv==6)
	{
		resetTime(&dt);
		dt.Hour = -1*(obs->timeZone + obs->daylightSavings);
		addTime(&(obs->tstart),dt);  /// Convert to UTC
	}
	else if (nconv==4)
	{
		resetTime(&dt);
		dt.Day = -1*obs->tstart.Year;       // the correct position in sscanf above
		dt.Hour = -1*obs->tstart.Month;     //   "
		dt.Minute = -1*obs->tstart.Day;     //   "
		dt.Second = -1*obs->tstart.Hour;    //   "
		obs->tstart.Second = obs->tnow.Second;
		obs->tstart.Minute = obs->tnow.Minute;
		obs->tstart.Hour = obs->tnow.Hour;
		obs->tstart.Day = obs->tnow.Day;
		obs->tstart.Month = obs->tnow.Month;
		obs->tstart.Year = obs->tnow.Year;
		addTime(&(obs->tstart),dt);
	}
	else 
	{
		printf("Incorrect time start - need either (all int):\n\tabsolute(year month day hour minute second) or delta(day hour minute second)\n");
		exit(1);
	}
	obs->tstart.time = obs->tstart.Hour/24.0 + obs->tstart.Minute/24.0/60.0 + obs->tstart.Second/24.0/60.0/60.0;

	// stop Year Month Day hour min second
	commentLine = 1;
	while (commentLine)
	{
		getline(fp,line,299);
		if (line[0] != '#')
			commentLine = 0;
	}
	nconv = sscanf(line,"%d %d %d %d %d %lf",&(obs->tstop.Year),&(obs->tstop.Month),&(obs->tstop.Day),&(obs->tstop.Hour),&(obs->tstop.Minute),&(obs->tstop.Second));
	if (nconv==6)
	{
		resetTime(&dt);
		dt.Hour = -1*(obs->timeZone + obs->daylightSavings);
		addTime(&(obs->tstop),dt);  /// Convert to UTC
	}
	else
	{
		if (nconv==4)  // otherwise use the same one as above
		{
			resetTime(&dt);
			dt.Day = obs->tstop.Year;
			dt.Hour = obs->tstop.Month;
			dt.Minute = obs->tstop.Day;
			dt.Second = obs->tstop.Hour;
		}
		else 
		{
			dt.Day*=-1;
			dt.Hour*=-1;
			dt.Minute*=-1;
			dt.Second*=-1.0;
		}
		obs->tstop.Second = obs->tnow.Second;
		obs->tstop.Minute = obs->tnow.Minute;
		obs->tstop.Hour = obs->tnow.Hour;
		obs->tstop.Day = obs->tnow.Day;
		obs->tstop.Month = obs->tnow.Month;
		obs->tstop.Year = obs->tnow.Year;
		addTime(&(obs->tstop),dt);
	}	
	obs->tstop.time = obs->tstop.Hour/24.0 + obs->tstop.Minute/24.0/60.0 + obs->tstop.Second/24.0/60.0/60.0;
	
	// step
	commentLine = 1;
	while (commentLine)
	{
		getline(fp,line,299);
		if (line[0] != '#')
			commentLine = 0;
	}
	sscanf(line,"%lf",&(obs->tstep));	

	// baseline
	commentLine = 1;
	while (commentLine)
	{
		getline(fp,line,299);
		if (line[0] != '#')
			commentLine = 0;
	}
	sscanf(line,"%lf %lf %lf",&(obs->Xlam),&(obs->Ylam),&(obs->Zlam));
	fclose(fp);

	nconv = obs->timeZone + obs->daylightSavings;
	printf("Observing from %s\n\tlongitude = %.5f\n\tlatitude = %.5f\n\talt = %.2f m\n\thorizon = %.2f deg\n\n",obs->Name,obs->lng,obs->lat,obs->alt,obs->Horizon);
	printf("Start (UTC@%d):       %02d/%02d/%02d  %02d:%02d:%02d\n",nconv,obs->tstart.Month,obs->tstart.Day,obs->tstart.Year,
		   obs->tstart.Hour,obs->tstart.Minute,(int)obs->tstart.Second);
	printf("Stop (UTC@%d):        %02d/%02d/%02d  %02d:%02d:%02d\n",nconv,obs->tstop.Month,obs->tstop.Day,obs->tstop.Year,
		   obs->tstop.Hour,obs->tstop.Minute,(int)obs->tstop.Second);
	printf("Step:                 %.3lf [min]\n",obs->tstep);
	printf("Footprint (UTC@%d):   %02d/%02d/%02d  %02d:%02d:%02d\n",nconv,obs->tnow.Month,obs->tnow.Day,obs->tnow.Year,
		   obs->tnow.Hour,obs->tnow.Minute,(int)obs->tnow.Second);
	printf("Baseline distances:   %lf, %lf, %lf [wavelengths]\n\n",obs->Xlam,obs->Ylam,obs->Zlam);
	fclose(fp);

	return 1;
}

int readLocation(struct observer *obs)
{
	int commentLine, lline;
	double a, b, c;
	char line[120], code[20], tmpLat[20], tmpLng[20];
	FILE *fp;
	
	fp = fopen("obs.dat","r");
	if (fp==NULL)
	{
		printf("obs.dat not found.\n");
		return 0;
	}
	
	// get to observatory listing
	commentLine = 1;
	while (commentLine)
	{
		getline(fp,line,299);
		if (line[0] == '!')
			commentLine = 0;
	}
	// find desired observatory
	lline = 50;
	while (!feof(fp) && lline>30)
	{
		getline(fp,line,299);
		lline = strlen(line);
		sscanf(line,"%s",code);
		if (!strcmp(code,obs->Code))
			lline = 0;
	}
	sscanf(line,"%s %s %s %s %lf %lf %d",obs->Code,obs->Name,tmpLat,tmpLng,&(obs->alt),&(obs->Horizon),&(obs->timeZone));
	sscanf(tmpLat,"%lf:%lf:%lf",&a,&b,&c);
	obs->lat = a + b/60.0 + c/3600.0;
	sscanf(tmpLng,"%lf:%lf:%lf",&a,&b,&c);
	obs->lng = a + b/60.0 + c/3600.0;
	fclose(fp);
	
	return 1;
}
	
void addTime(struct satTime *base_t, struct satTime diff_t)
{
	int leap, inc, yr, mn, dy;
	
	base_t->Second += diff_t.Second;
	base_t->Minute += diff_t.Minute;
	base_t->Hour   += diff_t.Hour;
	base_t->Day    += diff_t.Day;

	// fix up Seconds
	inc = 0;
	while (base_t->Second < 0.0) {
		base_t->Second+=60.0;
		inc-=1;
	}
	while (base_t->Second >= 60.0) {
		base_t->Second-=60.0;
		inc+=1;
	}
	base_t->Minute+=inc;
	
	// fix up Minutes
	inc = 0;
	while (base_t->Minute < 0) {
		base_t->Minute+=60;
		inc-=1;
	}
	while (base_t->Minute >= 60) {
		base_t->Minute-=60;
		inc+=1;
	}
	base_t->Hour+=inc;
	
	// fix up Hours
	inc = 0;
	while (base_t->Hour < 0) {
		base_t->Hour+=24;
		inc-=1;
	}
	while (base_t->Hour >= 24) {
		base_t->Hour-=24;
		inc+=1;
	}
	base_t->Day+=inc;
	
	// fix up Days Months Years ==> assume only one month boundary for now
	yr = base_t->Year;
	dy = base_t->Day;
	mn = base_t->Month;
	leap = !(yr%4) && yr%100  ||  !(yr%400);
	if ( base_t->Day < 1 )
	{
		if ( mn==2 || mn==4 || mn==6 || mn==8 || mn==9 || mn==11 )
		{
			base_t->Month-=1;
			base_t->Day+=31;
		}
		else if ( (mn==5 || mn==7 || mn==10 || mn==12) && dy==31)
		{
			base_t->Month-=1;
			base_t->Day+=30;
		}
		else if (mn==1)
		{
			base_t->Month=12;
			base_t->Day+=31;
			base_t->Year-=1;
		}
		else if (mn==3)
		{
			base_t->Month=2;
			base_t->Day+=(28+leap);
		}
	}
	else if  ( base_t->Day > (28+leap) )
	{
		if ( (mn==1 || mn==3 || mn==5 || mn==7 || mn==8 || mn==10) && dy>31 )
		{
			base_t->Month+=1;
			base_t->Day = base_t->Day - 31;
		}
		else if ( (mn==4 || mn==6 || mn==9 || mn==11) && dy>30 )
		{
			base_t->Month+=1;
			base_t->Day = base_t->Day - 30;
		}
		else if ( mn==2 &&  dy>28+leap)
		{
			base_t->Month=3;
			base_t->Day=base_t->Day - (28+leap);
		}
		else if (mn==12 && dy>31)
		{
			base_t->Month=1;
			base_t->Day = base_t->Day - 31;
			base_t->Year+=1;
		}
	}
	return;
}

void resetTime(struct satTime *t)
{
	t->Year = 0;
	t->Month = 0;
	t->Day = 0;
	t->Hour = 0;
	t->Minute = 0;
	t->Second = 0.0;
	t->time = 0.0;
	return;
}
					
int invxyz(double x, double y, double z, double *lng, double *lat, double *h)
{
	double htmp, lattmp, lat0, s[3];
	double rsat, rho, rhs;
	double rpo=6356.752, req=6378.137, e=0.08181919;

	*lng = atan2(y,x);

	rho = sqrt(SQ(x) + SQ(y));
	rsat = sqrt(SQ(x) + SQ(y) + SQ(z));
	lat0 = asin(z/rsat);
	lattmp = lat0;
	htmp = rsat/req - rpo/sqrt( SQ(rpo*cos(lattmp)) + SQ(req*sin(lattmp)) );
	*h = htmp*req;
	*lat = lattmp;

	s[2] = 1.0E8;
	rhs = x/(req*cos(*lng));
	for (lattmp=lat0 - 0.2; lattmp<lat0+0.2; lattmp+=0.001)
	{
		s[0] = z/(req*tan(lattmp)) + e*e*cos(lattmp)/sqrt(1.0 - SQ(e*sin(lattmp)));
		s[1] = fabs(s[0] - rhs);
		if (s[1] < s[2])
		{
			s[2] = s[1];
			*lat = lattmp;
		}
	}
	*h = x/(cos(*lat)*cos(*lng)) - req/sqrt(1.0 - SQ(e*sin(*lat)));
	return 1;
}

int otherTerms(struct observer obs, double jd, double *ro, struct observer *subsat,
			   double *Az, double *El, double *RA, double *Dec)
{
	int i,j,k;
	double gmst, rad, TmpDbl[3], HA,lst;
	double sslng, sslat, h, range, rsat;
	struct GeoPosition pos;
	double x, y, lmb;
	
	rad = 180.0/PI;
	
	gmst = gstime(jd);
	lst = Time2LST(obs.lng,-99.9,jd,i,j,k);
	
	//rotate to TEME to PEF
	TmpDbl[0] = ro[0];
	TmpDbl[1] = ro[1];
	ro[0] =  cos(gmst)*TmpDbl[0] + sin(gmst)*TmpDbl[1];
	ro[1] = -sin(gmst)*TmpDbl[0] + cos(gmst)*TmpDbl[1];
	
	rsat = sqrt(ro[0]*ro[0] + ro[1]*ro[1] + ro[2]*ro[2]);
	invxyz(ro[0],ro[1],ro[2],&sslng,&sslat,&h);
	
	subsat->lng = sslng*rad;
	subsat->lat = sslat*rad;
	subsat->h = h;
	subsat->alt = rsat;
	
	pos.Longitude = obs.lng;
	pos.Latitude = obs.lat;
	pos.Altitude = obs.alt;
	
	subazel(subsat->lng,subsat->lat,Az,El,&range,rsat,&pos);
	subsat->range = range;
	Eq2Hor(*Az,*El,&HA,Dec,obs.lat);
	*RA = 15.0*lst - HA;
	while (*RA < 0.0)
		*RA+=360.0;
	
	return 1;
}

int writeFootprint(double lngs, double lats, double rsat)
{
	double cg, g, latf, lngf, cll;
	FILE *fp;
	fp = fopen("footprint.out","w");
	
	lngs *= PI/180.0;
	lats *= PI/180.0;
	cg = 6378.0/rsat;
	g = acos(cg);
	printf("central angle = %lf\n",g*180.0/PI);
	for (latf=lats-g+1.0E-6; latf<lats+g; latf+=0.01)
	{
		cll = cg/cos(latf)/cos(lats) - tan(latf)*tan(lats);
		lngf = lngs - acos(cll);
		fprintf(fp,"%lf\t%lf\tup\n",lngf*180.0/PI,latf*180.0/PI);
	}
	for (latf=lats+g-1.0E-6; latf>lats-g; latf-=0.01)
	{
		cll = cg/cos(latf)/cos(lats) - tan(latf)*tan(lats);
		lngf = lngs + acos(cll);
		fprintf(fp,"%lf\t%lf\tdn\n",lngf*180.0/PI,latf*180.0/PI);
	}
	latf=lats-g+1.0E-6;
	cll = cg/cos(latf)/cos(lats) - tan(latf)*tan(lats);
	lngf = lngs - acos(cll);
	fprintf(fp,"%lf\t%lf\tup\n",lngf*180.0/PI,latf*180.0/PI);
	fclose(fp);
	return 1;
}
		
