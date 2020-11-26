/*****************************main.c***************************/

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
	int argSC=0, argTLE=0, i, j, k, n, valid;
	long mcnt, nup;
	int transit, printOutTransit;
    char TLEprefix[25], TLEFile[25], SCName[40], SCSearch[10], line[50];
    char outname[70], outsatpos[70];
	double RA, HA, Dec, Az, El, dEl, dAz, oldel, oldaz, olddel, period;
	double DOY, olddoy, jday1, jnow, timeconv[3], lngconv[3];
	double SCVector[15];
	double w, cw, cDcH, cDsH, sD;
	FILE *fp, *fpFringe, *fpTransit, *fpSubsat;
	struct observer obs, subsat, now;
	/////////////////This is cut-and-pasted in from testcpp.cpp///////////////
	char str[2];
	double ro[3];
	double vo[3];
	char typeinput, typerun, opsmode;
	gravconsttype  whichconst;
	int whichcon;
	FILE *outfile;
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
	
	printf("%s\n", SGP4Version);
    if (argc > 1)
    {
        argTLE = 1;
        strcpy(TLEprefix, argv[1]);
        printf("Using %s  ", TLEprefix);
    }
    if (argc > 2)
	{
		argSC = 1;
        sscanf(argv[2], "%d", &k);
        printf(" %d\n", k);
	}

	// Get satellite TLE's
	if (!argTLE)
	{
		printf("------TLE files------\n");
		system("ls -1 tle/*.tle");
		printf("\t--------> TLE Filename (leave off the tle's):  ");
		scanf("%s",TLEprefix);
		
	}
    sprintf(TLEFile,"tle/%s.tle",TLEprefix);
	fp=fopen(TLEFile,"r");
	if (fp==NULL)
    {
		printf("%s not found.\n",TLEFile);
		exit(1);
	}
    if (!argSC)
    {
        i=0;
        printf("\n\nSatellites...\n");
        while(!feof(fp))
        {
            getline(fp,SCName,29);
            getline(fp,longstr1,100);
            n=getline(fp,longstr2,100);
            ++i;
            if (n > 10)
                printf("\t%d %s\n",i,SCName);
        }
		rewind(fp);
        printf("Input spacecraft number:  ");
		scanf("%d",&k);
	}
	for(i=0;i<k;++i)
	{
		getline(fp, SCName, 30);
		getline(fp, longstr1, 110);
		getline(fp, longstr2, 110);
	}
	fclose(fp);
 
	for(i=strlen(SCName)-1;isspace(SCName[i]);--i)
		SCName[i]='\0';
	printf("\n\n%s\n%s\n%s\n\n",SCName,longstr1,longstr2);
	readObserver(&obs);
	jday( obs.tstart.Year,1,0,0,0,0.0, jday1 );
    jday( obs.tnow.Year,obs.tnow.Month,obs.tnow.Day,obs.tnow.Hour,obs.tnow.Minute,obs.tnow.Second,jnow );
	
	//opsmode = 'a' best understanding of how afspc code works
	//opsmode = 'i' improved sgp4 resulting in smoother behavior
	opsmode = 'a';
	//typeinput = 'm' input start stop mfe
	//typeinput = 'e' input start stop ymd hms
	//typeinput = 'd' input start stop yr dayofyr
	typeinput = 'e';
	typerun = 'c';  // use 'c' (see testcpp.cpp) to just read in tle's.  start/stop gets set below (so overwrites the 'c' settings in sgp4io)
	
	whichcon = 84;
	if (whichcon == 721) whichconst = wgs72old;
	if (whichcon == 72) whichconst = wgs72;
	if (whichcon == 84) whichconst = wgs84;
	
	getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
	
	// ---------------- setup files for operation ------------------
    sprintf(outname, "sp%s_%d.out", TLEprefix, k);
    outfile = fopen(outname, "w");
	fpSubsat = fopen("subsat.out","w");
	//fp   = fopen("predix.out","w");
	//fpFringe = fopen("proFringe.out","w");
	//fpTransit = fopen("transit.out","w");
    fprintf(outfile, "%d %s xx ", k, SCName);

	// convert the char string to sgp4 elements
	// includes initialization of sgp4
	twoline2rv( longstr1, longstr2, typerun, typeinput, opsmode, whichconst,
			   startmfe, stopmfe, deltamin, satrec );

	// This is pulled out of sgp4io.cpp to change the start/stop times
	jday( obs.tstart.Year,obs.tstart.Month,obs.tstart.Day,obs.tstart.Hour,obs.tstart.Minute,obs.tstart.Second, jdstart );
	jday( obs.tstart.Year,obs.tstop.Month,obs.tstop.Day,obs.tstop.Hour,obs.tstop.Minute,obs.tstop.Second, jdstop );
	startmfe = (jdstart - satrec.jdsatepoch) * 1440.0;
	stopmfe  = (jdstop - satrec.jdsatepoch) * 1440.0;
	deltamin = obs.tstep;
	
	fprintf(outfile, "%ld xx", satrec.satnum);
	// call the propagator to get the initial state vector value
	sgp4 (whichconst, satrec,  0.0, ro,  vo);

	jd = satrec.jdsatepoch;
	invjday( jd, year,mon,day,hr,min, sec );
	printf("stk.v.4.3 \n"); // must use 4.3...
	printf("ScenarioEpoch:  %3i %3s%5i%3i:%2i:%12.9f \n",day,monstr[mon],year,hr,min,sec );
	printf("CoordinateSystem:  TEME-to-PEF\n\n");
	
	// First do it for 'now'
	tsince = (jnow - satrec.jdsatepoch)*1440.0;
	sgp4 (whichconst, satrec,  tsince, ro,  vo);
	otherTerms(obs,jnow,ro,&now,&Az,&El,&RA,&Dec);
    period = 2.0*pi/satrec.no;
    fprintf(outfile, "  period = %12.9f min\n", period);
    printf("Not doing all otherTerms\n");
	printf("Current Position:\n\t(x,y,z,r):  %lf, %lf, %lf, %lf\n",ro[0],ro[1],ro[2],now.alt);
	printf("\t(sub_lng, sub_lat, h, range):  %lf, %lf, %lf %lf\n",now.lng,now.lat,now.h,now.range);
	printf("\t(Az, El, RA, Dec):  %lf, %lf, %lf, %lf\n",Az,El,RA/15.0,Dec);
	
	// initialize variables
	nup = 0;
	mcnt = 0;
	transit = 0;
	dEl = 0.0;
	dAz = 0.0;
	oldel = 0.0;
	oldaz = 0.0;
	olddel = 0.0;
	tsince = startmfe;
	// check so the first value isn't written twice
	if ( fabs(tsince) > 1.0e-8 )
		tsince = tsince - deltamin;
	lngconv[0] = 999.9;
	// ----------------- loop to perform the propagation ----------------
	while ((tsince < stopmfe) && (satrec.error == 0))
	{
		tsince = tsince + deltamin;
		if(tsince > stopmfe)
			tsince = stopmfe;
		sgp4 (whichconst, satrec,  tsince, ro,  vo);
		
		if (satrec.error > 0)
			printf("# *** error: t:= %f *** code = %3d\n",satrec.t, satrec.error);
		if (satrec.error == 0)
		{
			jd = satrec.jdsatepoch + tsince/1440.0;
			invjday( jd, year,mon,day,hr,min, sec );
			DOY = jd - jday1;
			otherTerms(obs,jd,ro,&subsat,&Az,&El,&RA,&Dec);
			lngconv[1] = subsat.lng;
			if (lngconv[0]!=999.9 && fabs(lngconv[1]-lngconv[0])>350.0) 
			{
				lngconv[1] += 360.0;
			}
			lngconv[0] = lngconv[1];
			fprintf(fpSubsat,"%lf\t%lf\n",lngconv[1],subsat.lat);
			if (El > obs.Horizon)
			{
				timeconv[0] = DOY + (obs.timeZone + obs.daylightSavings)/24.0;
				timeconv[1] = hr + min/60.0 + sec/3600.0 + (obs.timeZone + obs.daylightSavings);
				if (timeconv[1]>24.0) timeconv[1]-=24.0;
				if (timeconv[1]<0.0) timeconv[1]+=24.0;
				fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",timeconv[0],timeconv[1],Az,El,subsat.range,dAz,dEl);
				cDcH = cos(obs.lat/rad)*sin(El*pi/180.0) - sin(obs.lat/rad)*cos(El*pi/180.0)*cos(Az*pi/180.0);
				cDsH = -1.0*cos(El*pi/180.0)*sin(Az*pi/180.0);
				sD   = sin(obs.lat/rad)*sin(El*pi/180.0) + cos(obs.lat/rad)*cos(El*pi/180.0)*cos(Az*pi/180.0);
				w = obs.Xlam*cDcH - obs.Ylam*cDsH + obs.Zlam*sD;
				cw = cos(2.0*pi*w);
				//fprintf(fpFringe,"%lf\t%lf\t%lf\t%lf\n",timeconv[0],timeconv[1],w,cw);
				if (mcnt>0L)
				{
					dEl = (El-oldel)/(DOY-olddoy)/(60.0*24.0); // Convert to deg/min 
					dAz = (Az-oldaz)/(DOY-olddoy)/(60.0*24.0);
					if (dEl == 0.0) transit = 1;
					else if (olddel/dEl < 0.0) transit = 1;
					else transit = 0;
				}
				if ( transit )
				{
					if (printOutTransit>NOPRINT)
					{
						if (printOutTransit==PRINTLABEL)
						{
							printf("\nSpacecraft                       day-of-year  time     Az [deg]  El [deg]   range [km]  dAz [deg/min]  dEl [deg/min]\n");
							printOutTransit = YESPRINT;
						}
						printf("%20s transiting at %8.4lf  %7.4lf   %6.2lf     %5.2lf      %7.1lf     %.2lf          %.2lf\n",SCName,timeconv[0],timeconv[1],Az,El,subsat.range,dAz,dEl);
					}
					/*fprintf(fpTransit,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",timeconv[0],timeconv[1],Az,El,subsat.range,dAz,dEl);*/
				}
				++nup;
				oldaz = Az;
				oldel = El;
				olddel = dEl;
				olddoy = DOY;
			}
			fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f",
					tsince,ro[0],ro[1],ro[2],vo[0],vo[1],vo[2]);
            period = 2.0*pi/satrec.no;
			rv2coe(ro, vo, mu, p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper );
			fprintf(outfile, " %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f %5i%3i%3i %2i:%2i:%9.6f\n",
					a, ecc, incl*rad, node*rad, argp*rad, nu*rad,
					m*rad,year,mon,day,hr,min,sec);
			++mcnt;
		} // if satrec.error == 0
	}

	printf("\n%s is up %g%% (%ld/%ld)\n",SCName, 100.0*( (double)nup )/( (double) mcnt ), nup, mcnt );
	printf("\nOutput to:\n");
	//printf("\tpredix.out\n");
	//printf("\tproFringe.out\n");
	//printf("\ttransit.out\n");
	printf("\tsubsat.out\n");
	printf("\tsatpos.out\n");
	//fclose(fp);
	//fclose(fpFringe);
	//fclose(fpTransit);
	fclose(fpSubsat);
	fclose(outfile);
	printf("\n");
	writeFootprint(now.lng,now.lat,now.alt);
	fp = fopen("info.out","w");
	fprintf(fp,"%s | %s\n",obs.Name,SCName);
	fprintf(fp,"%lf %lf\n",obs.lng,obs.lat);
	fprintf(fp,"=%lf %lf\n",now.lng, now.lat);
	fclose(fp);
	//system("python pltSat.py");
	
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
	
	// Year Month Day hour min second (UTC)
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
		obs->tnow.Year = tm_now->tm_year;
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
    fclose(fp);
	sscanf(line,"%lf %lf %lf",&(obs->Xlam),&(obs->Ylam),&(obs->Zlam));

	nconv = obs->timeZone + obs->daylightSavings;
	printf("Observing from %s\n\tlongitude = %.5f\n\tlatitude = %.5f\n\talt = %.2f m\n\thorizon = %.2f deg\n\n",obs->Name,obs->lng,obs->lat,obs->alt,obs->Horizon);
	printf("Start (UTC@%d):       %02d/%02d/%02d  %02d:%02d:%02d\n",nconv,obs->tstart.Month,obs->tstart.Day,obs->tstart.Year,
		   obs->tstart.Hour,obs->tstart.Minute,(int)obs->tstart.Second);
	printf("Stop (UTC@%d):        %02d/%02d/%02d  %02d:%02d:%02d\n",nconv,obs->tstop.Month,obs->tstop.Day,obs->tstop.Year,
		   obs->tstop.Hour,obs->tstop.Minute,(int)obs->tstop.Second);
	printf("Step:                 %.3lf [min]\n",obs->tstep);
	printf("Footprint (UTC@%d):   %02d/%02d/%02d  %02d:%02d:%02d\n",nconv,obs->tnow.Month,obs->tnow.Day,obs->tnow.Year,
		   obs->tnow.Hour,obs->tnow.Minute,(int)obs->tnow.Second);
	//printf("Baseline distances:   %lf, %lf, %lf [wavelengths]\n\n",obs->Xlam,obs->Ylam,obs->Zlam);
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
	
	/* Commenting out 11/25/20 to get working...
    subazel(subsat->lng,subsat->lat,Az,El,&range,rsat,&pos);
	subsat->range = range;
	Eq2Hor(*Az,*El,&HA,Dec,obs.lat);
	*RA = 15.0*lst - HA;
	while (*RA < 0.0)
		*RA+=360.0;
	*/
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
		
