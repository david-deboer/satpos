/*****************************main.c***************************/
/* copy of satpos, but massively truncated to just get the periods of all satellites in database*/
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
	long mcnt, nup, finished, validFILE;
	int transit, printOutTransit;
	char TLEFile[55], SCName[40], SCSearch[40], line[100];
	FILE *fp, *fp2;
	/////////////////This is cut-and-pasted in from testcpp.cpp///////////////
	char str[2];
	double ro[3];
	double vo[3];
	char typeinput, typerun, opsmode;
	gravconsttype  whichconst;
	int whichcon;
	FILE *outfile;
	// ----------------------------  locals  -------------------------------
	double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper, period;
	double sec,  jd, rad, tsince, startmfe, stopmfe, deltamin, jdstart, jdstop;
	double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;
	int  year, mon, day, hr, min;
	char longstr1[130];
	typedef char str3[4];
	str3 monstr[13];
	char longstr2[130];
	elsetrec satrec;
	
	printf("%s\n",SGP4Version );
    system("ls -1 tle/*.tle > tlefiles.dat");
    opsmode = 'a';
    typeinput = 'e';
    typerun = 'c';  // use 'c' (see testcpp.cpp) to just read in tle's.  start/stop gets set below (so overwrites the 'c' settings in sgp4io)
    whichcon = 72;
	if (whichcon == 721) whichconst = wgs72old;
	if (whichcon == 72) whichconst = wgs72;
	if (whichcon == 84) whichconst = wgs84;
	getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

    outfile = fopen("periods.dat","w");
	fp=fopen("tlefiles.dat","r");
    validFILE = 1;
	while (validFILE)
	{
        getline(fp,TLEFile,55);
        if (strlen(TLEFile)<4)
        {
            validFILE = 0;
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
                    printf("%s: %s ==> %lf\n",TLEFile,SCName,period);
                    fprintf(outfile,"%s:%s | %lf\n",TLEFile,SCName,period);
                }
            }
            fclose(fp2);
        }
    }
    fclose(fp);

	return 1;
}
