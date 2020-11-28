from astropy.coordinates import EarthLocation
from astropy import units as u
import numpy as np
from argparse import Namespace


class Track:

    C = 3E8
    KB = 1.38E-23
    satpar = ['since', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'H']

    def __init__(self, fname):
        self.fname = fname
        for sp in self.satpar:
            setattr(self, sp, [])
        self.read_file()

    def read_file(self):
        self.header = ''
        with open(self.fname, 'r') as fp:
            for line in fp:
                if 'xx' in line:
                    self.header += (line.strip() + ' ')
                    if 'period' in line:
                        self.period = float(self.header.split('=')[1].split()[0])
                    continue
                data = line.split()
                df = [float(data[i]) for i in range(len(self.satpar))]
                for i, sp in enumerate(self.satpar):
                    getattr(self, sp).append(df[i])
        for sp in self.satpar:
            if sp == 'since':
                self.since = np.array(self.since) * 60.0  # convert to sec
                continue
            setattr(self, sp, np.array(getattr(self, sp))*1000.0)  # convert to m

    def calc(self, freq, loc=None):
        if loc is not None:
            self.location(loc)
        self.view()
        self.rates(freq)
        self.subsat()

    def location(self, name, lon=None, lat=None, alt=None):
        if isinstance(name, dict):
            self.observer, lon, lat, alt = name['name'], name['lon'], name['lat'], name['alt']
        else:
            self.observer = name
        self.loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=alt*u.m)

    def view(self):
        self.R = Namespace(x=(self.x-self.loc.x.value),
                           y=(self.y-self.loc.y.value),
                           z=(self.z-self.loc.z.value))
        self.D = np.sqrt(self.R.x**2 + self.R.y**2 + self.R.z**2)
        robs = np.sqrt(self.loc.x.value**2 + self.loc.y.value**2 + self.loc.z.value**2)
        cxyz = self.loc.x.value*self.R.x + self.loc.y.value*self.R.y + self.loc.z.value*self.R.z
        self.za = np.rad2deg(np.arccos(cxyz / (robs*self.D)))
        self.viewable = len(np.where(self.za < 88.0)[0]) > 1

    def vis(self, arr, val=0.0):
        varr = 1.0 * np.array(arr)
        varr[np.where(self.za > 90.0)] = val
        return varr

    def rates(self, f=982E6):
        self.f = f
        self.V = (self.vx*self.R.x + self.vy*self.R.y + self.vx*self.R.z) / self.D
        self.doppler = (np.array(self.V) / self.C) * f
        self.drift = [0.0]
        for i in range(1, len(self.doppler)):
            dt = self.since[i] - self.since[i-1]
            dd = self.doppler[i] - self.doppler[i-1]
            self.drift.append(dd/dt)
        self.drift[0] = self.drift[1]

    def subsat(self):
        self.lon = np.rad2deg(np.arctan2(self.y, self.x))
        self.lat = np.rad2deg(np.arcsin(self.z / self.D))

    def waterfall(self, pwr=4.0*np.pi, Tsys=50.0, BW=2.0):
        self.Rxfreq = self.f + self.doppler
        self.Rxpwr = self.vis(pwr / (4.0 * np.pi * self.D**2))
        ds2 = (self.doppler.max() - self.doppler.min()) / 2.0
        numch = int(np.ceil((self.Rxfreq.max() - self.Rxfreq.min() + 2.0*ds2)/BW))
        self.ch = np.linspace(self.Rxfreq.min()-ds2, self.Rxfreq.max()+ds2, numch)
        self.chnum = (self.Rxfreq - self.Rxfreq.min())/BW + 0.5
        self.wf = []
        int_time = self.since[1] - self.since[0]  # not necessarily a good approach
        sm = np.sqrt(BW * int_time)
        pw = self.KB * BW * Tsys
        for i in range(len(self.Rxfreq)):
            pn = np.random.rayleigh(pw/sm, numch)
            pn[int(self.chnum[i])] += self.Rxpwr[i]
            self.wf.append(list(np.log(pn)))
        self.wf = np.array(self.wf)

# int writeFootprint(double lngs, double lats, double rsat)
# {
# 	double cg, g, latf, lngf, cll;
# 	FILE *fp;
# 	fp = fopen("footprint.out","w");
#
# 	lngs *= PI/180.0;
# 	lats *= PI/180.0;
# 	cg = 6378.0/rsat;
# 	g = acos(cg);
# 	for (latf=lats-g+1.0E-6; latf<lats+g; latf+=0.01)
# 	{
# 		cll = cg/cos(latf)/cos(lats) - tan(latf)*tan(lats);
# 		lngf = lngs - acos(cll);
# 		//fprintf(fp,"%lf\t%lf\tup\n",lngf*180.0/PI,latf*180.0/PI);
#         fprintf(fp,"%lf\t%lf\n",lngf*180.0/PI,latf*180.0/PI);
# 	}
# 	for (latf=lats+g-1.0E-6; latf>lats-g; latf-=0.01)
# 	{
# 		cll = cg/cos(latf)/cos(lats) - tan(latf)*tan(lats);
# 		lngf = lngs + acos(cll);
# 		//fprintf(fp,"%lf\t%lf\tdn\n",lngf*180.0/PI,latf*180.0/PI);
#         fprintf(fp,"%lf\t%lf\n",lngf*180.0/PI,latf*180.0/PI);
# 	}
# 	latf=lats-g+1.0E-6;
# 	cll = cg/cos(latf)/cos(lats) - tan(latf)*tan(lats);
# 	lngf = lngs - acos(cll);
# 	//fprintf(fp,"%lf\t%lf\tup\n",lngf*180.0/PI,latf*180.0/PI);
#     fprintf(fp,"%lf\t%lf\n",lngf*180.0/PI,latf*180.0/PI);
# 	fclose(fp);
# 	return 1;
# }
