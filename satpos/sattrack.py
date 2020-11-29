import numpy as np
from argparse import Namespace
from my_ephem import ephem


C = 3E8
KB = 1.38E-23
REARTH = 6371.0E3


class Track:
    satpar = ['since', 'x', 'y', 'z', 'vx', 'vy', 'vz']

    def __init__(self, fname):
        self.fname = fname
        for sp in self.satpar:
            setattr(self, sp, [])
        self.read_file()

    def read_file(self):
        self.headerlines = []
        with open(self.fname, 'r') as fp:
            for line in fp:
                if line[0] == '#':
                    self.headerlines.append(line.strip()[1:])
                    continue
                data = line.strip().split()
                df = [float(data[i]) for i in range(len(self.satpar))]
                for i, sp in enumerate(self.satpar):
                    getattr(self, sp).append(df[i])
        for sp in self.satpar:
            if sp == 'since':
                self.since = np.array(self.since) * 60.0  # convert to sec
                continue
            setattr(self, sp, np.array(getattr(self, sp))*1000.0)  # convert to m
        self.parse_header()

    def parse_header(self):
        self.headers = {'scname': Namespace(text='spacecraft name', type=str, nval=0),
                        'satnum': Namespace(text='satnum', type=int, nval=0),
                        'period': Namespace(text='period', type=float, nval=0),
                        'sublon': Namespace(text='starting lon lat h', type=float, nval=0),
                        'sublat': Namespace(text='starting lon lat h', type=float, nval=1),
                        'height': Namespace(text='starting lon lat h', type=float, nval=2)}
        for hdrline in self.headerlines:
            for hdr, valns in self.headers.items():
                if valns.text in hdrline:
                    try:
                        data = hdrline.split(':')[1].strip().split()
                        setattr(self, hdr, valns.type(data[valns.nval]))
                        self.headers[hdr].value = valns.type(data[valns.nval])
                    except IndexError:
                        continue
                    break

    def calc(self, freq, loc):
        self.location(loc)
        self.view()
        self.rates(freq)
        self.subsat()

    def location(self, name, lon=None, lat=None, alt=None):
        self.location = ephem.Location(name, lon, lat, alt)
        self.loc = Namespace(x=self.location.loc.x.value,
                             y=self.location.loc.y.value,
                             z=self.location.loc.z.value)

    def view(self):
        self.R = Namespace(x=(self.x-self.loc.x),
                           y=(self.y-self.loc.y),
                           z=(self.z-self.loc.z))
        self.D = np.sqrt(self.R.x**2 + self.R.y**2 + self.R.z**2)
        robs = np.sqrt(self.loc.x**2 + self.loc.y**2 + self.loc.z**2)
        cxyz = self.loc.x*self.R.x + self.loc.y*self.R.y + self.loc.z*self.R.z
        self.za = np.rad2deg(np.arccos(cxyz / (robs*self.D)))
        self.viewable = len(np.where(self.za < 88.0)[0]) > 1

    def vis(self, arr, val=0.0):
        varr = 1.0 * np.array(arr)
        varr[np.where(self.za > 90.0)] = val
        return varr

    def rates(self, f=982E6):
        self.f = f
        self.V = (self.vx*self.R.x + self.vy*self.R.y + self.vx*self.R.z) / self.D
        self.doppler = (np.array(self.V) / C) * f
        self.drift = [0.0]
        for i in range(1, len(self.doppler)):
            dt = self.since[i] - self.since[i-1]
            dd = self.doppler[i] - self.doppler[i-1]
            self.drift.append(dd/dt)
        self.drift[0] = self.drift[1]

    def subsat(self):
        self.lon = np.rad2deg(np.arctan2(self.y, self.x))
        self.rsat = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        r = self.z / self.rsat
        self.lat = np.rad2deg(np.arcsin(r))

    def footprint(self, i):
        return np.array(footprint(self.lon[i], self.lat[i], self.rsat[i]))

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
        pw = KB * BW * Tsys
        for i in range(len(self.Rxfreq)):
            pn = np.random.rayleigh(pw/sm, numch)
            pn[int(self.chnum[i])] += self.Rxpwr[i]
            self.wf.append(list(np.log(pn)))
        self.wf = np.array(self.wf)


def footprint(lngs, lats, rsat):
    lngs = np.deg2rad(lngs)
    lats = np.deg2rad(lats)
    cg = REARTH/rsat
    g = np.arccos(cg)
    footprint = []
    for latf in np.arange(lats-g, lats+g, 0.01):
        cll = cg/np.cos(latf)/np.cos(lats) - np.tan(latf)*np.tan(lats)
        lngf = lngs - np.arccos(cll)
        footprint.append([np.rad2deg(lngf), np.rad2deg(latf)])
    for latf in np.arange(lats+g, lats-g, -0.01):
        cll = cg/np.cos(latf)/np.cos(lats) - np.tan(latf)*np.tan(lats)
        lngf = lngs + np.arccos(cll)
        footprint.append([np.rad2deg(lngf), np.rad2deg(latf)])
    latf = lats-g
    cll = cg/np.cos(latf)/np.cos(lats) - np.tan(latf)*np.tan(lats)
    lngf = lngs - np.arccos(cll)
    footprint.append([np.rad2deg(lngf), np.rad2deg(latf)])
    return footprint
