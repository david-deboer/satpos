from astropy.coordinates import EarthLocation
from astropy import units as u
import numpy as np


class Track:

    satpar = ['since', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'distance']

    def __init__(self, fname):
        self.fname = fname
        for sp in self.satpar:
            setattr(self, sp, [])

    def read_file(self):
        self.name = ''
        with open(self.fname, 'r') as fp:
            for line in fp:
                if 'xx' in line:
                    self.name += (line.strip() + ' ')
                    if 'Period' in line:
                        self.period = float(self.name.split('=')[1].split()[0])
                    continue
                data = line.split()
                df = [float(data[i]) for i in range(len(self.satpar))]
                for i, sp in enumerate(self.satpar):
                    getattr(self, sp).append(df[i])
        for sp in self.satpar:
            setattr(self, sp, np.array(getattr(self, sp)))

    def location(self, name, lon=None, lat=None, alt=None):
        if isinstance(name, dict):
            lon, lat, alt = name['lon'], name['lat'], name['alt']
        self.loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=alt*u.m)

    def separation_angle(self):
        robs = np.sqrt(self.loc.x.value**2 + self.loc.y.value**2 + self.loc.z.value**2)
        rxyz = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        cosxyz = self.loc.x.value * self.x + self.loc.y.value * self.y + self.loc.z.value * self.z
        self.ang = np.rad2deg(np.arccos(cosxyz / (robs * rxyz)))
        self.viewable = len(np.where(self.ang < 90.0)) > 0
