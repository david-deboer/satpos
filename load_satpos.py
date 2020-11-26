
class Track:

    satpar = ['since', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'distance']

    def __init__(self, name=''):
        self.name = name
        for sp in self.satpar:
            setattr(self, sp, [])


sat = Track()


with open("satpos.out", 'r') as fp:
    for line in fp:
        if 'xx' in line:
            sat.name += (line.strip() + ' ')
            period = float(sat.name.split('=')[1].split()[0])
            continue
        data = line.split()
        df = [float(data[i]) for i in range(len(sat.satpar))]
        for i, sp in enumerate(sat.satpar):
            getattr(sat, sp).append(df[i])
