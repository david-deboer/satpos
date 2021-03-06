from argparse import Namespace
from satpos import sattrack
from os.path import join


def filter_drift(loc, f=982e6, rng=[0.025, 0.04],
                 trackfilelist='viewable.csv', path='output', addtag='out'):
    outsats = open('outsats.csv', 'w')
    drifts = {}
    with open(trackfilelist, 'r') as satslist:
        for line in satslist:
            data = line.strip().split(',')
            if data[0] == 'file':
                continue
            fname = data[0]
            if addtag is not None:
                fname = f"{fname}.{addtag}"
            fname = join(path, fname)
            s = sattrack.Track(fname)
            s.view(loc)
            s.rates(f)
            for i, drift in enumerate(s.drift):
                if s.za[i] < 90.0 and drift > rng[0] and drift < rng[1]:
                    drifts.setdefault(s.satnum, Namespace(t=[], drift=[], za=[], period=[]))
                    drifts[s.satnum].t.append(s.since[i])
                    drifts[s.satnum].drift.append(drift)
                    drifts[s.satnum].za.append(s.za[i])
                    drifts[s.satnum].period.append(s.period)
                    print(f"{fname},{s.satnum},{s.since[i]},{drift},{s.x[i]},{s.y[i]},{s.z[i]},"
                          f"{s.za[i]},{s.period}", file=outsats)
    outsats.close()
    return drifts


def find_viewable(loc, rng=None, trackfilelist='ls.out', path='output'):
    viewable = open('viewable.csv', 'w')
    notviewable = open('notviewable.csv', 'w')
    hdrstr = "file,scname,satnum,orbit,period,sublon,zamin,zamax"
    print(hdrstr, file=viewable)
    print(hdrstr, file=notviewable)
    count = Namespace(leo=0, meo=0, geo=0, deep=0, other=0, viewable=0, notviewable=0)
    with open(trackfilelist, 'r') as fp:
        i = 0
        for line in fp:
            if rng is not None:
                if i < rng[0] or i > rng[1]-1:
                    continue
            fname = line.strip()
            s = sattrack.Track(join(path, fname))
            s.view(loc)
            if s.period > 1500.0:
                count.deep += 1
                orbit = 'deep'
            elif s.period < 1450.0 and s.period > 1420.0:
                count.geo += 1
                orbit = 'geo'
            elif s.period < 1000.0 and s.period > 500.0:
                count.meo += 1
                orbit = 'meo'
            elif s.period < 200.0:
                count.leo += 1
                orbit = 'leo'
            else:
                count.other += 1
                orbit = 'other'
            try:
                zamin = s.za.min()
                zamax = s.za.max()
            except ValueError:
                zamin = '!'
                zamax = '!'
            fnp = fname.split('.')[0]
            pline = f"{fnp},{s.scname},{s.satnum},{orbit},{s.period},{s.sublon},{zamin},{zamax}"
            if s.viewable:
                print(pline, file=viewable)
                count.viewable += 1
            else:
                print(pline, file=notviewable)
                count.notviewable += 1
            i += 1
    viewable.close()
    notviewable.close()

    print(f"LEO: {count.leo}")
    print(f"MEO: {count.meo}")
    print(f"GEO: {count.geo}")
    print(f"DEEP: {count.deep}")
    print(f"OTHER: {count.other}")
    print(f"VIEWABLE: {count.viewable}")
    print(f"NOTVIEWABLE: {count.notviewable}")


def generate_check_all(fname, tot):
    with open('check_all.sh', 'w') as fp:
        for i in range(tot):
            print(f"satpos {fname} {i+1}", file=fp)


def generate_complete_set(path='tle', fmname='master.dat'):
    import os.path
    satellites = {}
    sats_by_file = {}
    total_count = 0
    flistname = os.path.join(path, fmname)
    with open(flistname, 'r') as fp:
        for line in fp:
            fname = os.path.join(path, f"{line.strip().split(':')[0]}")
            sats_by_file[fname] = []
            with open(fname, 'r') as fptle:
                for tleline in fptle:
                    data = tleline.split()
                    if data[0] not in ['1', '2']:
                        scname = tleline.strip()
                        line0 = tleline.strip('\n')
                    elif data[0] == '1':
                        line1 = tleline.strip('\n')
                    elif data[0] == '2':
                        line2 = tleline.strip('\n')
                        key = data[1]
                        total_count += 1
                        satellites.setdefault(key, {'scname': scname, 'files': []})
                        satellites[key]['line0'] = line0
                        satellites[key]['line1'] = line1
                        satellites[key]['line2'] = line2
                        satellites[key]['files'].append(fname)
    satlist = list(satellites.keys())

    print("Total satellites listed: {}".format(total_count))
    print("Total unique satellites:  {}".format(len(satlist)))

    # Still do, even though aren't using
    for i in range(len(satlist)):
        this_sat = satlist.pop()
        for fname in sats_by_file.keys():
            if fname in satellites[this_sat]['files']:
                satdes = '{}:{}'.format(this_sat, satellites[this_sat]['scname'])
                sats_by_file[fname].append(satdes)
                break
    with open('tle/completeset.tle', 'w') as fp:
        for this_sat in satellites.keys():
            print(satellites[this_sat]['line0'], file=fp)
            print(satellites[this_sat]['line1'], file=fp)
            print(satellites[this_sat]['line2'], file=fp)
