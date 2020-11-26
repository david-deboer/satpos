# from my_ephem import base_ephem

satellites = {}
sats_by_file = {}
total_count = 0
with open('tle/master.dat', 'r') as fp:
    for line in fp:
        fname = f"tle/{line.split(':')[0]}"
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

with open('tle/complete_set.tle', 'w') as fp:
    for this_sat in satellites.keys():
        print(satellites[this_sat]['line0'], file=fp)
        print(satellites[this_sat]['line1'], file=fp)
        print(satellites[this_sat]['line2'], file=fp)

with open('check_all.sh', 'w') as fp:
    for i in range(len(satellites.keys())):
        print("./satpos {} tle/complete_set.tle".format(i+1), file=fp)
