import numpy as np

ElLimit = 10.0

Az = []
El = []
RA = []
Dec = []
fp = open('snapshotAERD.out','r')
for line in fp:
    data = line.split()
    Az.append(float(data[0]))
    El.append(float(data[1]))
    RA.append(float(data[2]))
    Dec.append(float(data[3]))
Az = np.array(Az)
El = np.array(El)
RA = np.array(RA)
Dec = np.array(Dec)

x = []
y = []
z = []
a = []
fp = open('snapshotXYZ.out','r')
for line in fp:
    data = line.split()
    x.append(float(data[0]))
    y.append(float(data[1]))
    z.append(float(data[2]))
    a.append(float(data[3]))
x = np.array(x)
y = np.array(y)
z = np.array(z)
a = np.array(a)

lng = []
lat = []
h = []
r = []
fp = open('snapshotLLH.out','r')
for line in fp:
    data = line.split()
    lng.append(float(data[0]))
    lat.append(float(data[1]))
    h.append(float(data[2]))
    r.append(float(data[3]))
lng = np.array(lng)
lat = np.array(lat)
h = np.array(h)
r = np.array(r)

    
