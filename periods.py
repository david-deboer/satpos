import matplotlib.pyplot as plt
import numpy as np

fp = open('periods.dat','r')

sats = {}
for line in fp:
    data = line.split('|')
    sats[data[0]] = float(data[1])
fp.close()

##for s in sats:
##    print '%s ==> %lf' % (s, sats[s])

sortedPeriods = []
for val in sorted(sats.values(),key=float):
    sortedPeriods.append(val)

pvals = np.linspace(0.0,sortedPeriods[-1],2000)
cumulative = np.zeros(len(pvals))

for i,cv in enumerate(pvals):
    for s in sats:
        if sats[s] < cv:
            cumulative[i]+=1.0

plt.figure(3)
plt.plot(pvals,cumulative)

deriv = np.zeros(len(pvals)-1)
dvals = np.zeros(len(deriv))
dx = pvals[1] - pvals[0]
for i in range(len(pvals)-1):
    dy = cumulative[i+1] - cumulative[i]
    deriv[i] = dy/dx
    dvals[i] = (pvals[i+1]+pvals[i])/2.0

plt.figure(4)
plt.plot(dvals,deriv)
plt.figure(5)
#rndderiv = np.around(deriv)
#plt.stem(dvals,rndderiv)
ceilderiv = np.ceil(deriv)
plt.plot(dvals,ceilderiv)
#plt.plot(dvals,ceilderiv,'.')

fp = open('leo.dat','w')
leo = {}
for s in sats:
    if sats[s] < 400.0:
        leo[s] = sats[s]
sortedLeo = []
for val in sorted(leo.keys(),key=str.lower):
    sortedLeo.append(val)
i = 1
for s in sortedLeo:
    sw = '%03d %-40s  %lf\n' % (i,s,leo[s])
    fp.write(sw)
    i+=1
fp.close()

fp = open('meo.dat','w')
meo = {}
for s in sats:
    if sats[s] >= 400.0 and sats[s]<=1300.0:
        meo[s] = sats[s]
sortedMeo = []
for val in sorted(meo.keys(),key=str.lower):
    sortedMeo.append(val)
i = 1
for s in sortedMeo:
    sw = '%03d %-40s  %lf\n' % (i,s,meo[s])
    fp.write(sw)
    i+=1
fp.close()

fp = open('geo.dat','w')
geo = {}
for s in sats:
    if sats[s] > 1300.0:
        geo[s] = sats[s]
sortedGeo = []
for val in sorted(geo.keys(),key=str.lower):
    sortedGeo.append(val)
i = 1
for s in sortedGeo:
    sw = '%03d %-40s  %lf\n' % (i,s,geo[s])
    fp.write(sw)
    i+=1
fp.close()

