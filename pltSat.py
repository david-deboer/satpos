import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

cursublng = 1000.0
cursublat = 1000.0
try:
    fp = open('info.out', 'r')
    line = fp.readline()
    data = line.split('|')
    obsName = data[0]
    scName = data[1]
    line = fp.readline()
    data = line.split()
    obsLng = float(data[0])
    obsLat = float(data[1])
    line = fp.readline()
    if line[0] == '=':
        data = line.strip('=').split()
        cursublng = float(data[0])
        cursublat = float(data[1])
        print('Current:', cursublng, cursublat)
except IOError:
    obsName = 'obs'
    scName = 'sc'
    obsLng = 0.0
    obsLat = 0.0

predix = True
try:
    fp = open('predix.out', 'r')
    x = []
    y1 = []
    y2 = []
    for line in fp:
        data = line.split()
        x.append(float(data[0]))
        y1.append(float(data[2]))
        y2.append(float(data[3]))
    fp.close()
except IOError:
    print('predix.out not found')
    predix = False

subsat = True
slong = []
slat = []
try:
    fp = open('subsat.out', 'r')
    slong = []
    slat = []
    for line in fp:
        data = line.split()
        slong.append(float(data[0]))
        slat.append(float(data[1]))
    fp.close()
except IOError:
    print('subsat.out not found')
    subsat = False

if predix:
    plt.figure(1)
    plt.subplot(121)
    plt.plot(x, y1, 'r.', label='Az')
    plt.plot(x, y2, 'b.', label='El')
    plt.xlabel('Time [days]')
    plt.ylabel('Degrees')
    plt.legend()
    plt.title('obs: %s, sc: %s' % (obsName, scName))
    plt.subplot(122)
    plt.plot(y1, y2, '.')
    plt.xlabel('Az [deg]')
    plt.ylabel('El [deg]')

if subsat:
    # ## Get observatory data
    useObs = ['SC', 'HN', 'NL', 'FD', 'LA', 'PT', 'KP', 'OV', 'BR', 'MK', 'UCB']
    obsfp = open('obs.dat', 'r')
    obsLats = []
    obsLongs = []
    obsCodes = []
    for line in obsfp:
        if line[0] == '#' or line[0] == '!':
            continue
        data = line.split()
        if data[0] in useObs:
            # print('Using '+data[1])
            obsCodes.append(data[0])
            tmp = data[2].split(':')
            obsLats.append(float(tmp[0]) + float(tmp[1])/60.0 + float(tmp[2])/3600.0)
            tmp = data[3].split(':')
            obsLongs.append(float(tmp[0]) + float(tmp[1])/60.0 + float(tmp[2])/3600.0)
    obsfp.close()
    # ## Draw map
    plt.figure(2)
    # m = Basemap(width=17000000,height=12000000,projection='lcc',
      #      resolution='c',lat_1=15.,lat_2=55,lat_0=50,lon_0=-100.)
    m = Basemap(lon_0=cursublng)
    m.drawmapboundary(fill_color='aqua')
    m.fillcontinents(color='coral',lake_color='aqua')
    parallels = np.arange(-75.,76,15.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[False,True,True,False])
    meridians = np.arange(int(-180.0+cursublng),int(181.0+cursublng),30.0)
    m.drawmeridians(meridians,labels=[True,False,False,True])
    xpt,ypt = m(obsLongs,obsLats)
    lonpt, latpt = m(xpt,ypt,inverse=True)
    m.plot(xpt,ypt,'bo')
    for i,c in enumerate(obsCodes):
        plt.text(xpt[i]+100000,ypt[i]+100000,'%s' % (c))
    xpt,ypt=m(slong,slat)
    m.plot(xpt,ypt,'-')
    m.drawcountries()
    plt.title(scName)
    if cursublng != 1000.0 and cursublat != 1000.0:
        xpt,ypt = m(cursublng,cursublat)
        print('Plotting current')
        m.plot(xpt,ypt,'ro')

    ffp = open('footprint.out','r')
    lngf = []
    latf = []
    for line in ffp:
        if line[0]=='#' or line[0]=='!':
            continue
        data = line.split()
        lngf.append(float(data[0]))
        latf.append(float(data[1]))
    ffp.close()
    xpt,ypt = m(lngf,latf)
    m.plot(xpt,ypt,'k-')

    plt.show()
