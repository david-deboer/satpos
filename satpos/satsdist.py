from argparse import Namespace
import numpy as np

view = Namespace(iarr=[], period=[], zmin=[], zmax=[])
notv = Namespace(iarr=[], period=[], zmin=[], zmax=[])
with open('viewable.out', 'r') as fp:
    for line in fp:
        data = line.split('|')
        view.iarr.append(int(data[1].split()[0]))
        view.period.append(float(data[1].split('=')[1].split()[0]))
        view.zmin.append(float(data[2]))
        view.zmax.append(float(data[3]))
view.iarr = np.array(view.iarr)
view.period = np.array(view.period)
view.zmin = np.array(view.zmin)
view.zmax = np.array(view.zmax)

with open('notviewable.out', 'r') as fp:
    for line in fp:
        data = line.split('|')
        notv.iarr.append(int(data[1].split()[0]))
        notv.period.append(float(data[1].split('=')[1].split()[0]))
        try:
            notv.zmin.append(float(data[2]))
            notv.zmax.append(float(data[3]))
        except ValueError:
            notv.zmin.append(0.0)
            notv.zmax.append(0.0)
notv.iarr = np.array(notv.iarr)
notv.period = np.array(notv.period)
notv.zmin = np.array(notv.zmin)
notv.zmax = np.array(notv.zmax)
