#! /usr/bin/env python
import os

os.system("wget http://www.celestrak.com/NORAD/elements/master.asp -O master.html")

tlefiles = {}

ignore = ['debris']

print 'Reading Celestrak master file'
master = open('master.html','r')
for line in master:
    data = line.split('"')
    for word in data:
        if '.txt' in word:
            description = line.split('>')[1].split('<')[0].strip('(').strip()
            tlefiles[word] = description
            break
master.close()

master = open('master.dat','w')
for lll in tlefiles:
    useThis = True
    for ig in ignore:
        if ig in tlefiles[lll].lower():
            useThis = False
    if useThis:
        a = lll.split('.')
        outfile = a[0]+'.tle'
        s = "wget http://www.celestrak.com/NORAD/elements/%s -O %s" % (lll,outfile)
        print 'Reading %s:  %s' % (lll,tlefiles[lll])
        os.system(s)
        s = "%s:  %s\n" % (outfile,tlefiles[lll])
        master.write(s)
master.close()
