#!/usr/bin/python

import os,sys


filename = sys.argv[1]

print 'version	0.1'

numVoxels = int(os.popen('wc -l '+filename).read().split()[0])
print '# total numVoxels:',numVoxels

print '# Assuming cubic'
sideLen = int(round(pow(numVoxels,(1.0/3.0))))
print 'numVoxels',sideLen,sideLen,sideLen

infile = open(filename)

firstline = infile.readline()
secondline = infile.readline()


voxelSize = abs(float(firstline.split()[2]) - float(secondline.split()[2]))
print 'voxelSize',voxelSize,voxelSize,voxelSize

print 'lowBounds',' '.join(firstline.split()[0:3])

print 'colloidPos', '1', -11.3795,0,0,0,'# First number is index'
print 'colloidPos', '2',  11.3795,0,0,0,'# First number is index'

print 'beginClCpCs zyxInc'
print ' '.join(firstline.split()[3:])
print ' '.join(secondline.split()[3:])
line=True
while line:
	line=infile.readline()
	if line: print ' '.join(line.split()[3:])
