import os
import csv
from array import array

cwd = os.getcwd()
outputsPath = cwd + "/../output/"
os.chdir(outputsPath)

print os.getcwd()

xFluidMarker = array('f')
yFluidMarker = array('f')
zFluidMarker = array('f')
xBoundaryMarker = array('f')
yBoundaryMarker = array('f')
zBoundaryMarker = array('f')

errorCounter = 0

with open('step.0.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        try:
            markerType = int(row[len(row)-1])
            if(markerType == -1):
                xFluidMarker.append(float(row[0]))
                yFluidMarker.append(float(row[1]))
                zFluidMarker.append(float(row[2]))
            elif(markerType == 0):
                xBoundaryMarker.append(float(row[0]))
                yBoundaryMarker.append(float(row[1]))
                zBoundaryMarker.append(float(row[2]))               
        except ValueError, e:
            errorCounter += 1
            print "Error",e,". errorCounter=",errorCounter

print "Number of Fluid markers = ", len(xFluidMarker)
print "Number of Boundary markers = ", len(xBoundaryMarker)

counter = 0
# Note that range stops at len(x)-1 not len(x)
for fluidIndex in range(0,len(xFluidMarker)):
    fluidX = xFluidMarker[fluidIndex]
    fluidY = yFluidMarker[fluidIndex]
    fluidZ = zFluidMarker[fluidIndex]
    for boundaryIndex in range(0,len(xBoundaryMarker)):
        boundaryX = xBoundaryMarker[boundaryIndex]
        boundaryY = yBoundaryMarker[boundaryIndex]
        boundaryZ = zBoundaryMarker[boundaryIndex]

        if((fluidX == boundaryX) & (fluidY == boundaryY) & (fluidZ == boundaryZ)):
            print "Equal fluid and boundary!"

        diffX = abs(fluidX - boundaryX) + 0.000001 
        diffY = abs(fluidY - boundaryY) + 0.000001
        diffZ = abs(fluidZ - boundaryZ) + 0.000001
        if((diffX < 0.05) & (diffY < 0.05) & (diffZ < 0.05)):
            print "Equal vals - ", diffX, ", ", diffY, ", ", diffZ, ", "
            counter +=1
            #print counter













