import os
import fnmatch
import sys


print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if(len(sys.argv) == 1):
    SimID = '0.05_4_1_0.0005_1001_1-1-1'
elif(len(sys.argv) == 2):
    SimID = sys.argv[1]
# else:
#     raise InputError()
SimDir = '/output/' + SimID + '/'
BoundaryDir = SimDir + 'boundary/'

cwd = os.getcwd()
FluidDir = cwd + SimDir + 'fluid/'
BoundaryDir = cwd + SimDir + 'boundary/'

print 'FluidDir = ', str(FluidDir)

os.chdir(FluidDir)
currStep = -1
totalNumFiles = 0
fileNameList = []
allFiles = os.listdir(FluidDir)

for fileName in allFiles:
    nameParts = fileName.split('.')
    step = nameParts[1]
    if(currStep != step):
        currStep = step
        currFileName = "fluid." + str(step)
        totalNumFiles += 1
        fileNameList.append(currFileName)
print fileNameList
firstFile = True
for fileName in fileNameList:
    currFiles = fnmatch.filter(allFiles, str(fileName + '.*'))
    print fileName
    with open(str(fileName+'.csv'), "wb") as outfile:
        for chareFile in currFiles:
            with open(chareFile, "rb") as infile:
                header = infile.readline()
                if(firstFile):
                    outfile.write(header)
                    firstFile = False
                outfile.write(infile.read())

    for f in currFiles:
        os.remove(FluidDir + f)
    firstFile = True


cwd = os.getcwd()
os.chdir(BoundaryDir)
currStep = -1
totalNumFiles = 0
fileNameList = []
allFiles = os.listdir(BoundaryDir)

for fileName in allFiles:
    nameParts = fileName.split('.')
    step = nameParts[1]
    if(currStep != step):
        currStep = step
        currFileName = "boundary." + str(step)
        totalNumFiles += 1
        fileNameList.append(currFileName)
print fileNameList
firstFile = True
for fileName in fileNameList:
    currFiles = fnmatch.filter(allFiles, str(fileName + '.*'))
    print fileName
    with open(str(fileName+'.csv'), "wb") as outfile:
        for chareFile in currFiles:
            with open(chareFile, "rb") as infile:
                header = infile.readline()
                if(firstFile):
                    outfile.write(header)
                    firstFile = False
                outfile.write(infile.read())

    for f in currFiles:
        os.remove(BoundaryDir + f)
    firstFile = True




