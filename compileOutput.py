import os
import fnmatch

cwd = os.getcwd()
outputsPath = cwd + "/output/fluid/"
os.chdir(outputsPath)
currStep = -1
totalNumFiles = 0
fileNameList = []
allFiles = os.listdir(outputsPath)

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
        os.remove(outputsPath + f)
    firstFile = True


cwd = os.getcwd()
outputsPath = cwd + "/../boundary/"
os.chdir(outputsPath)
currStep = -1
totalNumFiles = 0
fileNameList = []
allFiles = os.listdir(outputsPath)

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
        os.remove(outputsPath + f)
    firstFile = True




