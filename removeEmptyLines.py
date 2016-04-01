import os
import fnmatch

cwd = os.getcwd()
outputsPath = cwd + "/output2/"
os.chdir(outputsPath)
currStep = -1
totalNumFiles = 0
fileNameList = []
allFiles = os.listdir(outputsPath)

for fileName in allFiles:
    f = open(fileName, 'r')
    lines = f.readlines()
    f.close()
    print fileName
    f = open(fileName, 'w')
    for line in lines:
        if(line[0].isdigit() | (line[0] == 'x')):
            f.write(line)
    f.close()


