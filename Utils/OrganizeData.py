import json
from pprint import pprint

def s(i):
	return str(i)

def getSimID(numNodes, N, csm, h):
	n = numNodes * N
	SimID = 'n'+s(n)+'_N'+s(N)+'_csm'+s(csm)+'_h'+s(h) 
	return SimID 

numNodes = 1
h = 0.025
N = [1, 2, 4]
csm = [2, 4, 8]

TimingFilename = 'Timing_100.json'

SimID = getSimID(numNodes, N[0],csm[1],h)


with open('../output/'+SimID+'/Timing_100.json') as data_file:    
    data = json.load(data_file)

pprint(data)