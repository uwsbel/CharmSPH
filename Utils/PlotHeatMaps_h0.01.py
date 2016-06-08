import json
from pprint import pprint
import PlotHeatMap

def s(i):
	return str(i)

def getSimID(numNodes, N, csm, h):
	n = numNodes * N
	SimID = 'n'+s(n)+'_N'+s(N)+'_csm'+s(csm)+'_h'+s(h) 
	return SimID 

numNodes = [8]
h = 0.01
csm = [2, 4, 8]
csm_str = ['2h', '4h', '8h']
num_csm = 3
N = [1,2,4]
N_str = ['1 PE','2 PE','4 PE']
num_N = 3

# base_times = [0 for col in range(num_csm)]
# for i in range(0,num_csm):
# 	curr_numNodes = 1
# 	curr_N = 1
# 	curr_csm = csm[i]
# 	curr_id = getSimID(curr_numNodes, curr_N, curr_csm, h)
# 	with open('../h0.01/'+curr_id+'/Timing_100.json') as data_file: 
# 		curr_data = json.load(data_file)
# 		base_times[i] = curr_data['AvgTimePerStep'] / 1000




for k in range(0,1):
	curr_numNodes = numNodes[k]
	avgTimeData = [[0 for col in range(num_csm)] for row in range(num_N)]
	speedupData = [[0 for col in range(num_csm)] for row in range(num_N)]
	simIds = [['' for col in range(num_csm)] for row in range(num_N)]
	for i in range(0,num_N):
		curr_N = N[i]
		for j in range(0,num_csm):
			curr_csm = csm[j]
			curr_id = getSimID(curr_numNodes, curr_N, curr_csm, h)
			try:
				with open('../h0.01/'+curr_id+'/Timing_40.json') as data_file: 
					curr_data = json.load(data_file)
					# pprint(curr_data)
					simIds[i][j] = curr_id
					avgTimeData[i][j] = round(curr_data['AvgTimePerStep'] / 1000, 4)
					# speedupData[i][j] = round(base_times[j] / avgTimeData[i][j], 4)
			except FileNotFoundError as e:
				simIds[i][j] = curr_id
				avgTimeData[i][j] = round(100, 4)			

	plotTitle = 'Avg. Time Per Step (NumNodes: '+str(curr_numNodes)+', h: '+str(h) + ')'
	PlotHeatMap.PlotHeatMap(avgTimeData, csm_str, N_str, plotTitle)
	# plotTitle = 'Speedup Per Step (NumNodes: '+str(curr_numNodes)+', h: '+str(h) +')'
	# PlotHeatMap.PlotHeatMap(speedupData, csm_str, N_str, plotTitle)
