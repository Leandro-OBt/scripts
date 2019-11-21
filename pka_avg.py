import numpy as np
import pandas as pd
import os

reps = [ 1 , 2 , 3 , 4 , 5 ]
#reps = [ 1 ]


pwd = os.getcwd()

pkas = {}

for rep in reps:
	allfiles = os.listdir(pwd+'/'+str(rep))
	path = pwd + '/' + str(rep) + '/'
	
	
	for filename in allfiles:
		if( filename.split('.')[-1] == 'xvg' ):
			fullpath = path + filename
			
			filename = filename.replace('-','_').replace('.','_')
			filename = filename.split('_')
			
			res_name = str(filename[1])
			res_numb = str(filename[2])
			chain    = str(filename[3])
			
			residue  =  res_name + '-' + res_numb + '_' + chain
			
			if( residue not in pkas ):
				pkas[residue] = []
			
			xvg         = pd.read_csv( fullpath , sep='\t' , header=None)
			xvg.columns = ['time','pka']
			
			avg = xvg['pka'].mean()
			std = xvg['pka'].std(ddof=1)
			
			pkas[residue].append( avg )
			

residue_list = [x for x in pkas.keys()]
residue_list.sort()


with open('pka_avg.dat','w') as output_file:
	for residue in residue_list:
		vec = np.array( pkas[residue] )
		
		avg = vec.mean()
		err = vec.std() / np.sqrt( len(vec) )
		
		output_file.write('%s\t%.6f\t%.6f\n' % (residue , avg, err))
	
