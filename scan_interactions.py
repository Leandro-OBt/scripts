import subprocess
import multiprocessing as mp
import os
import numpy as np
import sys


def parse_ndx( target_grp='Protein' , input_file='residues.ndx' ):
	#      This function parses the groups of an input index file into a dictionary of residue names
	#      The input index file must contain:
	#           - The group of interest (target_grp)
	#           - Its decomposition into their residues
	#           - The decomposition must come AFTER the group itself
	
	
	# Creates the dictionary that will hold all residue names
	res_dict = {}
	
	with open( input_file , 'r') as index_file:
		lines = index_file.readlines()
		for line in lines:
			if( line[0]=='['):	# If the first character of this line is '[', we have a group that will be parsed
				
				words = line.strip().replace('[ ','').replace(' ]','').split('_')
				# words is now a list in the format [GROUP,RESN,RESI] or [GROUP]
				
				grp = words[0]	# The first element of words is always the group name
				
				if( grp==target_grp ):
					
					if( len(words)==3 ):	# If len(words)!=1, we have [GROUP,RESN,RESI]
						res_name   =      words[1] #.capitalize()
						res_number = int( words[2] )
						
						res_dict[res_number] = res_name
				
	return res_dict

def find_contacting_residues( grp_decomp='Protein_chain1', grp_whole='Protein_chain2' , rep=1 , cutoff=0.4, first_res_decomp=1, input_path='../../' , dt='1000', b='25000', ndx='residues.ndx' , input_prefix='md_non-water'):
	output_file = 'mindist_res-from-'+str(grp_decomp)+'_-_'+str(grp_whole)+'.'+str(rep)+'.xvg'
	tmp_file = 'temporary_mindist_res-from-'+str(grp_decomp)+'_-_'+str(grp_whole)+'.'+str(rep)+'.xvg'
	
	command = ['gmx','mindist',
		'-f',input_path + input_prefix + '.' + str(rep) + '.xtc',
		'-s',input_path + input_prefix + '.' + str(rep) + '.tpr',
		'-n',ndx,
		'-dt',str(dt),
		'-b',str(b),
		'-xvg','none',
		'-or',output_file,
		'-od',tmp_file]
	
	input_string = str(grp_decomp) + ' ' + str(grp_whole)
	
	p = subprocess.run( command , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL , input=input_string.encode('utf-8'))
	os.remove(tmp_file)
	
	res_list = []
	
	with open(output_file,'r') as mindist:
		lines = mindist.readlines()
		for line in lines:
			line = [x for x in line.strip().split(' ') if x]
			if( float(line[1]) <= cutoff ):
				res_list.append( int(line[0]) + (first_res_decomp-1) )	# This is to correct the residue numbering. mindist numbers the residues starting from 1
	return( res_list )

def hbond( group1='Protein', group2='Protein' , rep=1 , input_path='../../' , dt='100', b='25000', ndx='residues.ndx' , input_prefix='md_non-water', threads=1):
	out_filename = 'hbonds_'+str(group1)+'_-_'+str(group2)+'.'+str(rep)+'.xvg'
	
	command = ['gmx','hbond',
		'-f',input_path + input_prefix + '.' + str(rep) + '.xtc',
		'-s',input_path + input_prefix + '.' + str(rep) + '.tpr',
		'-n',ndx,
		'-dt',str(dt),
		'-b',str(b),
		'-xvg','none',
		'-nthreads',str(threads),
		'-num',out_filename]
	
	input_string = str(group1) + ' ' + str(group2)
	p = subprocess.run( command , stdout=subprocess.PIPE, stderr=subprocess.PIPE , input=input_string.encode('utf-8'))
	
	
	hb=[]
	try:
		with open(out_filename,'r') as tmp:
			hb = [ int(line.strip().split()[1]) for line in tmp.readlines() ]
		
	except:
		if( ('Nothing to be done' in p.stderr.decode()) and ('No Donors found' in p.stdout.decode()) ):
			hb = np.zeros((5,), dtype=int)
			with open(out_filename,'w') as tmp:
				tmp.write('0 0 0\n0 0 0\n')
		else:
			print('-- ERROR --')
			print('command =',command)
			print('input_string =', input_string)
			print('stderr = ',p.stderr.decode())
			print('stdout = ',p.stdout.decode())
			print('-- ERROR --')
	
	hb  = np.array(hb)
	occ = hb > 0
	
	hb_num_avg = hb.mean()
	hb_num_std = hb.std(ddof=1)
	hb_occ_avg = occ.mean()
	hb_occ_std = occ.std(ddof=1)
	
	return( group1 , group2 , hb_num_avg , hb_num_std , hb_occ_avg, hb_occ_std)

def other_group( group_list , group ):
	return [x for x in group_list if x!=group][0]

#'GLY','SER','ASP','GLU','ASN','GLN','LYS'

area = {
'ALA':['CB','HB1','HB2','HB3'],
'THR':['CG2','1HG2','2HG2','3HG2'],
'CYS':['HG','SG','CB','HB1','HB2'],
'VAL':['CB','HB','CG1','1HG1','2HG1','3HG1','CG2','1HG2','2HG2','3HG2'],
'LEU':['CB','HB1','HB2','CG','HG','CD1','1HD1','2HD1','3HD1','CD2','1HD2','2HD2','3HD2'],
'ILE':['CB','HB','CG1','1HG1','2HG1','CD','HD1','HD2','HD3','CG2','1HG2','2HG2','3HG2'],
'MET':['CB','HB1','HB2','CG','HG1','HG2','SD','CE','HE1','HE2','HE3'],
'PRO':['CA','HA','CB','HB1','HB2','CG','HG1','HG2','CD','HD1','HD2','N'],
'PHE':['CB','HB1','HB2','CG','CD1','HD1','CE1','HE1','CZ','HZ','CE2','HE2','CD2','HD2'],
'TYR':['CB','HB1','HB2','CG','CD1','HD1','CE1','HE1','CZ','CE2','HE2','CD2','HD2'],
'TRP':['CB','HB1','HB2','CG','CD2','CE3','HE3','CZ3','HZ3','CH2','HH2','CZ2','HZ2','CE2','NE1','HE1','CD1','HD1']
}

stacking = {
'ARG':['NE','HE','CZ','NH1','1HH1','2HH1','NH2','2HH1','2HH2'],
'HIS':['CG','ND1','HD1','CE1','HE1','NE2','HE2','CD2','HD2'],
'PHE':['CG','CD1','HD1','CE1','HE1','CZ','HZ','CE2','HE2','CD2','HD2'],
'TYR':['CG','CD1','HD1','CE1','HE1','CZ','CE2','HE2','CD2','HD2'],
'TRP':['CG','CD2','CE3','HE3','CZ3','HZ3','CH2','HH2','CZ2','HZ2','CE2','NE1','HE1','CD1','HD1']
}

# List with the groups of interest
grps = []
#grps.append( 'ChainA' )
#grps.append( 'ChainB' )
grps.append( sys.argv[1] )
grps.append( sys.argv[2] )


num_proc = 12
cutoff   = 0.4

reps = [1,2,3,4,5]
reps_total = len(reps)

if __name__=='__main__':
	residues  = {}
	res_first = {}
	res_last  = {}
	
	dict_hbs_byres = {}
	dict_hbs_pairs = {}
	dict_hbs_total = {}
	
	for grp in grps:
		residues[grp] = parse_ndx(target_grp=grp , input_file='residues.ndx')
		
		res_first[grp] = sorted( list(residues[grp].keys()) )[0] 
		res_last[grp]  = sorted( list(residues[grp].keys()) )[-1]
	
	for rep in reps:
		# grp0-grp1 hydrogen bonds
		hbs = hbond( group1=grps[0], group2=grps[1] , rep=rep , input_path='../../' , dt='100', b='25000', ndx='residues.ndx' , input_prefix='md_non-water', threads=num_proc)
		
		if( grps[0] not in dict_hbs_total ):
			dict_hbs_total[ grps[0] ] = {}
		if( grps[1] not in dict_hbs_total[ grps[0] ] ):
			dict_hbs_total[ grps[0] ][ grps[1] ] = {}
		
		if( 'num' not in dict_hbs_total[ grps[0] ][ grps[1] ] ):
			dict_hbs_total[ grps[0] ][ grps[1] ]['num'] = []
			dict_hbs_total[ grps[0] ][ grps[1] ]['occ'] = []
			
		group1_whole   = hbs[0]
		group2_whole   = hbs[1]
		hb_num_avg     = hbs[2]
		hb_num_std     = hbs[3]
		hb_occ_avg     = hbs[4]
		hb_occ_std     = hbs[5]
		
		dict_hbs_total[ group1_whole ][ group2_whole ]['num'].append( hb_num_avg )
		dict_hbs_total[ group1_whole ][ group2_whole ]['occ'].append( hb_occ_avg )
		
		
		
		
		
		# Hydrogen bonds decomposed by residue for each group
		# finding the interacting residues in each group
		contacting_residues = {}
		
		for grp in grps:
			contacting_residues[grp] = find_contacting_residues( grp_decomp=grp , grp_whole=other_group(grps,grp) , rep=rep , cutoff=cutoff, first_res_decomp=res_first[grp] , input_path='../../' )
			param_list=[]
			
			for res in contacting_residues[grp]:
				
				group1 = grp+'_'+residues[grp][res]+'_'+str(res)
				group2 = other_group(grps,grp)
				
				param_list.append( (group1 , group2 , rep , '../../' , '100', '25000', 'residues.ndx' , 'md_non-water', 1) )
				
			pool = mp.Pool(num_proc)
			hbs = pool.starmap( hbond , param_list )
			pool.close()
			pool.join()
			
			for hb in hbs:
				group1_residue = hb[0]
				group2_whole   = hb[1]
				hb_num_avg     = hb[2]
				hb_num_std     = hb[3]
				hb_occ_avg     = hb[4]
				hb_occ_std     = hb[5]
				
				if( group1_residue not in dict_hbs_byres ):
					dict_hbs_byres[group1_residue] = {}
				
				if( group2_whole not in dict_hbs_byres[group1_residue] ):
					dict_hbs_byres[group1_residue][group2_whole] = {}
					dict_hbs_byres[group1_residue][group2_whole]['num'] = []
					dict_hbs_byres[group1_residue][group2_whole]['occ'] = []
					
				dict_hbs_byres[group1_residue][group2_whole]['num'].append( hb_num_avg )
				dict_hbs_byres[group1_residue][group2_whole]['occ'].append( hb_occ_avg )
				
		
		
		
		
		
		
		
		# which residues from grps[1] interact with each residue from grps[0] ?
		param_list = []
		for res in contacting_residues[ grps[0] ]:
			grp_decomp = other_group(grps,grps[0])
			grp_whole  = grps[0]+'_'+residues[ grps[0] ][res]+'_'+str(res)
			
			param_list.append( (grp_decomp,grp_whole,rep,cutoff,res_first[other_group(grps,grps[0])],'../../','1000','25000','residues.ndx','md_non-water') )
		
		pool = mp.Pool(num_proc)
		contacting_residues_pairs = pool.starmap( find_contacting_residues , param_list )
		pool.close()
		pool.join()
		
		contacting_pairs = []
		i=0
		for res_grp0 in contacting_residues[ grps[0] ]:
			for res_grp1 in contacting_residues_pairs[i]:
				contacting_pairs.append( (int(res_grp0) , int(res_grp1)) )
			i+=1
		
		param_list = []
		for pair in contacting_pairs:
			res_number_grp0 = pair[0]
			res_number_grp1 = pair[1]
			res_name_grp0   = residues[grps[0]][res_number_grp0]
			res_name_grp1   = residues[grps[1]][res_number_grp1]
			
			input_grp1 = grps[0]+'_'+res_name_grp0+'_'+str(res_number_grp0)
			input_grp2 = grps[1]+'_'+res_name_grp1+'_'+str(res_number_grp1)
			
			param_list.append( (input_grp1,input_grp2,rep,'../../','100','25000','residues.ndx','md_non-water') )
		
		pool = mp.Pool(num_proc)
		hbs = pool.starmap( hbond , param_list )
		pool.close()
		pool.join()
	
		for hb in hbs:
			words1 = hb[0].split('_')
			words2 = hb[1].split('_')
			
			#grp1            = words1[0]
			#grp1_res_name   = words1[1]
			#grp1_res_number = words1[2]
	
			#grp2            = words2[0]
			#grp2_res_name   = words2[1]
			#grp2_res_number = words2[2]
			
			grp1    = hb[0]
			grp2    = hb[1]
			num_avg = hb[2]
			num_std = hb[3]
			occ_avg = hb[4]
			occ_std = hb[5]
			
#			if( grp1_res_number not in dict_hbs_pairs ):
#				dict_hbs_pairs[ grp1_res_number ] = {}
#				
#			if( grp2_res_number not in dict_hbs_pairs[grp1_res_number] ):
#				dict_hbs_pairs[grp1_res_number][grp2_res_number] = {}
#				dict_hbs_pairs[grp1_res_number][grp2_res_number]['num'] = []
#				dict_hbs_pairs[grp1_res_number][grp2_res_number]['occ'] = []
#			
#			dict_hbs_pairs[grp1_res_number][grp2_res_number]['num'].append( num_avg )
#			dict_hbs_pairs[grp1_res_number][grp2_res_number]['occ'].append( occ_avg )
	
			if( grp1 not in dict_hbs_pairs ):
				dict_hbs_pairs[grp1] = {}
				
			if( grp2 not in dict_hbs_pairs[grp1] ):
				dict_hbs_pairs[grp1][grp2] = {}
				dict_hbs_pairs[grp1][grp2]['num'] = []
				dict_hbs_pairs[grp1][grp2]['occ'] = []
			
			dict_hbs_pairs[grp1][grp2]['num'].append( num_avg )
			dict_hbs_pairs[grp1][grp2]['occ'].append( occ_avg )
	
	
	
	
	list_hbs_total = []
	for grp1 in dict_hbs_total:
		for grp2 in dict_hbs_total[grp1]:
			
			hb_num = dict_hbs_total[grp1][grp2]['num']
			reps_missing = reps_total - len(hb_num)
			hb_num = np.pad( hb_num , (0,reps_missing) , 'constant' )
			
			hb_occ = dict_hbs_total[grp1][grp2]['occ']
			reps_missing = reps_total - len(hb_occ)
			hb_occ = np.pad( hb_occ , (0,reps_missing) , 'constant' )
			
			dict_hbs_total[grp1][grp2]['num_avg'] = hb_num.mean()
			dict_hbs_total[grp1][grp2]['num_err'] = hb_num.std(ddof=0)
			
			dict_hbs_total[grp1][grp2]['occ_avg'] = hb_occ.mean()
			dict_hbs_total[grp1][grp2]['occ_err'] = hb_occ.std(ddof=0)
			
			list_hbs_total.append( (grp1 , grp2 , dict_hbs_total[grp1][grp2]['num_avg'], dict_hbs_total[grp1][grp2]['num_err'] , dict_hbs_total[grp1][grp2]['occ_avg'] , dict_hbs_total[grp1][grp2]['occ_err']) )
	
	
	
	
	
	list_hbs_byres = [ [] , [] ]
	
#	list_hbs_byres_grp0 = []
#	list_hbs_byres_grp1 = []
	
	for grp_decomp in dict_hbs_byres:
		for grp_whole in dict_hbs_byres[grp_decomp]:
#			print(grp_decomp,grp_whole, dict_hbs_byres[grp_decomp][grp_whole] )
			
			hb_num = dict_hbs_byres[grp_decomp][grp_whole]['num']
			reps_missing = reps_total - len(hb_num)
			hb_num = np.pad( hb_num , (0,reps_missing) , 'constant' )
			
			hb_occ = dict_hbs_byres[grp_decomp][grp_whole]['occ']
			reps_missing = reps_total - len(hb_occ)
			hb_occ = np.pad( hb_occ , (0,reps_missing) , 'constant' )
			
			dict_hbs_byres[grp_decomp][grp_whole]['num_avg'] = hb_num.mean()
			dict_hbs_byres[grp_decomp][grp_whole]['num_err'] = hb_num.std(ddof=0)
			
			dict_hbs_byres[grp_decomp][grp_whole]['occ_avg'] = hb_occ.mean()
			dict_hbs_byres[grp_decomp][grp_whole]['occ_err'] = hb_occ.std(ddof=0)
			
			hbs_data = (grp_decomp , grp_whole , dict_hbs_byres[grp_decomp][grp_whole]['num_avg'], dict_hbs_byres[grp_decomp][grp_whole]['num_err'] , dict_hbs_byres[grp_decomp][grp_whole]['occ_avg'] , dict_hbs_byres[grp_decomp][grp_whole]['occ_err'])
			
			if( grp_whole == grps[0] ):
				list_hbs_byres[1].append( hbs_data )
			else:
				list_hbs_byres[0].append( hbs_data )
			
			
	for i in range(0,2):
		list_hbs_byres[i].sort(key=lambda tup: tup[2] , reverse=True)
	
	
	
	
	list_hbs_pairs = []
	for grp1 in dict_hbs_pairs:
		for grp2 in dict_hbs_pairs[grp1]:
			
			hb_num = dict_hbs_pairs[grp1][grp2]['num']
			reps_missing = reps_total - len(hb_num)
			hb_num = np.pad( hb_num , (0,reps_missing) , 'constant' )
			
			hb_occ = dict_hbs_pairs[grp1][grp2]['occ']
			reps_missing = reps_total - len(hb_occ)
			hb_occ = np.pad( hb_occ , (0,reps_missing) , 'constant' )
			
			dict_hbs_pairs[grp1][grp2]['num_avg'] = hb_num.mean()
			dict_hbs_pairs[grp1][grp2]['num_err'] = hb_num.std(ddof=0)
			
			dict_hbs_pairs[grp1][grp2]['occ_avg'] = hb_occ.mean()
			dict_hbs_pairs[grp1][grp2]['occ_err'] = hb_occ.std(ddof=0)
			
			list_hbs_pairs.append( (grp1,grp2, float(dict_hbs_pairs[grp1][grp2]['num_avg']), dict_hbs_pairs[grp1][grp2]['num_err'], dict_hbs_pairs[grp1][grp2]['occ_avg'], dict_hbs_pairs[grp1][grp2]['occ_err']) )
	list_hbs_pairs.sort(key=lambda tup: tup[2] , reverse=True)
	
	
	
	
	
	
	print('OVERALL HBs')
	with open('hbs_'+str(grps[0])+'-'+str(grps[1])+'_overall.dat','w') as output_file:
		output_file.write('%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\n' % list_hbs_total[0] )
	print('%s %s\t%.2f +- %.2f\t%.2f +- %.2f' % list_hbs_total[0] )
	print()
	print()
	
	
	print('HBs BY RESIDUE')
	for i in range(0,2):
		averages = []
		with open('hbs_'+str(grps[0])+'-'+str(grps[1])+'_byres_'+str(grps[i])+'.dat','w') as output_file:
			for hb in list_hbs_byres[i]:
				output_file.write('%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\n' % hb )
				if( hb[4] >= 0.01 ):
					print('%s %s\t%.2f +- %.2f\t%.2f +- %.2f' % hb )
				averages.append( hb[2] )
			total_avg = np.sum(averages)
			total_err = np.sqrt(np.sum(np.power(averages,2)))
			print('Total = %.2f +- %.2f' % (total_avg,total_err))
			output_file.write('Total\t%.6f\t%.6f\n' % (total_avg,total_err))
		print()
	
#	print()
	
	
	print('HBs BY PAIRS')
	averages = []
	with open('hbs_'+str(grps[0])+'-'+str(grps[1])+'_pairs.dat','w') as output_file:
		for hb in list_hbs_pairs:
			output_file.write('%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\n' % hb )
			if( hb[4] >= 0.01 ):
				print('%s %s\t%.2f +- %.2f\t%.2f +- %.2f' % hb )
			averages.append( hb[2] )
		total_avg = np.sum(averages)
		total_err = np.sqrt(np.sum(np.power(averages,2)))
		print('Total = %.2f +- %.2f' % (total_avg,total_err))
		output_file.write('Total\t%.6f\t%.6f\n' % (total_avg,total_err))
	print()
	print()
	
