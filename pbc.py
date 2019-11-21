import subprocess
import os




def trjconv( prefix_trj='md' , format_trj='xtc' , prefix_tpr='md' , format_tpr='tpr' , prefix_output = 'md_output' , format_output = 'xtc' , replica=1 , commands=[] , input_string='System' ):
	command=(['gmx','trjconv'])
	
	command.extend([ '-f' , prefix_trj + '.' + str(replica) + '.' + format_trj ])
	command.extend([ '-s' , prefix_tpr + '.' + str(replica) + '.' + format_tpr ])
	
	command.extend( commands )
	
	command.extend([ '-o' , prefix_output + '.' + str(replica) + '.' + format_output ])

	p = subprocess.Popen( command , stdin=subprocess.PIPE , stdout=devnull, stderr=devnull )
	p.communicate( input=input_string.encode('utf-8') )



def tprconv( prefix_input='md', prefix_output='md_out' , replica='1' , input_string='Non-water' ):
	command=(['gmx','convert-tpr'])
	
	command.extend([ '-s' , prefix_input + '.' + str(replica) + '.tpr' ])
	command.extend([ '-o' , prefix_output + '.' + str(replica) + '.tpr' ])

	p = subprocess.Popen( command , stdin=subprocess.PIPE , stdout=devnull, stderr=devnull )
	p.communicate( input=input_string.encode('utf-8') )



def treat_pbc( replica=1 , grp_center='Protein' , grp_fit='Protein' , prefix_trj='md' , format_trj='xtc' , prefix_tpr='md' , format_tpr='tpr' , verbose=True ):
	if verbose:
		print('\n\t\tMaking all molecules whole ... ', end='', flush=True)
	trjconv(
		prefix_trj    = prefix_trj ,
		format_trj    = format_trj ,
		prefix_tpr    = prefix_tpr ,
		format_tpr    = format_tpr ,
		prefix_output = prefix_trj+'_whole' ,
		format_output = 'xtc' ,
		replica       = replica ,
		commands      = ['-pbc','whole'] ,
		input_string  = 'System'
	)
	if verbose:
		print('done')
		
		
		
	if verbose:
		print('\t\tExtracting the first frame ... ', end='', flush=True)
	trjconv(
		prefix_trj    = prefix_trj+'_whole',
		format_trj    = 'xtc',
		prefix_tpr    = prefix_tpr,
		format_tpr    = 'tpr',
		prefix_output = prefix_trj+'_whole_t0',
		format_output = 'gro',
		replica       = replica,
		commands      = ['-e','0'],
		input_string  = 'System'
	)
	if verbose:
		print('done')
		
		
		
	if verbose:
		print('\t\tRemoving jumps ... ', end='', flush=True)
	trjconv(
		prefix_trj    = prefix_trj+'_whole',
		format_trj    = 'xtc',
		prefix_tpr    = prefix_trj+'_whole_t0',
		format_tpr    = 'gro',
		prefix_output = prefix_trj+'_nojump',
		format_output = 'xtc',
		replica       = replica,
		commands      = ['-pbc','nojump'],
		input_string  = 'System'
	)
	if verbose:
		print('done')
		
		
		
	if verbose:
		print('\t\tCentering on group %s ... ' % grp_center, end='', flush=True)
	trjconv(
		prefix_trj    = prefix_trj+'_nojump',
		format_trj    = 'xtc',
		prefix_tpr    = prefix_tpr,
		format_tpr    = 'tpr',
		prefix_output = prefix_trj+'_center',
		format_output = 'xtc',
		replica       = replica,
		commands      = ['-center'],
		input_string  = grp_center+' System'
	)
	if verbose:
		print('done')
		
		
		
	if verbose:
		print('\t\tReconstructing a compact box ... ', end='', flush=True)
	trjconv(
		prefix_trj    = prefix_trj+'_center',
		format_trj    = 'xtc',
		prefix_tpr    = prefix_tpr,
		format_tpr    = 'tpr',
		prefix_output = prefix_trj+'_pbc',
		format_output = 'xtc',
		replica       = replica,
		commands      = ['-pbc','mol','-ur','compact'],
		input_string  = 'System'
	)
	if verbose:
		print('done')
		
		
		
	if verbose:
		print('\t\tFitting to group %s ... ' % grp_fit , end='', flush=True)
	trjconv(
		prefix_trj    = prefix_trj+'_pbc',
		format_trj    = 'xtc',
		prefix_tpr    = prefix_tpr,
		format_tpr    = 'tpr',
		prefix_output = prefix_trj+'_fit',
		format_output = 'xtc',
		replica       = replica,
		commands      = ['-fit','rot+trans'],
		input_string  = grp_fit+' System'
	)
	if verbose:
		print('done')
		
		
		
	if verbose:
		print('\t\tCreating non-water trajectory ... ', end='', flush=True)
	trjconv(
		prefix_trj    = prefix_trj+'_fit',
		format_trj    = 'xtc',
		prefix_tpr    = prefix_tpr,
		format_tpr    = 'tpr',
		prefix_output = prefix_trj+'_non-water',
		format_output = 'xtc',
		replica       = replica,
		commands      = [],
		input_string  = 'Non-water'
	)
	if verbose:
		print('done')
		
		
	if verbose:
		print('\t\tCreating non-water tpr ... ', end='', flush=True)
	tprconv(
		prefix_input  = 'md',
		prefix_output = prefix_trj+'_non-water',
		replica       = replica,
		input_string  = 'Non-water'
	)
	if verbose:
		print('done')
		
		
		
	if verbose:
		print('\t\tCreating non-water trajectory for visualization ... ', end='', flush=True)
	trjconv(
		prefix_trj    = prefix_trj+'_non-water',
		format_trj    = 'xtc',
		prefix_tpr    = prefix_tpr+'_non-water',
		format_tpr    = 'tpr',
		prefix_output = prefix_trj+'_non-water_dt1ns',
		format_output = 'xtc',
		replica       = replica,
		commands      = ['-dt','1000'],
		input_string  = 'System'
	)
	if verbose:
		print('done')
		
		
		
	if verbose:
		print('\t\tExtracting the first frame of the non-water visualization trajectory ... ', end='', flush=True)
	trjconv(
		prefix_trj    = prefix_trj+'_non-water',
		format_trj    = 'xtc',
		prefix_tpr    = prefix_tpr+'_non-water',
		format_tpr    = 'tpr',
		prefix_output = prefix_trj+'_non-water_dt1ns',
		format_output = 'gro',
		replica       = replica,
		commands      = ['-e','0'],
		input_string  = 'System'
	)
	if verbose:
		print('done')



	os.remove( prefix_trj+'_whole.'    + str(replica) + '.xtc' )
	os.remove( prefix_trj+'_whole_t0.' + str(replica) + '.gro' )
	os.remove( prefix_trj+'_nojump.'   + str(replica) + '.xtc' )
	os.remove( prefix_trj+'_center.'   + str(replica) + '.xtc' )
	os.remove( prefix_trj+'_pbc.'      + str(replica) + '.xtc' )



def get_replicas():
	replicas=[int(file.split('.')[-2]) for file in os.listdir() if file.split('.')[-1]=='tpr' and file.split('.')[0]=='md']
	replicas.sort()
	
	return replicas





# ----------------------------------------------------------------------------

devnull = open(os.devnull, 'w')

replicas = get_replicas()

print("Removing PBC artifacts for %d trajectories" % len(replicas) )

for replica in replicas:
	print("\tReplica %d ... " % int(replica) , end='' )
	treat_pbc( replica=replica )
	print('\tOK\n')

devnull.close()
