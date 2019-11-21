import sys
from pymol import cmd

def rmsd_by_residue( pdbid_mobile , pdbid_target ):
	
	alignment_object = pdbid_mobile + '_aligned_to_' + pdbid_target
	
	alignment=cmd.cealign( target=pdbid_target+' and chain A' , mobile=pdbid_mobile + ' and (chain A)' , target_state=1 , mobile_state=1 , object=alignment_object , gap_max=10 , window=3 )
	print alignment

	cmd.select( name=pdbid_target + '_alignment' , selection='(' + pdbid_target + ') and ' + alignment_object )
	cmd.select( name=pdbid_mobile + '_alignment' , selection='(' + pdbid_mobile + ') and ' + alignment_object )

	stored.target=[]
	stored.mobile=[]

	cmd.iterate( selection=pdbid_target+'_alignment' , expression='stored.target.append([resn,resi,name,chain])' )
	cmd.iterate( selection=pdbid_mobile+'_alignment' , expression='stored.mobile.append([resn,resi,name,chain])' )


	output_file_txt = open('rmsd-by-residue_target-'+pdbid_target+'_mobile-'+pdbid_mobile+'.txt','w')
	output_file_xvg = open('rmsd-by-residue_target-'+pdbid_target+'_mobile-'+pdbid_mobile+'.xvg','w')

	output_file_txt.write("n\ttarget("+pdbid_target+")\tmobile("+pdbid_mobile+")\tRMSD[A]\n")

	for i in range(0,len(stored.target)):
		target_res_name = stored.target[i][0]
		target_res_num  = stored.target[i][1]
		target_atm_name = stored.target[i][2]
		target_chain    = stored.target[i][3]
		
		mobile_res_name = stored.mobile[i][0]
		mobile_res_num  = stored.mobile[i][1]
		mobile_atm_name = stored.mobile[i][2]
		mobile_chain    = stored.mobile[i][3]
		
		distance = cmd.get_distance(
			atom1 = '(name '+target_atm_name+' ) and (resid '+target_res_num+') and ('+pdbid_target+') and (chain '+target_chain+') and ((alt A) or (alt ""))' ,
			atom2 = '(name '+mobile_atm_name+' ) and (resid '+mobile_res_num+') and ('+pdbid_mobile+') and (chain '+mobile_chain+') and ((alt A) or (alt ""))'
		)
		
		print("%s\t%s-%s.%s:%s\t%s-%s.%s:%s\t%.3f\n" % (i+1, target_res_name, target_res_num , target_atm_name , target_chain, mobile_res_name, mobile_res_num , mobile_atm_name , mobile_chain, distance ))
		
		output_file_txt.write("%s\t%s-%s.%s:%s\t%s-%s.%s:%s\t%.3f\n" % (i+1, target_res_name, target_res_num , target_atm_name , target_chain, mobile_res_name, mobile_res_num , mobile_atm_name , mobile_chain, distance ))
		
		output_file_xvg.write("%d\t%.3f\n" % (int(target_res_num) , distance) )


	output_file_txt.write("\n\nOverall RMSD = %.3f A\nAlignment length = %d\n\n" % (alignment['RMSD'] , alignment['alignment_length'] ) )

	output_file_txt.close()
	output_file_xvg.close()
