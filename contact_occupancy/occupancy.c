#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char *argv[]){
	
	typedef struct {
		int res_num;	
		float x, y, z;	// coordinates
	} atom ;
	
	int i;
	int frame, total_frames, rep;
	char input_filename[50], output_filename[50];
	int atm, protein_atm, protein_atm_total, ligand_atm, ligand_atm_total, protein_res_total, res;
	float distance, distance_cutoff;
	int *contacts_atm, *contacts_res, *contacts_res_frame, total_contacts_atm;
	float *occupancy, *density, *relative_frequency;
	FILE *input, *output, *parameters;
	atom *protein, *ligand;
	
	sscanf( argv[1] , "%f", &distance_cutoff );
	sscanf( argv[2] , "%d", &rep );
	
	parameters = fopen( "parameters.dat" , "r" );
	
	fscanf( parameters, "%d", &total_frames );
	fscanf( parameters, "%d", &protein_atm_total );
	fscanf( parameters, "%d", &ligand_atm_total );
	fscanf( parameters, "%d", &protein_res_total );
	
	fclose(parameters);
	
//	total_frames=30003;
//	protein_atm_total=2038;
//	ligand_atm_total=12;
//	protein_res_total=132;
//	distance_cutoff = 3.0 ; // Angstrom
	
	
	
	protein = (atom *) malloc( (protein_atm_total+1) * sizeof(atom) );
	ligand =  (atom *) malloc( (ligand_atm_total+1) * sizeof(atom) );
	
	contacts_atm = malloc( (protein_res_total+1) * sizeof(int) );
	contacts_res = malloc( (protein_res_total+1) * sizeof(int) );
	contacts_res_frame = malloc( (protein_res_total+1) * sizeof(int) );
	
	occupancy = malloc( (protein_res_total+1) * sizeof(float));	// contact occupancy
	density = malloc( (protein_res_total+1) * sizeof(float));	// contact density
	relative_frequency = malloc( (protein_res_total+1) * sizeof(float));	// relative contact frequency
	
	
	for( i=0; i<=protein_res_total ; i++ ){ 
		contacts_atm[i] = 0 ;
		contacts_res[i] = 0 ;
		contacts_res_frame[i] = 0 ;
		occupancy[i]= 0.0 ;
		density[i] = 0.0 ;
		relative_frequency[i] = 0.0 ;
	}
	total_contacts_atm=0;
	
	

	for( frame = 0 ; frame < total_frames ; frame++ ){	// for each frame ...
		
		
		sprintf( input_filename , "frames/%d/frame.%d.pdb", rep, frame );	// Build the input filename
		
		input = fopen( input_filename , "r" );	// Opens the input file for reading
		
		
		// Read the input file coordinates for protein atoms and then for ligand atoms
		for( atm=1 ; atm<=protein_atm_total ; atm++ ) fscanf( input , "%d %f %f %f", &protein[atm].res_num, &protein[atm].x, &protein[atm].y, &protein[atm].z );
		for( atm=1 ; atm<=ligand_atm_total ; atm++ ) fscanf( input , "%d %f %f %f", &ligand[atm].res_num, &ligand[atm].x, &ligand[atm].y, &ligand[atm].z );
		
		
		fclose( input );
			

		for( ligand_atm = 1 ; ligand_atm <= ligand_atm_total ; ligand_atm++ ){	// for each ligand atom...
			for( protein_atm = 1 ; protein_atm <= protein_atm_total ; protein_atm++ ){	// for each protein atom ...
				
				// distance between each pair of protein/ligand atoms
				distance = sqrt(	pow( ligand[ligand_atm].x - protein[protein_atm].x , 2) + 
							pow( ligand[ligand_atm].y - protein[protein_atm].y , 2) +
							pow( ligand[ligand_atm].z - protein[protein_atm].z , 2)
						);
				
				if ( distance <= distance_cutoff ){
//					printf("%f\t",distance);

					contacts_atm[ protein[protein_atm].res_num ]++ ;	// for contact density calculation
					contacts_res_frame[ protein[protein_atm].res_num ] = 1 ;	// for occupancy calculation
					
					total_contacts_atm++;					
				}
				
				
			}
		}
		for( res=0; res<=protein_res_total ; res++ ){
			contacts_res[res] += contacts_res_frame[res] ;
			contacts_res_frame[res] = 0;
		}
	}
	
	
	for( res=0; res<=protein_res_total ; res++ ){
		occupancy[res] = contacts_res[res] / (float) total_frames ;
		density[res]   = contacts_atm[res] / (float) total_frames ;
		relative_frequency[res] = contacts_atm[res] / (float) total_contacts_atm ;
	}
	
	// occupancy = In which percentage of the frames there is contact between the protein and the ligand. For each residue, min=0.0 (no contact in the whole simulation) max=1.0 (contact in every frame of the simulation).
	// density = How many contacts each residue does with the ligand in average
	// relative_frequency = Contribution of each residue to the total number of contacts between the protein and the ligand. Must sum 1 for all residues.
	
	sprintf( output_filename , "contact_info_%.1fA.%d.xvg", distance_cutoff, rep );
	output = fopen( output_filename , "w");
	for( res=0; res<=protein_res_total ; res++ ){
		//fprintf( output , "%d\t%f\t%f\t%f\n",res,occupancy[res],density[res],relative_frequency[res]);
		fprintf( output , "%d\t%f\n",res,occupancy[res]);
	}	
	fclose(output);
	
/*	sprintf( output_filename , "bf_occupancy_%.1fA.dat", distance_cutoff );
	output = fopen( output_filename , "w");
	fprintf( output , "%d\n", protein_res_total );
	for( res=1; res<=protein_res_total ; res++ ){
		fprintf( output , "%d\t%f\n",res,occupancy[res]);
	}
	fclose(output);
	
	sprintf( output_filename , "bf_density_%.1fA.dat", distance_cutoff );
	output = fopen( output_filename , "w");
	fprintf( output , "%d\n", protein_res_total );
	for( res=1; res<=protein_res_total ; res++ ){
		fprintf( output , "%d\t%f\n",res,density[res]);
	}
	fclose(output);
	
	sprintf( output_filename , "bf_relative_frequency_%.1fA.dat", distance_cutoff );
	output = fopen( output_filename , "w");
	fprintf( output , "%d\n", protein_res_total );
	for( res=1; res<=protein_res_total ; res++ ){
		fprintf( output , "%d\t%f\n",res,relative_frequency[res]);
	}
	fclose(output);
*/


return 0;
}
