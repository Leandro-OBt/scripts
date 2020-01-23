#include <stdio.h>
#include <stdlib.h>
#include <math.h>



int main(int argc, const char * argv[]){
	// argv[1] = input files prefix
	// argv[2] = total atoms
	// argv[3] = total frames
	// e.g. for the files traj_ca_ligand_eigenfreq.0.xvg ... traj_ca_ligand_eigenfreq.401.xvg use "./a.out traj_ca_ligand_eigenfreq 0 401"
	
	char filename[100], ch;
	int n;
	
	int frame,total_frames;	// for reading data calculated for all frames
	int total_atoms;
	FILE *input, *output;
	
	double eigenfrequency ;
	
	
	const double h = 6.626e-34 ;	// J.s	Planck constant
	const double kB = 1.380e-23 ;	// J/K	Boltzmann constant
	const double N = 6.02e23 ;	// 		Avogadro constant
	const double c = 3e8 ;		// m/s	Vaccum speed of light
	const double T = 309 ;		// K		Absolute temperature
	const double beta = 1/(kB*T) ;// 1/J	Thermodynamic beta
	double exponential, first_term , second_term, S, S_kcalmolK, TS ;	// will help to calculate the entropy in a more easy-to-read way
	double *S_frame;
	
	
	S = 0.0 ;
	
	sscanf( argv[2], "%d", &total_atoms ) ;
	sscanf( argv[3], "%d", &total_frames ) ;
	
	S_frame = (double *) malloc( (total_frames+1) * sizeof(double) );
	for( n=0 ; n<=total_frames ; n++ ) S_frame[n] = 0.0 ;
	
	
	sprintf( filename , "%s_eigenfreq.dat", argv[1] );
	input = fopen( filename , "r" );
	
	for( frame=1 ; frame<=total_frames ; frame++ ){
		
		for( n=1 ; n<=6 ; n++ ) fscanf( input , "%lf" , &eigenfrequency );	// Reads and ignores the first 6 eigenfrequencies because thei correspond to the translational and rotational degrees of freedom. They could be removed from the input file instead, but I prefer to keep all output from gromacs because it may be useful for future checks.
		for( n=7 ; n<=(total_atoms*3) ; n++ ){	// For all 3n degrees of freedom...
			
			fscanf( input , "%lf" , &eigenfrequency );	// read each eigenfrequency
			
			eigenfrequency *= 100 ;	// converts the eigenfrequency from 1/cm to 1/m
			eigenfrequency *= c ;	// converts the eigenfrequency from 1/m to 1/s
			
			exponential = exp( (-1)*beta*h*eigenfrequency );
	//		printf("%e\n",exponential);
			
			S_frame[frame] += ((((h*eigenfrequency)/T)*(exponential/(1-exponential))) - (kB*log( 1-exponential ))) ;
			
		}
		
		
	}	// The vibrational entropy for all frames were calculated
	fclose(input);
	
	
	
	sprintf( filename, "%s_Svib.dat", argv[1]);
	output = fopen( filename ,"w");
	for( frame=1 ; frame<=total_frames ; frame++ ) fprintf( output, "%e\n",S_frame[frame]);
	fclose( output );
	
	printf("\n\nVibrational entropy (S_vib) calculation.\n\nUsage:\n\t./entropy_vibrational input_file_prefix total_atoms total_frames\n\nThe output file %s was writen. Each line contains the vibrational entropy for the corresponding frame. The output values are in J/K\n\nThe actual input file must be called \"$input_file_prefix\"\"_eigenfreq.dat\" and is a total_frames x total_atoms*3 matrix containing all the eigenfrequencies for a given frame in each line.\n\n" , filename);
	
/*	
	R code to perform bootstrap of the standard error with 200 resamplings
	
list_complex <- read.table("traj_ca_complex_Svib.dat")
data_complex <- as.numeric( unlist(list_complex) )
resamples_complex <- lapply(1:200, function(i) sample(data_complex, replace = T))
averages_complex <- sapply(resamples_complex, mean)

list_receptor <- read.table("traj_ca_receptor_Svib.dat")
data_receptor <- as.numeric( unlist(list_receptor) )
resamples_receptor <- lapply(1:200, function(i) sample(data_receptor, replace = T))
averages_receptor <- sapply(resamples_receptor, mean)

list_ligand <- read.table("traj_ca_ligand_Svib.dat")
data_ligand <- as.numeric( unlist(list_ligand) )
resamples_ligand <- lapply(1:200, function(i) sample(data_ligand, replace = T))
averages_ligand <- sapply(resamples_ligand, mean)

mean(data_complex)
sqrt(var(averages_complex))
mean(data_receptor)
sqrt(var(averages_receptor))
mean(data_ligand)
sqrt(var(averages_ligand))


*/
	
	free(S_frame);
	
	return 0;
}