#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, const char * argv[]){
	
	int total_atm,atm;
	int total_frames,frame;
	float unit_volume, bulk_occ;
	char input_filename_prefix[100], output_filename_prefix[100], input_filename[100], output_filename[100];
	float x_min,y_min,z_min,x_max,y_max,z_max,step,x_avg,y_avg,z_avg;	// Angstrom
	float protein_com_x, protein_com_y, protein_com_z, box_size ;
	int x_num,y_num,z_num;
	float x,y,z;
	int a,b,c;
	float ***occ;
	FILE *input, *output;
	float occ_avg, occ_stdev;
	int n_box_elements;
	int number_of_stdev , n;
	float R, T, K, dG;
	
	R = 0.001987 ; // kcal/mol.K
	T = 309  ;
	
	sscanf( argv[1], "%s", &input_filename_prefix[0]);
	sscanf( argv[2], "%d", &total_atm);
	sscanf( argv[3], "%d", &total_frames);
	sscanf( argv[4], "%s", &output_filename_prefix[0]);
//	total_frames=2501;
//	total_atm=17023;
	
	protein_com_x=60.743; // A
	protein_com_y=60.986; // A
	protein_com_z=28.915; // A
	box_size=40.0; // A
	
	
	x_max = protein_com_x + box_size ;
	y_max = protein_com_y + box_size ;
	z_max = protein_com_z + box_size ;
	x_min = protein_com_x - box_size ;
	y_min = protein_com_y - box_size ;
	z_min = protein_com_z - box_size ;
	
	
	
	step  = 1.0 ; // A
	unit_volume = step * step * step ;	// cubic A. The volume of each box element
	bulk_occ = 0.0333; // The density of 0.0333... OW per cubic Angstrom corresponds to the usual 1 kg/L desnityity of water
	
	x_num = ((x_max - x_min)/step) + 1;
	y_num = ((y_max - y_min)/step) + 1;
	z_num = ((z_max - z_min)/step) + 1;
	
	
//	printf("\n\n%d\t%d\t%d\n\n", x_num, y_num, z_num );
	
	// occ[x][y][z]
	occ = (float ***) malloc( (x_num+1) * sizeof(float **) );
	for( a=0 ; a<=x_num ; a++ ){
		occ[a] = (float **) malloc( (y_num+1) * sizeof(float *) );
		for( b=0 ; b<=y_num ; b++ ){
			occ[a][b] = (float *) malloc( (z_num+1) * sizeof(float) );
		}
	}

	for( a=0 ; a<x_num ; a++ ){
		for( b=0 ; b<y_num ; b++ ){
			for( c=0 ; c<z_num ; c++ ){
				occ[a][b][c]=0.0;
//				printf("%d\t%d\t%d\t%f\n", a, b, c, occ[a][b][c]);
			}
		}
	}
		
	
	
	for( frame=1 ; frame<=total_frames ; frame++ ){
		
		sprintf( input_filename , "%s_%d.dat", input_filename_prefix, frame );
		printf("%s ... ",input_filename);
		
		input = fopen( input_filename , "r" );
		printf("ok\n");
	
		for( atm=1 ; atm<=total_atm ; atm++ ){

			fscanf( input , "%f %f %f", &x, &y, &z) ;
			
			if( (x>=x_min) && (x<=x_max) && (y>=y_min) && (y<=y_max) && (z>=z_min) && (z<=z_max) ){
			
				x_num = ((x - x_min)/step) + 1;
				y_num = ((y - y_min)/step) + 1;
				z_num = ((z - z_min)/step) + 1;
			
				//printf("%d\t%d\t%d\n", x_num, y_num, z_num );
				occ[ x_num ][ y_num ][ z_num ] += 1.0 ;
			}
		}
		
		fclose(input);
	}
	
	
	
	
	
	
	occ_avg=0.0 ;
	occ_stdev=0.0 ;
	n_box_elements=0.0 ;
	
	x_num = ((x_max - x_min)/step) + 1 ;
	y_num = ((y_max - y_min)/step) + 1 ;
	z_num = ((z_max - z_min)/step) + 1 ;
	
	for( a=0 ; a<x_num ; a++ ){
		for( b=0 ; b<y_num ; b++ ){
			for( c=0 ; c<z_num ; c++ ){
				
				occ[a][b][c] /= total_frames ; // Occupancy of OW atoms in each box element
				
				occ_avg += occ[a][b][c] ;
				n_box_elements++ ;
				
			}
		}
	}
	
	occ_avg /= n_box_elements ;
	
//	for( a=0 ; a<x_num ; a++ ){
//		for( b=0 ; b<y_num ; b++ ){
//			for( c=0 ; c<z_num ; c++ ){
//				if ( occ[a][b][c] > 0.000 ){
//					occ_stdev += pow( (occ[a][b][c] - occ_avg), 2 );
//				}
//			}
//		}
//	}
//	occ_stdev /= (n_box_elements-1) ;
//	occ_stdev = sqrt( occ_stdev );
//	printf("avg = %f\nstdev = %f\nn_box= %d\n\n", occ_avg, occ_stdev, n_box_elements );
	printf("avg = %f\nn_box= %d\n\n", occ_avg, n_box_elements );
	
	
	
	sprintf( output_filename, "%s.pdb", output_filename_prefix);
	output = fopen( output_filename ,"w");
	
	n=1 ;
	for( x=x_min ; x<=x_max ; x+=step ){
		for( y=y_min ; y<=y_max ; y+=step ){
			for( z=z_min ; z<=z_max ; z+=step ){
				
				x_num = ((x - x_min)/step) + 1;
				y_num = ((y - y_min)/step) + 1;
				z_num = ((z - z_min)/step) + 1;
				
				x_avg = x + (step/2);
				y_avg = y + (step/2);
				z_avg = z + (step/2);
				
				
				if( occ[ x_num ][ y_num ][ z_num ] > 0.000 ){
					
					K = occ[ x_num ][ y_num ][ z_num ] / occ_avg;
					
					dG = (-1)*R*T*log(K) ;
					
					if( abs(dG) >= 0.6 ){
						fprintf( output , "ATOM %6d  OW  SOL%6d    %8.3f%8.3f%8.3f  1.00%6.3f\n", n, n , x_avg, y_avg, z_avg , dG ) ;
					}
					
					//n++;
				}
				
			}
		}
	}
	
	fclose(output);

	
	return 0;
}
