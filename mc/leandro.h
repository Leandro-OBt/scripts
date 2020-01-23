#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
	char res_name[5];
	int res_num;	
	float x, y, z;	// coordinates
	char type;	// atom type
	int atm_num;
	char element[2];
	char atm_name[5];
} atom ;









// Writes the model number. This must be the first line of each model in the file
void output_pdb_model( FILE *output,int n_model ){
	fprintf(output,"MODEL %d\n",n_model);
}

// Writes the full coordinate line for a given atom "atm" in the file.
void output_pdb_coordinates( FILE *output , int atm , atom *p ){
	fprintf(output,"ATOM  %5d %4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00           %s\n", atm , p[atm].atm_name , p[atm].res_name, p[atm].res_num , p[atm].x , p[atm].y , p[atm].z , p[atm].element);	
}

// Writes "ENDMDL". This must be the last line of each model in the file
void output_pdb_endmdl( FILE *output ){
	fprintf(output,"ENDMDL\n");
}


void output_pdb( FILE *output , atom *p , int total_atm, int model ){
	int i;
	
	output_pdb_model( output , i );
	for( i=1 ; i<=total_atm ; i++ ) output_pdb_coordinates(output, i, p);
	output_pdb_endmdl(output);

}
	







// Reads one coordinate line of the file.
// !! This function relies in the full PDB file format, with chain and element information !!
void input_pdb_coordinates(FILE *input, int atm , atom *p){
	char garbage[20];
	fscanf(input,"%s%s%s%s%s%s%f%f%f%s%s%s",&garbage[0],&garbage[0],&p[atm].atm_name[0],&p[atm].res_name[0],&garbage[0],&garbage[0],&p[atm].x,&p[atm].y,&p[atm].z,&garbage[0],&garbage[0],&p[atm].element[0]);
	
	//fscanf(input,"%s%s%s%s%s%f%f%f%s%s%s",&garbage[0],&garbage[0],&p[atm].atm_name[0],&p[atm].res_name[0],&garbage[0],&p[atm].x,&p[atm].y,&p[atm].z,&garbage[0],&garbage[0],&p[atm].element[0]);
	//fscanf(input,"%s%s%s%s%s%f%f%f%s%s",&garbage[0],&garbage[0],&p[atm].atm_name[0],&p[atm].res_name[0],&garbage[0],&p[atm].x,&p[atm].y,&p[atm].z,&garbage[0],&garbage[0]);
	
//	fscanf(input,"%s%s%f%f%f" , &p[atm].atm_name[0] , &p[atm].res_name[0] , &p[atm].x , &p[atm].y , &p[atm].z );
	
	
//	ATOM     22  HD1 HIS     2      39.600  61.520  29.570  1.00  0.00            
}



/*
void input_pdb_coordinates(FILE *input, int atm , atom *p){
	char garbage[20];	//fscanf(input,"%s%s%s%s%s%s%f%f%f%s%s%s",&garbage[0],&garbage[0],&p[atm].atm_name[0],&p[atm].res_name[0],&garbage[0],&garbage[0],&p[atm].x,&p[atm].y,&p[atm].z,&garbage[0],&garbage[0],&p[atm].element[0]);
	
	//fscanf(input,"%s%s%s%s%s%f%f%f%s%s%s",&garbage[0],&garbage[0],&p[atm].atm_name[0],&p[atm].res_name[0],&garbage[0],&p[atm].x,&p[atm].y,&p[atm].z,&garbage[0],&garbage[0],&p[atm].element[0]);
	fscanf(input,"%s%s%s%s%s%f%f%f%s%s",&garbage[0],&garbage[0],&p[atm].atm_name[0],&p[atm].res_name[0],&garbage[0],&p[atm].x,&p[atm].y,&p[atm].z,&garbage[0],&garbage[0]);
	
	
	
//	ATOM     22  HD1 HIS     2      39.600  61.520  29.570  1.00  0.00            
}
*/


void read_pdb(char filename[200] , int total_atm,  atom *p){
	int atm , res;	// atom and residue numbers
	char command[500];
	
	res = 0;
	
	sprintf( command, "grep ATOM %s > temporary_pdb_atom ; cut -c 12-21 temporary_pdb_atom > temporary_atm ; cut -c 30-55 temporary_pdb_atom > temporary_coord ; paste temporary_atm temporary_coord > temporary_read_PDB_coordinates", filename);
	system(command);
	
	FILE *input;
	input=fopen("temporary_read_PDB_coordinates","rb");
	for( atm=1 ; atm<=total_atm ; atm++ ){
		input_pdb_coordinates( input , atm, p );
		if ( strcmp(p[atm].atm_name,"N")==0 ) res++;
		p[atm].atm_num = atm;
		p[atm].res_num = res;
	}
	
	
	fclose(input);
	system("rm temporary_pdb_atom temporary_atm temporary_coord temporary_read_PDB_coordinates");
}








float distance(atom *p, int A, int B){
	return sqrt( ((p[A].x-p[B].x)*(p[A].x-p[B].x)) + ((p[A].y-p[B].y)*(p[A].y-p[B].y)) + ((p[A].z-p[B].z)*(p[A].z-p[B].z)) );
}




float distance_inter(atom *a, int A, atom *b, int B){
	return sqrt( pow((a[A].x-b[B].x),2) + pow((a[A].y-b[B].y),2) + pow((a[A].z-b[B].z),2) );
}





// Returns an integer random number between 0 and "max"
int random_int(int max){
	int r;
	r = rand() % (max+1);
return r;
}

// Returns a float random number between 0 and "max". The precision can be tuned by the constant "C"
float random_float(float max){
	float r;
	int C;
	int temp;
	
	C=10000;
	
	temp = (int) (max * C);
	
	r = rand() % (temp+1);
	r /= C;
return r;	
}

void random_seed(){
	srand(time(NULL));
}










void copy_structure(atom *from , atom *to , int total_atm){
	int i;
	
	for( i=1 ; i<=total_atm ; i++ ){
		to[i].x = from[i].x ;
		to[i].y = from[i].y ;
		to[i].z = from[i].z ;
	}
}
















// This function computes the contact map from the coordinates of "atom *p". The contact matrix is stored in "int **matrix"
// Only the CA atoms are considered
//
// cutoff is the distance criterion for above which a pair of residues are considered as contacting
// alpha is the number of neighbor residues which must be ignored.
//
// March 04, 2014
int contact_map_ca( int total_atm, atom *p, int alpha, float cutoff, int **matrix){
	int i,j;
	int total_contacts;
	
	total_contacts=0;
	
	for( i=1 ; i<=total_atm ; i++){
		for( j=(i+1); j<=total_atm ; j++ ){
			if( strcmp(p[i].atm_name,"CA")==0 && strcmp(p[j].atm_name,"CA")==0 ){
				if( p[j].res_num > (p[i].res_num + alpha) ){
					if( distance(p, i ,j) <= cutoff ){
						total_contacts++;
						matrix[ p[i].res_num ][ p[j].res_num ]++;
					}
				}
			}
		}
	}	
	
	return total_contacts;
}



// This function computes the contact map from the coordinates of "atom *p". The contact matrix is stored in "int **matrix"
//
// cutoff is the distance criterion for above which a pair of residues are considered as contacting
// alpha is the number of neighbor residues which must be ignored.
//
// March 04, 2014
int contact_map( int total_atm, atom *p, int alpha, float cutoff, int **matrix){
	int i,j;
	int total_contacts;
	
	total_contacts=0;
	
	for( i=1 ; i<=total_atm ; i++){
		for( j=(i+1); j<=total_atm ; j++ ){
			if( p[j].res_num > (p[i].res_num + alpha) ){
				if( distance(p, i ,j) <= cutoff ){
					total_contacts++;
					matrix[ p[i].res_num ][ p[j].res_num ]++;
				}
			}
		}
	}	
	
	return total_contacts;
}






// Function for dynamic memory allocation of an array of the "atom" struct
//
// March 04, 2014
atom *allocate_atom( atom *p , int total_atm ){
	p = (atom *) malloc( (total_atm+1) * sizeof(atom) );
return p;
}

// Function for memory liberation of a previously allocated array of the "atom" struct
//
// March 04, 2014
void free_atom( atom *p , int total_atm){
	free(p);
}


// Function for dynamic memory allocation of a 1D vector of integers
//
// May 25, 2014
int *allocate_int( int *a , int size){

	a = (int *) malloc( (size+1) * sizeof(int) );
	
	return a;
}

// Function for memory liberation of a previously allocated 1D vector of integers
//
// May 25, 2014
void free_int( int *a , int size ){
	free( a );
}

// Function for initializing a 1D vector of integers
//
// May 25, 2014
void initialize_int( int *a , int size , int value){
	int i;
	
	for( i=0 ; i<=size ; i++ ){
		a[i] = value ;
	}
	
}


// Function for dynamic memory allocation of a 1D vector of floats
//
// May 25, 2014
float *allocate_float( float *a , int size){

	a = (float *) malloc( (size+1) * sizeof(float) );
	
	return a;
}

// Function for memory liberation of a previously allocated 1D vector of integers
//
// May 25, 2014
void free_float( float *a , int size ){
	free( a );
}

// Function for initializing a 1D vector of integers
//
// May 25, 2014
void initialize_float( float *a , int size , float value){
	int i;
	
	for( i=0 ; i<=size ; i++ ){
		a[i] = value ;
	}
	
}


// Function for dynamic memory allocation of a 2D matrix of integers
//
// March 04, 2014
int **allocate_int_2d( int **a , int rows, int columns ){
	int i;
	
	a = (int **) malloc( (rows+1) * sizeof(int *) );
	for( i=0 ; i<=rows ; i++ ) a[i] = (int *) malloc( (columns+1) * sizeof(int) );
	
	return a;
}

// Function for memory liberation of a previously allocated 2D matrix of integers
//
// March 04, 2014
void free_int_2d( int **a , int rows, int columns ){
	int i;
	
	for( i=0 ; i<=rows ; i++ ) free( a[i] );
	free( a );
}


// Function for initializing a 2D matrix of integers
//
// March 04, 2014
void initialize_int_2d( int **a , int rows, int columns , int value){
	int i, j;
	
	for( i=0 ; i<=rows ; i++ ){
		for( j=0 ; j<=columns ; j++ ){
			a[i][j] = value ;
		}
	}
	
}



// Function for dynamic memory allocation of a 3D matrix of integers
//
// May 02, 2014
int ***allocate_int_3d( int ***a , int x, int y, int z){
	int i,j;
	
	a = (int ***) malloc( (x+1) * sizeof(int) );
	for( i=0 ; i<=x ; i++ ) a[i] = (int **) malloc( (y+1) * sizeof(int) );
	for( i=0 ; i<=x ; i++ ){
		for( j=0 ; j<=y ; j++ ) a[i][j] = (int *) malloc( (z+1) * sizeof(int) );
	}
	
	return a;
}

// Function for memory liberation of a previously allocated 2D matrix of integers
//
// May 02, 2014
void free_int_3d( int ***a , int x, int y, int z ){
	int i,j;
	
	for( i=0 ; i<=x ; i++ ){
		for( j=0 ; j<=y ; j++ ) free( a[i][j] );
	}	
	for( i=0 ; i<=x ; i++ ) free( a[i] );
	free( a );
}


// Function for initializing a 3D matrix of integers
//
// May 02, 2014
void initialize_int_3d( int ***a , int x, int y , int z, int value){
	int i,j,k;
	
	for( i=0 ; i<=x ; i++ ){
		for( j=0 ; j<=y ; j++ ){
			for( k=0 ; k<=z ; k++ ){
				a[i][j][k] = value ;
				printf("%d\t%d\t%d\n",i,j,k);
			}
		}
	}
	
}
