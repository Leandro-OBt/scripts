#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "leandro.h"






// *********************************************************************
// *                  Functions for writing PDB files                  *
// *********************************************************************




//void output_pdb_conect(FILE *output,int total_atm){
//	int i;
	
//	for( i=1 ; i<=total_atm ; i++ ){
//		fprintf(output,"CONECT%5d%5d\n",i,i+1);
//	}
//}





// *********************************************************************
// *                  Functions for reading PDB files                  *
// *********************************************************************


// Reads the coordinats from a PDB file
//
// >>>>> This relies that the PDB file is "clean", i.e. has only ATOM entries.
void read_initial_configuration(char filename[30] , int total_atm,  atom *p){
	int atm , res;	// atom and residue numbers
	
	res = 0;
	
	FILE *input;
	input=fopen(filename,"rb");
	for( atm=1 ; atm<=total_atm ; atm++ ){
		input_pdb_coordinates( input , atm, p );
		if ( strcmp(p[atm].atm_name,"N")==0 ) res++;
		p[atm].atm_num = atm;
		p[atm].res_num = res;
	}
	
	
	fclose(input);
}






// *********************************************************************
// *                Functions for random number generation             *
// *********************************************************************












void dihedral_rotation(atom *p , int atmB, int atmC, int atmD, float dih){
	
	atom orign, bc, u;
	
	float bc_norm;
	float x,y,z;
	
	// The coordinate frame orign lies in the atom B for the ABCD dihedral torsion
	orign.x = p[atmB].x ;
	orign.y = p[atmB].y ;
	orign.z = p[atmB].z ;
	
	bc.x = p[atmC].x - orign.x;
	bc.y = p[atmC].y - orign.y;
	bc.z = p[atmC].z - orign.z;
				
	bc_norm = sqrt( (bc.x*bc.x)+(bc.y*bc.y)+(bc.z*bc.z) );


	u.x = bc.x / bc_norm ;
	u.y = bc.y / bc_norm ;
	u.z = bc.z / bc_norm ;
	
	
	x = p[atmD].x - orign.x;
	y = p[atmD].y - orign.y;
	z = p[atmD].z - orign.z;
	
	p[atmD].x =  ((u.x*u.x*(1-cos(dih)))+( 1 *cos(dih)))*x + ((u.x*u.y*(1-cos(dih)))-(u.z*sin(dih)))*y + ((u.x*u.z*(1-cos(dih)))+(u.y*sin(dih)))*z;
	p[atmD].y =  ((u.x*u.y*(1-cos(dih)))+(u.z*sin(dih)))*x + ((u.y*u.y*(1-cos(dih)))+( 1 *cos(dih)))*y + ((u.y*u.z*(1-cos(dih)))-(u.x*sin(dih)))*z;
	p[atmD].z =  ((u.x*u.z*(1-cos(dih)))-(u.y*sin(dih)))*x + ((u.z*u.y*(1-cos(dih)))+(u.x*sin(dih)))*y + ((u.z*u.z*(1-cos(dih)))+( 1 *cos(dih)))*z;
	
	p[atmD].x += orign.x ;
	p[atmD].y += orign.y ;
	p[atmD].z += orign.z ;

}



// Calculates the distance between atoms A and B.
//float distance(atom *p, int A, int B){
//	return sqrt( ((p[A].x-p[B].x)*(p[A].x-p[B].x)) + ((p[A].y-p[B].y)*(p[A].y-p[B].y)) + ((p[A].z-p[B].z)*(p[A].z-p[B].z)) );
//}


float energy(atom *p, int total_atm){
	FILE *file;
	int i;
	float E;

	
	file = fopen("new.pdb","wb");
	for( i=1 ; i<=total_atm ; i++ )	output_pdb_coordinates(file,i,p);	// save the structure as new.pdb
	fclose(file);
	
	
	
	system("./energy.sh");	// script for calculating the potential energy of new.pdb and saves it in potential.dat
	
	
	
	file = fopen("potential.dat","rb");
	fscanf(file,"%f",&E);
	fclose(file);
	
	
return E;
}



float rmsd_native(atom *p, int total_atm){
	FILE *file;
	int i;
	float rmsd;
	
	file = fopen("new.pdb","wb");
	for( i=1 ; i<=total_atm ; i++ )	output_pdb_coordinates(file,i,p);	// save the structure as new.pdb
	fclose(file);
	
	
	
	system("./rmsd_native.sh");
	
	
	
	file = fopen("rmsd_ca-ca_native.dat","rb");
	fscanf(file,"%f",&rmsd);
	fclose(file);
	
	
return rmsd;
}



float rmsd_previous(atom *p, int total_atm){
	FILE *file;
	int i;
	float rmsd;
	
	file = fopen("new.pdb","wb");
	for( i=1 ; i<=total_atm ; i++ )	output_pdb_coordinates(file,i,p);	// save the structure as new.pdb
	fclose(file);
	
	system("./rmsd_previous.sh");
	
	file = fopen("rmsd_ca-ca_previous.dat","rb");
	fscanf(file,"%f",&rmsd);
	fclose(file);
	
return rmsd;
}



void build_initial_model(atom *p, int total_atm){
	int i,j;
	float step;
	FILE *sequence;
	char residue;
	
	sequence=fopen("sequence.txt","rb");
	i=0;

	for( i=1 ; i<=total_atm ; i++ ) fscanf(sequence,"%c",&p[i].type);
	
	for( i=0 ; i<=total_atm ; i++ ){
		p[i].x=0.0;
		p[i].y=0.0;
		p[i].z=0.0;
	}
	
	step=3.8;
	
	for( i=1 ; i<=total_atm ; i++ ){
		if( i%2 == 0 ){
			p[i].x = p[i-1].x + step;
			p[i].y = p[i-1].y ;
		}else{
			p[i].x = p[i-1].x ;
			p[i].y = p[i-1].y + step;
		}
	}
	
	
	fclose(sequence);
}


void reject(atom *prev, atom *p, int total_atm, float *e){
	copy_structure( prev , p , total_atm );	// Rejection implies copying the previous structure to the current
	e[1] = e[0] ;
}


void accept(atom *prev, atom *p, int total_atm, float *e){
	copy_structure( p , prev , total_atm );	// Acceptance implies in updating the previous structure with the new one
	e[0] = e[1] ;
}


int number_of_atoms( ){
	int number;
	FILE *sequence;
	char residue;

	sequence=fopen("sequence.txt","rb");
	
	number = 0;
	
	while( fscanf(sequence,"%c",&residue) != EOF ) number++;
	
	number--;	// Correction: "while" loops n+1 times

	fclose(sequence);
	
	return number;
}



int number_of_dihedrals(atom *p , int residue){
	int dih, first;
	
	first=1;
	while( p[first].res_num != residue ) first++;	// "first" is now the first atom of the (residue)th residue
	
	// Number of chi angles...
	if( strcmp(p[first].res_name,"ALA")==0 ) dih=0;
	if( strcmp(p[first].res_name,"CYS")==0 ) dih=1;
	if( strcmp(p[first].res_name,"ASP")==0 ) dih=2;
	if( strcmp(p[first].res_name,"GLU")==0 ) dih=3;
	if( strcmp(p[first].res_name,"PHE")==0 ) dih=2;
	if( strcmp(p[first].res_name,"GLY")==0 ) dih=0;
	if( strcmp(p[first].res_name,"HIS")==0 ) dih=2;
	if( strcmp(p[first].res_name,"ILE")==0 ) dih=2;
	if( strcmp(p[first].res_name,"LYS")==0 ) dih=4;
	if( strcmp(p[first].res_name,"LEU")==0 ) dih=2;
	if( strcmp(p[first].res_name,"MET")==0 ) dih=3;
	if( strcmp(p[first].res_name,"ASN")==0 ) dih=2;
	if( strcmp(p[first].res_name,"PRO")==0 ) dih=0;
	if( strcmp(p[first].res_name,"GLN")==0 ) dih=3;
	if( strcmp(p[first].res_name,"ARG")==0 ) dih=4;
	if( strcmp(p[first].res_name,"SER")==0 ) dih=1;
	if( strcmp(p[first].res_name,"THR")==0 ) dih=1;
	if( strcmp(p[first].res_name,"VAL")==0 ) dih=1;
	if( strcmp(p[first].res_name,"TRP")==0 ) dih=2;
	if( strcmp(p[first].res_name,"TYR")==0 ) dih=2;
	
	// Sums phi, psi and omega
	dih += 3;
	
	return dih;
	
}




void protein_dihedral_mainchain(atom *p, int total_atm, int res, char atmA[5], char atmB[5], char atmC[5], char atmD[5], float angle){
	
	int first, last, i;
	int B, C, D;
	
	first=1;
	while( p[first].res_num != res ) first++;	// "first" is now the first atom of the (res)th residue
	last = first;
	while( p[last].res_num == res ) last++;
	last--;						// "last" is now the last atom of the (res)th residue
//	printf("%d\t%d\n",first,last);	

	// Finds the atom numbers for atmB and atmC, which are fixed.
	for( i=first ; i<=last ; i++){	// scans the residue
		if( strcmp(atmB,p[i].atm_name)==0 ) B=i;
		if( strcmp(atmC,p[i].atm_name)==0 ) C=i;
	}
	// If atmC is "Nn", them this is a omega dihedral
	if( strcmp(atmC,"Nn")==0 ) C=last+1;

	
	// If atmC is CA, then it is a PHI dihedral move
	// In this case, all atoms which are not N nor CA (H, O and side chain) must be moved
	if( strcmp(atmC,"CA")==0 ){
		for( i=first ; i<=last ; i++){
			// If the ith atom is not N, CA nor C...
			if( (strcmp(p[i].atm_name,"N")!=0) && (strcmp(p[i].atm_name,"CA")!=0) && (strcmp(p[i].atm_name,"H")!=0) ){
				// it must be moved
				dihedral_rotation(p,B,C,i,angle);	
			}
		}		
	}
	
	

	// If atmC is C, then it is a PSI dihedral move
	// In this case, only O must be moved
	if( strcmp(atmC,"C")==0 ){
		for( i=first ; i<=last ; i++){
			// If the ith atom is O...
			if( strcmp(p[i].atm_name,"O")==0 ){
				// it must be moved
				dihedral_rotation(p,B,C,i,angle);	
			}
		}		
	}

	
	
	// Finally, moves the atoms of the forward residues

	for( i=(last+1) ; i<=total_atm ; i++ ) dihedral_rotation(p,B,C,i,angle);	
}



void protein_dihedral(atom *p , int first , int last , char atmA[5] , char atmB[5] , char atmC[5] , char atmD[5] , float angle){
	int i;
	int B,C,D;
	
	// finds atmB, atmC and atmD
	for( i=first ; i<=last ; i++){	// scans the whole residue
		if( strcmp(p[i].atm_name,atmB)==0 ) B=i;
		if( strcmp(p[i].atm_name,atmC)==0 ) C=i;
		if( strcmp(p[i].atm_name,atmD)==0 ) D=i;	
	}

	
	dihedral_rotation(p,B,C,D,angle);
}

void protein_dihedral_sidechain(atom *p, int total_atm, int res, int chi, float angle){

	int first, last;
	char a[4],b[4],c[4];
	
	first=1;
	while( p[first].res_num != res ) first++;	// "first" is now the first atom of the (res)th residue
	last = first;
	while( p[last].res_num == res ) last++;
	last--;						// "last" is now the last atom of the (res)th residue
	
	// 0 chi
	// GLY
	// ALA
		
//	if( strcmp(atmA,"Cp")==0 ) ... phi (Cp,N,CA,C)
		
//	if( strcmp(atmC,"Nn")==0 ) ... psi (N,CA,C,Nn)
		
//	if( strcmp(atmD,"CAn")==0) .. omega (CA,C,Nn,CAn)
		
	
//	if(dih==1){
//		protein_dihedral(p,first,last,"Cp","N","CA","C",angle);
//	}
//	if(dih==2) protein_dihedral(p,first,last,"N","CA","C","O",angle);		// ok. psi move apenas O do resíduo atual
//	if(dih==3) protein_dihedral(p,first,last,"CA","C","Nn","CAn",angle);
//	if(dih>3){
	
	// 1 chi, linear
	if( strcmp(p[first].res_name,"SER")==0 ){
		strcpy( a , "N" );
		strcpy( b , "CA" );
		strcpy( c , "CB" );	
		
		protein_dihedral(p,first,last,a,b,c,"OG",angle);
		
		protein_dihedral(p,first,last,a,b,c,"HB1",angle);
		protein_dihedral(p,first,last,a,b,c,"HB2",angle);
		protein_dihedral(p,first,last,a,b,c,"HG",angle);
	}
	
	
	if( strcmp(p[first].res_name,"CYS")==0 ) protein_dihedral(p,first,last,"N","CA","CB","SG",angle);
	
	// 1 chi, branched
	if( strcmp(p[first].res_name,"THR")==0 ){
		protein_dihedral(p,first,last,"N","CA","CB","OG1",angle);
		protein_dihedral(p,first,last,"N","CA","CB","CG2",angle);
	}
	if( strcmp(p[first].res_name,"VAL")==0 ){
		protein_dihedral(p,first,last,"N","CA","CB","CG1",angle);
		protein_dihedral(p,first,last,"N","CA","CB","CG2",angle);
	}
	
	// 2 chi
	if( strcmp(p[first].res_name,"ASP")==0 ){
		if( chi==1 ){
			strcpy( a , "N" );
			strcpy( b , "CA" );
			strcpy( c , "CB" );
			
			protein_dihedral(p,first,last,a,b,c,"CG",angle);
			protein_dihedral(p,first,last,a,b,c,"OD1",angle);
			protein_dihedral(p,first,last,a,b,c,"OD2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HB1",angle);
			protein_dihedral(p,first,last,a,b,c,"HB2",angle);
		}
		if( chi==2 ){
			strcpy( a , "CA" );
			strcpy( b , "CB" );
			strcpy( c , "CG" );
			
			protein_dihedral(p,first,last,a,b,c,"OD1",angle);
			protein_dihedral(p,first,last,a,b,c,"OD2",angle);
		}	
	}
	
	if( strcmp(p[first].res_name,"ASN")==0 ){
		if( chi==1 ){
			strcpy( a , "N" );
			strcpy( b , "CA" );
			strcpy( c , "CB" );
			
			protein_dihedral(p,first,last,a,b,c,"CG",angle);
			protein_dihedral(p,first,last,a,b,c,"OD1",angle);
			protein_dihedral(p,first,last,a,b,c,"ND2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HB1",angle);
			protein_dihedral(p,first,last,a,b,c,"HB2",angle);
			protein_dihedral(p,first,last,a,b,c,"1HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HD2",angle);
			
		}
		if( chi==2 ){
			strcpy( a , "CA" );
			strcpy( b , "CB" );
			strcpy( c , "CG" );
			
			protein_dihedral(p,first,last,a,b,c,"OD1",angle);
			protein_dihedral(p,first,last,a,b,c,"ND2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"1HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HD2",angle);
		}	
	}
	
	if( strcmp(p[first].res_name,"ILE")==0 ){
		if( chi==1 ){
			strcpy( a , "N" );
			strcpy( b , "CA" );
			strcpy( c , "CB" );
			
			protein_dihedral(p,first,last,a,b,c,"CG1",angle);
			protein_dihedral(p,first,last,a,b,c,"CD1",angle);
			protein_dihedral(p,first,last,a,b,c,"CG2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HB",angle);
			protein_dihedral(p,first,last,a,b,c,"1HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"2HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"1HG2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HG2",angle);
			protein_dihedral(p,first,last,a,b,c,"3HG2",angle);
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HD3",angle);
			
		}
		if( chi==2 ){
			strcpy( a , "CA" );
			strcpy( b , "CB" );
			strcpy( c , "CG1" );
			
			protein_dihedral(p,first,last,a,b,c,"CD1",angle);
			
			protein_dihedral(p,first,last,a,b,c,"1HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"2HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HD3",angle);
		}	
	}
	
	if( strcmp(p[first].res_name,"LEU")==0 ){
		if( chi==1 ){
			strcpy( a , "N" );
			strcpy( b , "CA" );
			strcpy( c , "CB" );
			
			protein_dihedral(p,first,last,a,b,c,"CG",angle);
			protein_dihedral(p,first,last,a,b,c,"CD1",angle);
			protein_dihedral(p,first,last,a,b,c,"CD2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HB1",angle);
			protein_dihedral(p,first,last,a,b,c,"HB2",angle);
			protein_dihedral(p,first,last,a,b,c,"HG",angle);
			protein_dihedral(p,first,last,a,b,c,"1HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"2HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"3HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"1HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"3HD2",angle);
		}
		if( chi==2 ){
			strcpy( a , "CA" );
			strcpy( b , "CB" );
			strcpy( c , "CG" );
			
			protein_dihedral(p,first,last,a,b,c,"CD1",angle);
			protein_dihedral(p,first,last,a,b,c,"CD2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HG",angle);
			protein_dihedral(p,first,last,a,b,c,"1HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"2HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"3HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"1HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"3HD2",angle);
			
		}
	}
	
	/*if( strcmp(p[first].res_name,"PRO")==0 ){
		if( chi==1 ){
			protein_dihedral(p,first,last,"N","CA","CB","CG",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CD",angle);
		}
		if( chi==2 ){
			protein_dihedral(p,first,last,"CA","CB","CG","CD",angle);
		}	
	}*/
	
	
	
	if( strcmp(p[first].res_name,"PHE")==0 ){
		if( chi==1 ){
			protein_dihedral(p,first,last,"N","CA","CB","CG",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CD1",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CD2",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CE1",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CE2",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CZ",angle);
		}
		if( chi==2 ){
			protein_dihedral(p,first,last,"CA","CB","CG","CD1",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","CD2",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","CE1",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","CE2",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","CZ",angle);
		}	
	}
	
	if( strcmp(p[first].res_name,"HIS")==0 ){
		if( chi==1 ){
			protein_dihedral(p,first,last,"N","CA","CB","CG",angle);
			protein_dihedral(p,first,last,"N","CA","CB","ND1",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CD2",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CE1",angle);
			protein_dihedral(p,first,last,"N","CA","CB","NE2",angle);
		}
		if( chi==2 ){
			protein_dihedral(p,first,last,"CA","CB","CG","ND1",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","CD2",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","CE1",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","NE2",angle);
		}	
	}
	
	if( strcmp(p[first].res_name,"TYR")==0 ){
		if( chi==1 ){
			strcpy( a , "N" );
			strcpy( b , "CA" );
			strcpy( c , "CB" );			
			
			protein_dihedral(p,first,last,a,b,c,"CG",angle);
			protein_dihedral(p,first,last,a,b,c,"CD1",angle);
			protein_dihedral(p,first,last,a,b,c,"CD2",angle);
			protein_dihedral(p,first,last,a,b,c,"CE1",angle);
			protein_dihedral(p,first,last,a,b,c,"CE2",angle);
			protein_dihedral(p,first,last,a,b,c,"CZ",angle);
			protein_dihedral(p,first,last,a,b,c,"OH",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HB1",angle);
			protein_dihedral(p,first,last,a,b,c,"HB2",angle);
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HE1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE2",angle);
			protein_dihedral(p,first,last,a,b,c,"HH",angle);		
		}
		if( chi==2 ){
			strcpy( a , "CA" );
			strcpy( b , "CB" );
			strcpy( c , "CG" );
			
			protein_dihedral(p,first,last,a,b,c,"CD1",angle);
			protein_dihedral(p,first,last,a,b,c,"CD2",angle);
			protein_dihedral(p,first,last,a,b,c,"CE1",angle);
			protein_dihedral(p,first,last,a,b,c,"CE2",angle);
			protein_dihedral(p,first,last,a,b,c,"CZ",angle);
			protein_dihedral(p,first,last,a,b,c,"OH",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HE1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE2",angle);
			protein_dihedral(p,first,last,a,b,c,"HH",angle);
			
		}	
	}

	if( strcmp(p[first].res_name,"TRP")==0 ){
		if( chi==1 ){
			strcpy( a , "N" );
			strcpy( b , "CA" );
			strcpy( c , "CB" );
			
			protein_dihedral(p,first,last,a,b,c,"CG",angle);
			protein_dihedral(p,first,last,a,b,c,"CD1",angle);
			protein_dihedral(p,first,last,a,b,c,"CD2",angle);
			protein_dihedral(p,first,last,a,b,c,"NE1",angle);
			protein_dihedral(p,first,last,a,b,c,"CE3",angle);
			protein_dihedral(p,first,last,a,b,c,"CE2",angle);
			protein_dihedral(p,first,last,a,b,c,"CZ3",angle);
			protein_dihedral(p,first,last,a,b,c,"CZ2",angle);
			protein_dihedral(p,first,last,a,b,c,"CH2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HB1",angle);
			protein_dihedral(p,first,last,a,b,c,"HB2",angle);
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE3",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ3",angle);
			protein_dihedral(p,first,last,a,b,c,"HH2",angle);
			
		}
		if( chi==2 ){
			strcpy( a , "CA" );
			strcpy( b , "CB" );
			strcpy( c , "CG" );
			
			protein_dihedral(p,first,last,a,b,c,"CD1",angle);
			protein_dihedral(p,first,last,a,b,c,"CD2",angle);
			protein_dihedral(p,first,last,a,b,c,"NE1",angle);
			protein_dihedral(p,first,last,a,b,c,"CE3",angle);
			protein_dihedral(p,first,last,a,b,c,"CE2",angle);
			protein_dihedral(p,first,last,a,b,c,"CZ3",angle);
			protein_dihedral(p,first,last,a,b,c,"CZ2",angle);
			protein_dihedral(p,first,last,a,b,c,"CH2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE3",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ3",angle);
			protein_dihedral(p,first,last,a,b,c,"HH2",angle);
		}	
	}
	
	
	
	// 3 chis
	if( strcmp(p[first].res_name,"MET")==0 ){
		if( chi==1 ){
			protein_dihedral(p,first,last,"N","CA","CB","CG",angle);
			protein_dihedral(p,first,last,"N","CA","CB","SD",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CE",angle);
		}
		if( chi==2 ){
			protein_dihedral(p,first,last,"CA","CB","CG","SD",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","CE",angle);
		}
		if( chi==3 ){
			protein_dihedral(p,first,last,"CB","CG","SD","CE",angle);
		}
	}
	
	if( strcmp(p[first].res_name,"GLN")==0 ){
		if( chi==1 ){
			strcpy( a , "N" );
			strcpy( b , "CA" );
			strcpy( c , "CB" );
			
			protein_dihedral(p,first,last,a,b,c,"CG",angle);
			protein_dihedral(p,first,last,a,b,c,"CD",angle);
			protein_dihedral(p,first,last,a,b,c,"OE1",angle);
			protein_dihedral(p,first,last,a,b,c,"NE2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HB1",angle);
			protein_dihedral(p,first,last,a,b,c,"HB2",angle);
			protein_dihedral(p,first,last,a,b,c,"HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"HG2",angle);
			protein_dihedral(p,first,last,a,b,c,"1HE2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HE2",angle);
			
		}
		if( chi==2 ){
			strcpy( a , "CA" );
			strcpy( b , "CB" );
			strcpy( c , "CG" );
			
			protein_dihedral(p,first,last,a,b,c,"CD",angle);
			protein_dihedral(p,first,last,a,b,c,"OE1",angle);
			protein_dihedral(p,first,last,a,b,c,"NE2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"HG2",angle);
			protein_dihedral(p,first,last,a,b,c,"1HE2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HE2",angle);
		}
		if( chi==3 ){
			strcpy( a , "CB" );
			strcpy( b , "CG" );
			strcpy( c , "CD" );
			
			protein_dihedral(p,first,last,a,b,c,"OE1",angle);
			protein_dihedral(p,first,last,a,b,c,"NE2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"1HE2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HE2",angle);
		}
	}
	
	if( strcmp(p[first].res_name,"GLU")==0 ){
		if( chi==1 ){
			protein_dihedral(p,first,last,"N","CA","CB","CG",angle);
			protein_dihedral(p,first,last,"N","CA","CB","CD",angle);
			protein_dihedral(p,first,last,"N","CA","CB","OE1",angle);
			protein_dihedral(p,first,last,"N","CA","CB","OE2",angle);
		}
		if( chi==2 ){
			protein_dihedral(p,first,last,"CA","CB","CG","CD",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","OE1",angle);
			protein_dihedral(p,first,last,"CA","CB","CG","OE2",angle);
		}
		if( chi==3 ){
			protein_dihedral(p,first,last,"CB","CG","CD","OE1",angle);
			protein_dihedral(p,first,last,"CB","CG","CD","OE2",angle);
		}
	}
	
	
	
	// 4 chi
	if( strcmp(p[first].res_name,"LYS")==0 ){
		if( chi==1 ){
			strcpy( a , "N" );
			strcpy( b , "CA" );
			strcpy( c , "CB" );
			
			protein_dihedral(p,first,last,a,b,c,"CG",angle);
			protein_dihedral(p,first,last,a,b,c,"CD",angle);
			protein_dihedral(p,first,last,a,b,c,"CE",angle);
			protein_dihedral(p,first,last,a,b,c,"NZ",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HB1",angle);
			protein_dihedral(p,first,last,a,b,c,"HB2",angle);
			protein_dihedral(p,first,last,a,b,c,"HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"HG2",angle);
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HE1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ1",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ3",angle);
		}
		if( chi==2 ){
			strcpy( a , "CA" );
			strcpy( b , "CB" );
			strcpy( c , "CG" );
			
			protein_dihedral(p,first,last,a,b,c,"CD",angle);
			protein_dihedral(p,first,last,a,b,c,"CE",angle);
			protein_dihedral(p,first,last,a,b,c,"NZ",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"HG2",angle);
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HE1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ1",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ3",angle);
		}
		if( chi==3 ){
			strcpy( a , "CB" );
			strcpy( b , "CG" );
			strcpy( c , "CD" );
			
			protein_dihedral(p,first,last,a,b,c,"CE",angle);
			protein_dihedral(p,first,last,a,b,c,"NZ",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HE1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ1",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ3",angle);
		}
		if( chi==4 ){
			strcpy( a , "CG" );
			strcpy( b , "CD" );
			strcpy( c , "CE" );
			
			protein_dihedral(p,first,last,a,b,c,"NZ",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HE1",angle);
			protein_dihedral(p,first,last,a,b,c,"HE2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ1",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ2",angle);
			protein_dihedral(p,first,last,a,b,c,"HZ3",angle);
		}
	}
	
	
	// 5 chi
	if( strcmp(p[first].res_name,"ARG")==0 ){
		if( chi==1 ){
			strcpy( a , "N" );
			strcpy( b , "CA" );
			strcpy( c , "CB" );
			
			protein_dihedral(p,first,last,a,b,c,"CG",angle);
			protein_dihedral(p,first,last,a,b,c,"CD",angle);
			protein_dihedral(p,first,last,a,b,c,"NE",angle);
			protein_dihedral(p,first,last,a,b,c,"CZ",angle);
			protein_dihedral(p,first,last,a,b,c,"NH1",angle);
			protein_dihedral(p,first,last,a,b,c,"NH2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HB1",angle);
			protein_dihedral(p,first,last,a,b,c,"HB2",angle);
			protein_dihedral(p,first,last,a,b,c,"HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"HG2",angle);
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HE",angle);
			protein_dihedral(p,first,last,a,b,c,"1HH1",angle);
			protein_dihedral(p,first,last,a,b,c,"2HH1",angle);
			protein_dihedral(p,first,last,a,b,c,"1HH2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HH2",angle);			
		}
		if( chi==2 ){
			strcpy( a , "CA" );
			strcpy( b , "CB" );
			strcpy( c , "CG" );
			
			protein_dihedral(p,first,last,a,b,c,"CD",angle);
			protein_dihedral(p,first,last,a,b,c,"NE",angle);
			protein_dihedral(p,first,last,a,b,c,"CZ",angle);
			protein_dihedral(p,first,last,a,b,c,"NH1",angle);
			protein_dihedral(p,first,last,a,b,c,"NH2",angle);
			
			
			protein_dihedral(p,first,last,a,b,c,"HG1",angle);
			protein_dihedral(p,first,last,a,b,c,"HG2",angle);
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HE",angle);
			protein_dihedral(p,first,last,a,b,c,"1HH1",angle);
			protein_dihedral(p,first,last,a,b,c,"2HH1",angle);
			protein_dihedral(p,first,last,a,b,c,"1HH2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HH2",angle);
		}
		if( chi==3 ){
			strcpy( a , "CB" );
			strcpy( b , "CG" );
			strcpy( c , "CD" );
			
			protein_dihedral(p,first,last,a,b,c,"NE",angle);
			protein_dihedral(p,first,last,a,b,c,"CZ",angle);
			protein_dihedral(p,first,last,a,b,c,"NH1",angle);
			protein_dihedral(p,first,last,a,b,c,"NH2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HD1",angle);
			protein_dihedral(p,first,last,a,b,c,"HD2",angle);
			protein_dihedral(p,first,last,a,b,c,"HE",angle);
			protein_dihedral(p,first,last,a,b,c,"1HH1",angle);
			protein_dihedral(p,first,last,a,b,c,"2HH1",angle);
			protein_dihedral(p,first,last,a,b,c,"1HH2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HH2",angle);
		}
		if( chi==4 ){
			strcpy( a , "CG" );
			strcpy( b , "CD" );
			strcpy( c , "NE" );
			
			protein_dihedral(p,first,last,a,b,c,"CZ",angle);
			protein_dihedral(p,first,last,a,b,c,"NH1",angle);
			protein_dihedral(p,first,last,a,b,c,"NH2",angle);
			
			protein_dihedral(p,first,last,a,b,c,"HE",angle);
			protein_dihedral(p,first,last,a,b,c,"1HH1",angle);
			protein_dihedral(p,first,last,a,b,c,"2HH1",angle);
			protein_dihedral(p,first,last,a,b,c,"1HH2",angle);
			protein_dihedral(p,first,last,a,b,c,"2HH2",angle);
		}		
	}
	
//	}
}




int main(){

// Variable declaration

	int step, atm, i, j;	// loop variables
	int r;			// random number
	float rr;
	int total_atm, total_steps, total_res, total_dihedral;		// total number of atoms
	int res , dihedral ; 		// random residue, dihedral and 
	float angle;
	int rejected;
	float rmsd_caca_native , rmsd_caca_previous;
	
	float prob_phi, prob_psi;
	
	FILE *arq, *file;		// file pointer
	atom *p, *prev;		// data structure
	
	float x,y,z;
	float orign_x, orign_y, orign_z;
	float u_x, u_y, u_z, bc_x, bc_y, bc_z, bc_norm ;
	float e_current, e_previous, *e, e_diff;
	float torsion_amplitude, torsion_direction, pi;
	
	float k_b,T,prob,T_ini,R;
	int freq;
	
	pi=3.1415;
	k_b = 1.3806488E-23 ; //m2kg/s2K
	R = 8.3144621E-3; // kJ/Kmol
	T=309;	// K
	
// Initializations
	srand(time(NULL));	// random number generator
//	total_atm=number_of_atoms();

//	1l2y
	total_atm=304;
	total_res=20;
	
//	1le0
//	total_atm=217;
//	total_res=12;
	
//	total_atm=504;e0
//	total_res=28;
	
	
	
	
	total_steps=1E3;
	freq=1;
	
	torsion_amplitude=30.0;	// degrees
	torsion_amplitude /= (180.00/3.1415);
	
	
// Memmory allocation
	p=(atom *) malloc( (total_atm+1) * sizeof(atom));
	prev=(atom *) malloc( (total_atm+1) * sizeof(atom));
	e=(float *) malloc( 1 * sizeof(float) );
	
// Reads the initial configuration
	read_initial_configuration("1l2y.pdb" , total_atm, p);
	read_initial_configuration("1l2y.pdb" , total_atm, prev);

	//for(i=1;i<=10;i++) printf("%f %f %f\n",p[i].x,p[i].y,p[i].z);
	
	
//	build_initial_model(p,total_atm);
//	build_initial_model(prev,total_atm);
	
	arq=fopen("saida.pdb","wb");
	
	output_pdb_model( arq , 0 );
	for( atm=1 ; atm<=total_atm ; atm++ ) output_pdb_coordinates(arq,atm,p);
	output_pdb_endmdl(arq);
	
	e[0]=energy(p,total_atm);	// Potential energy of the initial structure
	
	printf("# step\taccepted_E\talternative_E\tE_difference\taccept_prob\trmsd_native\trmsd_previous\n");
	printf("# \t[kJ/mol]\t[kJ/mol]\t[kJ/mol]\t\t\t[A]\t\t[A]\n");

	rejected=0;
	// Main loop
	for( step=1 ; step<=total_steps ; step++ ){
		res = random_int(total_res-1)+1;				// Takes a random residue
		total_dihedral = number_of_dihedrals(p,res);		// How many dihedral angle does this residue has?
		dihedral = random_int(total_dihedral-1)+1;		// Takes a random dihedral. 1=phi; 2=psi; 3=omega; 4+ = chi1,2 ...
		
			while( (dihedral==1) && (res==1) )         dihedral = random_int(total_dihedral-1)+1;	// If res==N-term, can't use phi
			while( (dihedral==2) && (res==total_res) ) dihedral = random_int(total_dihedral-1)+1;	// If res==C-term, can't use psi
			while( (dihedral==3) && (res==total_res) ) dihedral = random_int(total_dihedral-1)+1;	// If res==C-term, can't use omega
		
		angle = random_float(torsion_amplitude);	// Takes a random angle [-torsion_amplitude, +torsion_amplitude]
		torsion_direction = random_float(1);		// Takes a random direction for the torsion, i.e. sign of the angle
		if( torsion_direction <= 0.5000 ) angle *= -1;
		
		// Makes the torsions
		if( dihedral == 1 ) protein_dihedral_mainchain(p, total_atm, res,"Cp","N","CA","C",angle);	// phi
		if( dihedral == 2 ) protein_dihedral_mainchain(p, total_atm, res,"N","CA","C","Nn",angle);	// psi
		if( dihedral == 3 ) protein_dihedral_mainchain(p, total_atm, res,"CA","C","Nn","CAn",angle);	// omega		
		if( dihedral > 3 )  protein_dihedral_sidechain(p, total_atm, res, (dihedral-3), angle);		// chis
		
		e[1] = energy(p,total_atm);	// Calculates the energy of the new structure
		e_diff=e[1]-e[0];			// What is the difference between the new energy and the previous one?
		
		printf("%d\t%15e\t%15e",step,e[0],e[1]);
		
		// Checks the new structure acceptance
		file = fopen("old.pdb","wb");
		for( i=1 ; i<=total_atm ; i++ ) output_pdb_coordinates(file,i,prev);	// save the previous structure as old.pdb
		fclose(file);
		
		if( e[1]>e[0] ){		// If the energy of the new structure is higher than of the previous...
			prob=exp( (-1)*((e[1]-e[0])/(R*T) ) );	// Metropolis criterion
			rr=random_float(1) + 0.000001;	// correction for not allowing 0.00000000 to be chosen
			
			if( rr <= prob ){
				accept(prev, p, total_atm , e);
				rejected=0;
			}else{
				reject(prev, p, total_atm , e);
				rejected++;
				e_diff=0.0;
			}
		}else{
			prob=1.000;
			accept(prev, p, total_atm , e);
			rejected=0;
		}
		
		printf("\t%10.5f\t%8.5f",e_diff,prob);
		
		
		
		// PDB output
		if( step % freq == 0 ){
			output_pdb_model( arq , step );
			for( atm=1 ; atm<=total_atm ; atm++ ){	
				output_pdb_coordinates(arq,atm,p);
			}
			output_pdb_endmdl(arq);
		}
		
		rmsd_caca_native = rmsd_native( p , total_atm );
		rmsd_caca_previous = rmsd_previous( p, total_atm );
		
		printf("\t%.2f\t%.2f\t%d\n", rmsd_caca_native , rmsd_caca_previous, rejected);
		
	}
	
//	output_pdb_conect(arq,total_atm);
	fclose(arq);
	free(e);
	free(prev);
	free(p);
	
return 0;
}