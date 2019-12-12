#include <stdio.h>
#include <stdlib.h>
#include <gtk/gtk.h>
#include <math.h>
#include <sys/time.h>
#include "structSet.h"
#include <mcheck.h>

#define PI 3.141592
#define gamma 6.67384e-11 
#define si2mg 1e5

int InfiNet			  = 0; // krotnosc rozciągnięcia siatki *InfiNet*DimX
int areFieldMGData     = 0; // Flag: are Gravi field data read?
int areFwdSolvers     = 0; // Flag: are the fwd solvers done?
int isShiftCalculated = 0;
int isRegular         = 1; // Flag: is mesh regular? 0-no 1-yes
int MGIsIn 			  = 1;

// ******* Functions definition
void Magneto_new(double *p_model_params, double *misfit);
void alloc_2D_array(double ***p_array, int p_dim_x, int p_dim_y);
double anorm(int n, double *p_d);
void calculate_mag(block *p_net, model* calc_model, double **MG_results_temp, int p_MG_sites, double **p_MG_data);
int count_sites(char *p_file);
void get_input_files(char *model_file, char *g_file, char *mg_file);
double Magneto(model *calc_model, block *p_net, int p_MG_sites, double **p_MG_data, double **p_MG_results);
void open_MG_data_auto(char *g_file, double **p_MG_data);
double tab_min(double **tab, int nelems);
double tab_max(double **tab, int nelems);
void read_line(char *p_line, FILE *p_plik);
void read_line_id(char *p_line);
block *read_model_auto(model *p_model1, char *model_file, block *p_net, GSList *p_ListOf_net_z_range);
void run_fwd_solvers(model *p_model1, block *p_net, int p_MG_sites, double *p_MGmisfit, double **p_MG_data, double **p_MG_results);
void run_rasterization(model *p_model1, block *p_net, GSList *p_ListOf_net_z_range);
void thicken_mesh(model *p_model1, GSList *p_entThMeshParams, block *p_net, GSList *p_ListOf_net_z_range); 

// ******* From ModelRaster
int ccw(vertex a, vertex b, vertex c);
int check_edges_intersection(vertex *vrx, int p_NumVrx, int l);
long det_matrix(vertex x, vertex y, vertex z);
int fill_mesh_with_model(model *p_model1, block *p_net);
int fun_max(double a, double b);
int fun_min(double a, double b);
gpointer g_list_id_get_data(model *modeltemp, intptr_t id, int flag);
int get_points_count(model *modeltemp);
int if_belongs(vertex x, vertex y, vertex z);
int layer_vertices_counter (model *modeltemp, int LayerId);
vertex *layers_to_polygons(model *modeltemp, vertex *vrx, int l, int *p_NumVrx, double *p_Rhob, double *p_Res, double *p_Inc, double *p_Ints, double *p_AngM, double *p_Suscept/*, colour *fillColour*/);
block *net_maker(model *p_model1, int isRegular, block *p_net, GSList *p_ListOf_net_z_range);
int sign(long a);
// *******

void alloc_2D_array(double ***p_array, int p_dim_x, int p_dim_y)
{
	int i = 0;

	(*p_array) = (double**)malloc(p_dim_x*sizeof(double*)); // alokacja pamieci na dane polowe
	for(i=0; i<p_dim_x; i++)
		(*p_array)[i]=(double*)malloc(p_dim_y*sizeof(double));
					
	if(p_array == NULL) 
	{ 
		g_print("\n Przydzielenie pamięci dla tablicy nie było możliwe"); 
		getchar(); 
	}	

	return;
}

// function checks if points a, b, c form a counterclockwise angle 
// when triangle made from a, b, c has positive area, a->b->c are counterclockwise - return 1
// when triangle made from a, b, c has negative area, a->b->c are clockwise - return -1
// when triangle made from a, b, c has zero area, a->b->c are collinear - return 0
int ccw(vertex a, vertex b, vertex c)
{
	return (b.x - a.x) * (c.z - a.z) - (c.x - a.x) * (b.z - a.z);
}

// Funkcja sprawdzajaca czy punkt lezy na krawedzi i czy krawedzie sie przecinaja
// 1 - przecinaja sie lub wierzcholek na krawedzi
// 0 - nie przecinaja sie i nie ma wierzcholkow na krawedziach
int check_edges_intersection(vertex *vrx, int p_NumVrx, int l)
{
	int i, k = 0;
	 
	for(i=0; i<p_NumVrx; i++)  // Sprawdzanie, czy punkt vrx[k] lezy na odcinku |vrx[i]vrx[i+1]|
	{		
		if( (i+1) == p_NumVrx ) // ostatnia krawedz
			k = 1;
		else
			k = 0;
					
		for(; k<p_NumVrx; k++)
		{							 	
			if( if_belongs(vrx[i],vrx[i+1],vrx[k]) == 0 && k!=i && k!=i+1 )																	
				return 1;	
		}
	} 
	for(i=0; i<p_NumVrx; i++)  // sprawdzenie czy krawedzie sie przecinaja
	{		
		for(k=i+2; k<=p_NumVrx-1; k++)
		//for(k=i+2; k<p_NumVrx-1; k++)
		{  													   
			if( (sign (det_matrix(vrx[i],vrx[i+1],vrx[k])) != sign (det_matrix(vrx[i],vrx[i+1],vrx[k+1]))) &&
				( (sign (det_matrix(vrx[k],vrx[k+1],vrx[i])) != sign (det_matrix(vrx[k],vrx[k+1],vrx[i+1])))) )
				return 1; 
			
		}
	}   
	return 0;
}

long det_matrix(vertex a, vertex b, vertex c) 
{	
	return(a.x*b.z + b.x*c.z + c.x*a.z - c.x*b.z - a.x*c.z - b.x*a.z);		
}

int fill_mesh_with_model(model *p_model1, block *p_net)
{	
	struct rasterization // zrzutowane krawedzie na siatke
	{
		int index; 
		int flag;
	};
	
	vertex *vrx = NULL;
	struct rasterization *PolySeg; 	
	struct rasterization tempPolySeg;
		
	int *CommPolyZ;
	int NumPoly;        // liczba wielokatow  
	int NumVrx;         // ilosc wierzcholkow poszczegolnego wielokata
	double x1,z1,x2,z2; // wsp. poszczegolnego wierzcholka
	double Rhob = 0.0;    // gestosc danego wielokata    
	double Res = 0.0;     // opornosc danego wielokata
	double Inc = 0.0;
	double Ints = 0.0;
	double AngM = 0.0;
	double Suscept = 0.0;
	
	int flag2 = 0;
	int k=0; 
	int i=0;
	int j=0;
	int n=0;
	int l=0;
	int i_net=0; 
	int NumNetX = p_model1->NumNetX;  // ilość oczek liczac po X  
	int NumNetZ = p_model1->NumNetZ;  // ilość oczek liczac po Z
	int NumNet = NumNetX * NumNetZ;    // ilość oczek siatki

	double a=0.0;
	double b=0.0;

	NumPoly = g_slist_length(p_model1->listOfLayersId);  							
	for(l=0; l<NumPoly; l++) // rasteryzacja kazdego wielokata  (warstwy) 
	{ 			
		vrx = layers_to_polygons(p_model1, vrx, l, &NumVrx, &Rhob, &Res, &Inc, &Ints, &AngM, &Suscept/*, fillColour*/);			
		i_net = 0; // index do zapisu oczek na ktore zostaly zrzutowane krawedzie                       
		PolySeg = (struct rasterization*)malloc(NumNet*sizeof(struct rasterization)); // alokacja pamieci na zbior zrasteryzowanych indeksow siatki 			
		if( check_edges_intersection(vrx, NumVrx, l) )
			return 1;
					
		//****************************************************************************************************************************	
		for(k=0; k<NumVrx; k++) // rasteryzacja każdej krawędzi danego wielokata 
		{
			x1 = vrx[k].x;	// wierzcholek 1 
			x2 = vrx[k+1].x;  // wierzcholek 2
			z1 = vrx[k].z;	// wierzcholek 1
			z2 = vrx[k+1].z;	// wierzcholek 2			
														
			if( ( x1 > p_model1->xRange || x1 < 0 ) || ( x2 > p_model1->xRange || x2 < 0 ) || ( z1 > p_model1->zRange || z1 < 0 ) || ( z2 > p_model1->zRange || z2 < 0 )  )
			{									  
				g_print ("\n Pkt. poza obszarem modelu - krawedz: x1=%lf, z1=%lf, x2=%lf, z2=%lf \n",x1,z1,x2,z2);
				return 1;
			}
						
			if(x1 == p_model1->xRange) // zabezpieczenie zeby nie leżał na krawędzi siatki - doraźne  
				x1 = x1 - 0.5;
			if(x2 == p_model1->xRange)
				x2 = x2 - 0.5;		
			if(z1 == p_model1->zRange)
				z1 = z1 - 0.5;	
			if(z2 == p_model1->zRange)
				z2 = z2 - 0.5;	
			// lokalny ukl. wspolrzednych: x1,z1 srodek ukladu
			
			if( x1 < x2 && z1 >= z2 )  // 0 (360) <= alf < 90 ..... cwiartka 1
				flag2=1;
			if( x1 >= x2 &&  z1 > z2 )  // 90 <= alf < 180 ..... cwiartka 2
				flag2=2;
			if( x1 > x2 &&  z1 <= z2 )  // 180 <= alf < 270 ..... cwiartka 3
				flag2=3;
			if( x1 <= x2 && z1 < z2 )  // 270 <= alf < 360 ..... cwiartka 4
				flag2=4;      
					                                                                 
			for(i=0; i<NumNet; i++) // petla przez cala siatke							
			{	
				if (p_net[i].x2> x1 && x1 >= p_net[i].x1 && p_net[i].z2 > z1 && z1>=p_net[i].z1 ) // znalezienie oczka w ktorym znajduje sie poczatkowy wierzcholek  	
				{   
					p_net[i].rhob = Rhob;
					p_net[i].res = Res;
					p_net[i].inc  = Inc;
					p_net[i].ints = Ints;
					p_net[i].angm = AngM;
					p_net[i].suscept = Suscept;
					p_net[i].flag = flag2;
											
					PolySeg[i_net].index= i;
					PolySeg[i_net].flag = flag2;	
					i_net++;
					
					break;
				}
			}
			
			switch(flag2) // wybranie w zaleznosci od cwiartki w ktorej znajduje sie krawedz i wykonanie po calej krawedzi 
			{
				case 1:					
					while( (0 == (p_net[i].x2 > x2 && x2 >= p_net[i].x1 && p_net[i].z2>z2 && z2>=p_net[i].z1)) ) // wykonywanie az do konca krawedzi   
					{				
						a = ( (p_net[i].z1-z2)-(p_net[i].z1-z1) ) / ( (x2-(p_net[i].x2))-(x1-(p_net[i].x2)) );
						b = (p_net[i].z1-z1) - a*(x1-(p_net[i].x2));
																			
						if(b > 0)
						{
							i = i-NumNetX;
							p_net[i].rhob = Rhob;
							p_net[i].res = Res;
							p_net[i].inc = Inc;
							p_net[i].ints = Ints;
							p_net[i].angm = AngM;
							p_net[i].suscept = Suscept;
							p_net[i].flag = flag2;											
							PolySeg[i_net].index = i;
							PolySeg[i_net].flag = flag2;	
							i_net++;							
						}	
		       			else if (b < 0)
						{
							i = i+1;
							p_net[i].rhob = Rhob;
							p_net[i].res = Res;
							p_net[i].inc = Inc;
							p_net[i].ints = Ints;
							p_net[i].angm = AngM;
							p_net[i].suscept = Suscept;
							p_net[i].flag = flag2;
							PolySeg[i_net].index = i;
							PolySeg[i_net].flag = flag2;	
							i_net++;							
						}
						else if (b == 0)					
						{
							i = i + 1;
							p_net[i].rhob = Rhob;
							p_net[i].res = Res;
							p_net[i].inc = Inc;
							p_net[i].ints = Ints;
							p_net[i].angm = AngM;
							p_net[i].suscept = Suscept;
							p_net[i].flag = flag2;
							PolySeg[i_net].index = i;
							PolySeg[i_net].flag = flag2;	
							i_net++;							
						}
					}			
				break; // cala jedna krawedz zrasteryzowana
	
				case 2:
						while( (0 == (p_net[i].x2 > x2 && x2 >= p_net[i].x1 && p_net[i].z2>z2 && z2 >= p_net[i].z1)) )   // wykonywanie az do konca krawedzi   
						{																					
							if(x2 == x1)
							{
								i = i - NumNetX;	
								p_net[i].rhob = Rhob;
								p_net[i].res = Res;
								p_net[i].inc = Inc;
								p_net[i].ints = Ints;
								p_net[i].angm = AngM;
								p_net[i].suscept = Suscept; 
								p_net[i].flag = flag2;
								PolySeg[i_net].index = i;
								PolySeg[i_net].flag = flag2;	
								i_net++;						
							}
						   	else 	
							{
								a = ((p_net[i].z1-z2) - (p_net[i].z1-z1) ) / ( (x2-p_net[i].x1) - (x1-p_net[i].x1) );
								b =  (p_net[i].z1-z1) - a*(x1-p_net[i].x1);															
									
								if(b > 0)
								{
									i = i - NumNetX;	
									p_net[i].rhob = Rhob;
									p_net[i].res = Res;
									p_net[i].inc = Inc;
									p_net[i].ints = Ints;
									p_net[i].angm = AngM;
									p_net[i].suscept = Suscept;
									p_net[i].flag = flag2;
									PolySeg[i_net].index = i;
									PolySeg[i_net].flag = flag2;	
									i_net++;									
								}
								else if(b == 0)
								{
									i = i-NumNetX-1;	
									p_net[i].rhob = Rhob;
									p_net[i].res = Res;
									p_net[i].inc = Inc;
									p_net[i].ints= Ints;
									p_net[i].angm = AngM;
									p_net[i].suscept=Suscept;
									p_net[i].flag = flag2;
									PolySeg[i_net].index = i;
									PolySeg[i_net].flag = flag2;	
									i_net++;
								}																				
								else 
								{
									i = i - 1;	
									p_net[i].rhob = Rhob;
									p_net[i].res = Res;
									p_net[i].inc = Inc;
									p_net[i].ints= Ints;
									p_net[i].angm = AngM;
									p_net[i].suscept = Suscept; 
									p_net[i].flag = flag2;
									PolySeg[i_net].index = i;
									PolySeg[i_net].flag = flag2;	
									i_net++;									
								}
							}													
						}						
				break; // cala jedna krawedz zrasteryzowana

				case 3:
						while( (0 == (p_net[i].x2 > x2 && x2 >= p_net[i].x1 && p_net[i].z2 >= z2 && z2>=p_net[i].z1)) )   // wykonywanie az do konca krawedzi   
						{							
								a = ((p_net[i].z2-z2) - (p_net[i].z2-z1) ) / ( (x2-p_net[i].x1) - (x1-p_net[i].x1) );
								b =  (p_net[i].z2-z1) - a*(x1-p_net[i].x1);
											
								if(b > 0)
								{
									i = i - 1;	
									p_net[i].rhob = Rhob;
									p_net[i].res = Res;
									p_net[i].inc = Inc;
									p_net[i].ints = Ints;
									p_net[i].angm = AngM;
									p_net[i].suscept=Suscept;
									p_net[i].flag = flag2;
									PolySeg[i_net].index = i;
									PolySeg[i_net].flag = flag2;	
									i_net++;									
								}
								else if(b==0)
								{
									i = i + NumNetX - 1;	
									p_net[i].rhob = Rhob;
									p_net[i].res = Res;
									p_net[i].inc = Inc;
									p_net[i].ints = Ints;
									p_net[i].angm = AngM;
									p_net[i].suscept=Suscept; 
									p_net[i].flag = flag2;
									PolySeg[i_net].index = i;
									PolySeg[i_net].flag = flag2;	
									i_net++;									
								}
								else 
								{
									i = i + NumNetX;	
									p_net[i].rhob = Rhob;
									p_net[i].res = Res;
									p_net[i].inc = Inc;
									p_net[i].ints = Ints;
									p_net[i].angm = AngM;
									p_net[i].suscept=Suscept;
									p_net[i].flag = flag2;
									PolySeg[i_net].index = i;
									PolySeg[i_net].flag = flag2;	
									i_net++;									
								}
						}
						//w przypadku gdy nie doszlo do  koncowego wierzcholka krawedzi (pkt wypadl na wezle)
												
						if(0 == ( p_net[PolySeg[i_net-1].index].x2 > x1 && x1 >= p_net[PolySeg[i_net-1].index].x1 && p_net[PolySeg[i_net-1].index].z2 > z1 && z1 >= p_net[PolySeg[i_net-1].index].z1 ))					
							if( ((i+NumNetX) < NumNet) && (p_net[i+NumNetX].x1 == x2 ) && (p_net[i+NumNetX].z1 == z2 ) )
							{
								i = i + NumNetX;	
								p_net[i].rhob = Rhob;
								p_net[i].res = Res;
								p_net[i].inc = Inc;
								p_net[i].ints= Ints;
								p_net[i].angm = AngM;
								p_net[i].suscept=Suscept;
								p_net[i].flag = flag2;
								PolySeg[i_net].index= i;
								PolySeg[i_net].flag = flag2;	
								i_net++;	
							}										
				break; // cala jedna krawedz zrasteryzowana
				
				case 4:
						while( (0 == (p_net[i].x2 > x2 && x2 >= p_net[i].x1 && p_net[i].z2 > z2 && z2 >= p_net[i].z1)) )   // wykonywanie az do konca krawedzi   
						{										
							if(x2 == x1)		
							{
								i = i + NumNetX;	
								p_net[i].rhob = Rhob;
								p_net[i].res = Res;
								p_net[i].inc = Inc;
								p_net[i].ints = Ints;
								p_net[i].angm = AngM;
								p_net[i].suscept=Suscept;
								p_net[i].flag = flag2;
								PolySeg[i_net].index = i;
								PolySeg[i_net].flag = flag2;	
								i_net++;								
							}
							else
							{				
								a = ( (p_net[i].z2-z2)-(p_net[i].z2-z1) ) / ( (x2-(p_net[i].x2)) - (x1-(p_net[i].x2)) );
								b = (p_net[i].z2-z1) - a*(x1 - (p_net[i].x2));
																			
								if(b > 0)
								{
									i = i + 1;	
									p_net[i].rhob = Rhob;
									p_net[i].res = Res;
									p_net[i].inc = Inc;
									p_net[i].ints= Ints;
									p_net[i].angm = AngM;
									p_net[i].suscept=Suscept;
									p_net[i].flag = flag2;
									PolySeg[i_net].index= i;
									PolySeg[i_net].flag = flag2;	
									i_net++;									
								}									
								else if(b==0)
								{
									i = i + NumNetX + 1;	
									p_net[i].rhob = Rhob;
									p_net[i].res = Res;
									p_net[i].inc = Inc;
									p_net[i].ints = Ints;
									p_net[i].angm = AngM;
									p_net[i].suscept = Suscept; 
									p_net[i].flag = flag2;
									PolySeg[i_net].index = i;
									PolySeg[i_net].flag = flag2;	
									i_net++;									
								}								
								else 
								{
									i = i + NumNetX;	
									p_net[i].rhob = Rhob;
									p_net[i].res = Res;
									p_net[i].inc = Inc;
									p_net[i].ints = Ints;
									p_net[i].angm = AngM;
									p_net[i].suscept = Suscept;
									p_net[i].flag = flag2;
									PolySeg[i_net].index = i;
									PolySeg[i_net].flag = flag2;	
									i_net++;									
								}
							}	
					}
				break; // cala jedna krawedz zrasteryzowana				 
			} // koniec switch i rasteryzacji jednej krawedzi
	
		} // krawedzie wielokata zrasteryzowane
		//****************************************************************************************************************************					
		PolySeg[i_net] = PolySeg[0]; 		
		
		for(i=0; i<i_net; i++) // petla przez wsztystkie zrasteryzowane oczka 
		{
			if( p_net[PolySeg[i].index].z1 == p_net[PolySeg[i+1].index].z1 )
			{
				if((PolySeg[i].flag == PolySeg[i+1].flag)			  || (PolySeg[i].flag == 1 && PolySeg[i+1].flag == 2 ) ||  // redukuje 2 pkt na laczeniu 
					(PolySeg[i].flag == 2 && PolySeg[i+1].flag == 1 ) || (PolySeg[i].flag == 2 && PolySeg[i+1].flag == 2 ) ||
					(PolySeg[i].flag == 3 && PolySeg[i+1].flag == 4 ) || (PolySeg[i].flag == 4 && PolySeg[i+1].flag == 3 ) )
				{	
					PolySeg[i].flag = -1; 	
				}			
			}
		}	
 
		CommPolyZ = (int*)malloc(i_net*sizeof(int)); // zbior indeksow o takiej samej wsp. Z 

        for(i=0; i<i_net-1; i++) // sortowanie babelkowe zamarkowanych punktow ze wzgledu na Z  
		{
			for(n=0; n<i_net-1-i; n++)
			{
				if(p_net[PolySeg[n].index].z1 > p_net[PolySeg[n+1].index].z1)
				{
					tempPolySeg  = PolySeg[n+1];
					PolySeg[n+1] = PolySeg[n];
					PolySeg[n]  = tempPolySeg;
				}
			}
		}	 
		    
		for(j=0; j<i_net-1; ) // wypelnianie wielokata  
		{
			k = 0; // k liczba punktow wspolnych 			
			 															
			while( (j<i_net-1) && (p_net[PolySeg[j].index].z1 == p_net[PolySeg[j+1].index].z1)) // znalezenie wszystkich zamarkowanych oczek o wspolnym Z i przypisanie do CommPolyZ, w kolejnych petlach nizsze Z
			{				
				CommPolyZ[k] = j;
				j++;
				k++;							
			}		

			CommPolyZ[k] = j;			
			j++;
			k++;																			
					
			for(i=0; i<k-1; i++) // sortowanie babelkowe ze wzgledu na X dla pkt. o tej samej wsp. Z 
			{
				for(n=0; n<k-1-i; n++)
				{
					if(p_net[PolySeg[CommPolyZ[n]].index].x1 > p_net[PolySeg[CommPolyZ[n+1]].index].x1)
					{
						tempPolySeg  = PolySeg[CommPolyZ[n+1]];
						PolySeg[CommPolyZ[n+1]] = PolySeg[CommPolyZ[n]];
						PolySeg[CommPolyZ[n]]  = tempPolySeg;
					}
				}
			}		
												
			for(i=0; i<k; i++) // przechodzenie po lini od X min do X max
			{
				if( PolySeg[CommPolyZ[i]].flag != -1 ) // znalezenie pkt. liczacego sie (zamarkowanego) 
				{	
					for(n=i+1; n<k; n++)
					{
						if(PolySeg[CommPolyZ[n]].flag != -1) // znalezenie 2-giego pkt. liczacego sie
						{									
							while (PolySeg[CommPolyZ[i]].index != PolySeg[CommPolyZ[n]].index) //wypelnienie obszeru pomiedzy 2 pkt. liczacymi
							{
								p_net[PolySeg[CommPolyZ[i]].index + 1].rhob = p_net[PolySeg[CommPolyZ[i]].index].rhob;
								p_net[PolySeg[CommPolyZ[i]].index + 1].res = p_net[PolySeg[CommPolyZ[i]].index].res;
								p_net[PolySeg[CommPolyZ[i]].index + 1].inc = p_net[PolySeg[CommPolyZ[i]].index].inc;
								p_net[PolySeg[CommPolyZ[i]].index + 1].ints = p_net[PolySeg[CommPolyZ[i]].index].ints;
								p_net[PolySeg[CommPolyZ[i]].index + 1].angm = p_net[PolySeg[CommPolyZ[i]].index].angm;
								p_net[PolySeg[CommPolyZ[i]].index + 1].suscept = p_net[PolySeg[CommPolyZ[i]].index].suscept;
																								
								PolySeg[CommPolyZ[i]].index ++;
							}
							i = n;
							break;
						}
					}
				}
			}		
		}      	
		g_free(CommPolyZ);
        g_free(PolySeg);
		g_free(vrx); 
	}  // koniec rasteryzacji wielokata	       
	
	l = NumNet - 1;
	 
  	for(k=0; k<NumNet; k += NumNetX ) // extend grid, left side 
	{
	    l++;
		p_net[l].rhob = p_net[k].rhob; 
		p_net[l].res  = p_net[k].res; 
		p_net[l].inc = p_net[k].inc; 
		p_net[l].ints  = p_net[k].ints; 
		p_net[l].angm  = p_net[k].angm; 
		p_net[l].suscept = p_net[k].suscept;
	}

	for(k=NumNetX-1; k<=NumNet; k+= NumNetX) // extend grid, right side
	{
		l++;
        p_net[l].rhob = p_net[k].rhob; 
		p_net[l].res  = p_net[k].res;
		p_net[l].inc  = p_net[k].inc;
		p_net[l].ints  = p_net[k].ints;
		p_net[l].angm  = p_net[k].angm;
		p_net[l].suscept = p_net[k].suscept;
	} 

  	return 0; 
} 

int fun_max(double a, double b)
{
	if (a > b)
		return a;
	else
		return b;
} 

int fun_min(double a, double b)
{
	if (a < b)
		return a;
	else 
		return b;
}

gpointer g_list_id_get_data(model *modeltemp, intptr_t id, int flag) // 1 - vertex, 2 - boundary, 3 - layer
{
	int i;
	switch(flag) 
	{
		case 1:	
				for(i = 0; i < g_slist_length(modeltemp->listOfVertex); i++)
				{	
					if( (((vertex*)(g_slist_nth_data(modeltemp->listOfVertex, i)))->id) == id ) 
						return  (g_slist_nth_data(modeltemp->listOfVertex, i)); 					
				}
				break;
		case 2:		
				for(i = 0; i < g_slist_length(modeltemp->listOfBoundaries); i++)
				{
					if( (((boundary*)(g_slist_nth_data(modeltemp->listOfBoundaries, i)))->id) == id ) 
						return  (g_slist_nth_data(modeltemp->listOfBoundaries, i)); 
				}
				break;
		case 3:	
					
				for(i = 0; i < g_slist_length(modeltemp->listOfLayers); i++)
				{
					if( (((layer*)(g_slist_nth_data(modeltemp->listOfLayers, i)))->id) == id ) 
						return  (g_slist_nth_data(modeltemp->listOfLayers, i)); 
				}
				break;		
	}		
	g_print("************* Brak elementu %d na liscie z %d (1-vertex, 2-boundary, 3-layer) ****************** \n", (int)id, flag);
	
	return NULL; 
}

// function to get count of vertex
int get_points_count(model *modeltemp) 
{
	int PointsCount = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	int BoundaryId;
	int VertexId;

	layer* layertemp; 	
	
	for(i=0; i<g_slist_length(modeltemp->listOfLayersId); i++)  // loop for layers
	{		
		layertemp = (layer*)(g_list_id_get_data(modeltemp, (intptr_t)(g_slist_nth(modeltemp->listOfLayersId, (i))->data), 3));
		
		for(j=0; j<g_slist_length(layertemp->listOfBoundariesId); j++)  // loop for boundries
		{
			BoundaryId = (intptr_t)(g_slist_nth(layertemp->listOfBoundariesId, j)->data); 
			((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->used = 0; // clear used marks	
			
			for(k=0; k<g_slist_length( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId); k++)
			{
				VertexId = (intptr_t)(g_slist_nth( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2))) ->listOfVerticesId, k)->data);						
				((vertex*)(g_list_id_get_data(modeltemp,VertexId,1)))->used = 0;			
			}
		}		
		layertemp = NULL;
	}
	
	for(i=0; i<g_slist_length(modeltemp->listOfLayersId); i++) // loop for layers
	{		
		layertemp = (layer*)(g_list_id_get_data(modeltemp, (intptr_t)(g_slist_nth(modeltemp->listOfLayersId, (i))->data), 3));
		
		for(j=0; j<g_slist_length(layertemp->listOfBoundariesId); j++) // loop for boundries
		{		
			BoundaryId = (intptr_t)(g_slist_nth(layertemp->listOfBoundariesId, j)->data); 	
			
			if( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->fixed == 0) // check if boundary can be modify
			{				
				if( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->used == 0) // check if boundary was used before
				{
					for(k=0; k<g_slist_length( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId); k++)
					{
						VertexId = (intptr_t)(g_slist_nth( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId, k)->data);
					
					  	if( ((vertex*)(g_list_id_get_data(modeltemp,VertexId,1)))->used == 0 )
					    {	
							if( ((vertex*)(g_list_id_get_data(modeltemp,VertexId,1)))->fixed_x == 0)
								PointsCount ++; 
						  	if( ((vertex*)(g_list_id_get_data(modeltemp,VertexId,1)))->fixed_z == 0)
								PointsCount ++;		
						   
							((vertex*)(g_list_id_get_data(modeltemp,VertexId,1)))->used = 1;  
					    }														    
				    }
				    
					((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2))) ->used = 1; // mark used boundary
				}
			}			
		}	
		layertemp = NULL;		
	}

	for(i=0; i<g_slist_length(modeltemp->listOfLayersId); i++)  // loop for layers
	{
		layertemp = (layer*)(g_list_id_get_data(modeltemp, (intptr_t)(g_slist_nth(modeltemp->listOfLayersId, (i))->data), 3));
		
		for(j=0; j<g_slist_length(layertemp->listOfBoundariesId); j++)  // loop for boundries
		{
			BoundaryId = (intptr_t)(g_slist_nth(layertemp->listOfBoundariesId, j)->data); 
			((boundary*)(g_list_id_get_data(modeltemp,BoundaryId, 2))) ->used = 0; // clear used marks	
			
			for(k=0; k<g_slist_length( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId, 2)))->listOfVerticesId); k++)
			{
				VertexId = (intptr_t)(g_slist_nth( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2))) ->listOfVerticesId, k)->data);						
				((vertex*)(g_list_id_get_data(modeltemp,VertexId,1)))->used = 0;
			}
		}
		layertemp = NULL;			
	}

	// gdy do inwersji włączone zostaja param. warstwty (opornosc + gestosc)
	PointsCount += 2*g_slist_length(modeltemp->listOfLayersId);

	return PointsCount;
}

int if_belongs(vertex x, vertex y, vertex z)
{
	long det; // wyznacznik macierzy
	det = det_matrix(x,y,z);	
	if(det != 0) // det == 0 collinear points
		return 1; 
	else
	{
		if((fun_min(x.x,y.x)<=z.x)&&(z.x<=fun_max(x.x,y.x)) && (fun_min(x.z,y.z)<=z.z)&&(z.z<=fun_max(x.z,y.z)))
			return 0; 
		else
			return 1; // wierzcholek nie lezy na krawedzi
	}	
}

int layer_vertices_counter(model *modeltemp, int LayerId) 
{
	int i = 0;
	int counter = 0;
	int BoundaryId = 0;

	for(i=0; i < g_slist_length( ((layer*)(g_list_id_get_data(modeltemp,LayerId,3)))->listOfBoundariesId); i++) // ilosc granic 
	{
		BoundaryId = (intptr_t)(g_slist_nth_data(( ((layer*)(g_list_id_get_data(modeltemp,LayerId,3)))->listOfBoundariesId),i));
		counter += g_slist_length ( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2))) -> listOfVerticesId );   
	}

	return(counter);
}

// Funkcja przepisujaca warstwe na poligon. 
// Parametry:
// int l        - nr warstwy
// vertex *vrx  - struktura przechowujaca wspolrzedne wszystkich wierzcholkow w warstwie
// int p_NumVrx - ilosc wierzcholkow w warstwie
vertex *layers_to_polygons(model *modeltemp, vertex *vrx, int l, int *p_NumVrx, double *p_Rhob, double *p_Res, double *p_Inc, double *p_Ints, double *p_AngM, double *p_Suscept/*, colour *fillColour*/)	
{
	polygon *polygontemp;
	
	int i = 0;
	int j = 0;
	int LayerId    = 0;
	int BoundaryId = 0;
	int VertexId   = 0;
	
	// ******************** Przepisanie warstwy na struct polygon //***************************************************************
		VertexId = 0;
		polygontemp= malloc(sizeof(polygon)); 
		polygontemp->listOfVerticesId=NULL;	// struktura na ktora przepisujemy 	
		
		LayerId = (intptr_t)(g_slist_nth_data(modeltemp->listOfLayersId,l));		
		polygontemp->density = ((layer*)(g_list_id_get_data(modeltemp,LayerId,3)))->density;
		polygontemp->resistivity = ((layer*)(g_list_id_get_data(modeltemp,LayerId,3)))->resistivity;
		polygontemp->inclination = ((layer*)(g_list_id_get_data(modeltemp,LayerId,3)))->inclination;
		polygontemp->intensity = ((layer*)(g_list_id_get_data(modeltemp,LayerId,3)))->intensity;
		polygontemp->angle_mag = ((layer*)(g_list_id_get_data(modeltemp,LayerId,3)))->angle_mag;
		polygontemp->susceptibility=((layer*)(g_list_id_get_data(modeltemp,LayerId,3)))->susceptibility;		

		for(j=0; j<g_slist_length( ((layer*)(g_list_id_get_data(modeltemp,LayerId,3))) -> listOfBoundariesId ); j++) // pętla po wszystkich granicach 
		{
			BoundaryId = (intptr_t)(g_slist_nth( (((layer*)(g_list_id_get_data(modeltemp,LayerId,3)))->listOfBoundariesId), j )->data);
	
			if(j == 0) // dla pierwszej granicy 
			{
				for(i=0; i<g_slist_length(((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2))) -> listOfVerticesId ); i++) // pętla po wszystkich wierzcholkach
				{
					VertexId   = (intptr_t)(g_slist_nth( (((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId ), i )->data);
					polygontemp->listOfVerticesId = g_slist_append (polygontemp->listOfVerticesId, GINT_TO_POINTER(VertexId));
				}
		    }
		    else 
			{	/*
				przepisanie wierzcholkow z granicy w odpowiedniej kolejnosci zaleznej od ostatniego wierzcholka w poprzedniej granicy tzn.,
				jezeli ostatni wierzcholek z poprzedzajacej granicy pokrywa sie z pierwszym wierzcholki wczytywanej pobieranie w normalniej kolejnosci 
				jezeli ostatni wierzcholek z poprzedzajacej granicy pokrywa sie z ostatnim wierzcholki wczytywanej pobieranie w odwroconej kolejnosci 				  
				*/
				if( VertexId == (intptr_t)(g_slist_nth( (((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId ), 0 )->data) ) 
			    {
					for(i=0; i<g_slist_length( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2))) -> listOfVerticesId ); i++) // pętla po wszystkich wierzcholkach
					{		
						if( VertexId != (intptr_t)(g_slist_nth( (((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId ), i )->data)  )	
						{
							VertexId = (intptr_t)(g_slist_nth( (((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId ), i )->data);				
							polygontemp->listOfVerticesId = g_slist_append (polygontemp->listOfVerticesId, GINT_TO_POINTER(VertexId)); 	
						}
					}
				}
				else if( VertexId == (intptr_t)(g_slist_last( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId )->data) )
				{
					for(i=g_slist_length( ((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2))) -> listOfVerticesId )-1; i>=0; i--) // pętla po wszystkich wierzcholkach,odwrocona kolejnosc
					{		
						if( VertexId != (intptr_t)(g_slist_nth( (((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId ), i )->data)  )	
						{
							VertexId =(intptr_t)(g_slist_nth( (((boundary*)(g_list_id_get_data(modeltemp,BoundaryId,2)))->listOfVerticesId ), i )->data);				
							polygontemp->listOfVerticesId = g_slist_append(polygontemp->listOfVerticesId, GINT_TO_POINTER (VertexId)); 	
						}
					}
				}
				else
				{
					g_print("\n\n ****Error: MainForm %d. LayerID:%d BoundaryID:%d VertexID:%d \n\n", __LINE__, LayerId,BoundaryId,VertexId); 
				}
			}
		} 
	// ******************** KONIEC Przepisanie warstwy na struct polygon //***************************************************************
 
	*p_NumVrx = g_slist_length(polygontemp->listOfVerticesId);		
	
	*p_Rhob = polygontemp->density;
	*p_Res  = polygontemp->resistivity;	
	*p_Inc = polygontemp->inclination;
	*p_Ints = polygontemp->intensity;
	*p_AngM = polygontemp->angle_mag;
	*p_Suscept = polygontemp->susceptibility;
						
	vrx = (vertex*)malloc((*p_NumVrx)*sizeof(vertex)); // alokacja pamieci dla wierzcholkow 	  
								
	for(i=0, j=0; i<(*p_NumVrx); i++, j++) // pobranie wsp. wierzcholkow
	{
		vrx[j].x = ((vertex*) (g_list_id_get_data(modeltemp, (intptr_t)(g_slist_nth_data(polygontemp->listOfVerticesId,i)),1)) )->x;
		vrx[j].z = ((vertex*) (g_list_id_get_data(modeltemp, (intptr_t)(g_slist_nth_data(polygontemp->listOfVerticesId,i)),1)) )->z;	
	}		

	*p_NumVrx = j - 1;

	g_free(polygontemp); 
	
	return vrx;
}

block *net_maker(model *p_model1, int isRegular, block *p_net, GSList *p_ListOf_net_z_range)
{
	//FILE *net_data;
	int i=0; 
	int j=0;
	int l=0;
	int k=0;       
	int num_horiz=0;  // num_horiz*2+1 = liczba wezlow po x 
	int num_vertic=0; // liczba wezlow po y  
	int size=0; 
	long int extend=0; 
	
	int temp_mid = p_model1->widMidNod;   // szer. srodkowego oczka siatki
	int temp_first = p_model1->thick1Lay; // miazszasc pierwszej warsty   
	double x_base = p_model1->ratioX; 	 // podstawa wykładnika siatki X
	double z_base = p_model1->ratioZ;  	 // podstawa wykładnika siatki X
	double temp_step_x = 0.0; 			 	 // rozrzucenie niedomiaru/nadmiaru po x
	double temp_step_z = 0.0;  				 // rozrzucenie niedomiaru/nadmiaru po z
	double temp_last_brick = 0.0;
	
	if(!isRegular)
	{
		//******************** X  **************************************************************************************************	
			// rozszerzamy siatkę od srodka - jako ze pracujemy na int sprawdzamy parzystosc 
			if( (temp_mid % 2) != 0 )	
				g_print("Error: width of middle brick is odd"); 
			if( (p_model1->xRange%2 ) != 0 )	
				g_print("Error: length of profile is odd"); 

			temp_step_x = p_model1->xRange / 2 + temp_mid / 2; // wart. poczatkowa do wyliczenia sumy   

			for (i=1; temp_step_x < (p_model1->xRange); i++) // zliczenie il. oczek od prawej granicy srodkowego oczka				
				temp_step_x += ceil(temp_mid*pow(x_base,i));  			   
			
			i = i-1;	

			temp_last_brick = temp_step_x - ceil(temp_mid*pow(x_base, i));
			
			if((p_model1->xRange - temp_last_brick) < ceil(temp_mid*pow(x_base,i-1))*0.4 )  // rozrzuci po oczkach nadmiaru gdy zabraknie  < 0.4 szerokosci ostatniego oczka do prawej granicy  
			{ 
				num_horiz = i-1; // -1 zwiazane z działaniem petli for 		 
				temp_step_x = 0;
				for(i=1; i<=num_horiz; i++) 
					temp_step_x += pow(x_base,i);
				temp_step_x = (p_model1->xRange - temp_last_brick) / temp_step_x; // nadmiar zalezny od pow(x_base,i)
			}
			else  // dodatnie po X dodatkowego oczka i rozrzucenie niedomiaru 
			{
				num_horiz = i; 				   // bez -1, dodanie oczka
				temp_last_brick = temp_step_x; // przez co zmiana polozenia granicy ostatniego oczka		
				temp_step_x = 0;
				for(i=1; i<=num_horiz; i++)  
					temp_step_x += pow(x_base, i);
				temp_step_x = (p_model1->xRange - temp_last_brick) / temp_step_x; // niedomiar zalezny od pow(z_base,i)
			}

		//******************** Z  **************************************************************************************************
			for(i=0; temp_step_z < (p_model1->zRange); i++)  // wyliczenie max. potegi  
				temp_step_z += ceil(temp_first*pow(z_base, i));
			i = i - 1;
				
			temp_last_brick = temp_step_z - ceil(temp_first*pow(z_base, i));

			// rozrzuci po oczkach nadmiaru gdy zabraknie  < 0.4 szerokosci ostatniego oczka do dolnej granicy
			if( (p_model1->zRange - temp_last_brick) < ceil(temp_first*pow(z_base, i-1))*0.4 )    
			{ 
				num_vertic = i-1; // -1 zwiazane z działaniem petli for 
				
				temp_step_z = 0;
				for(i=1; i<=num_vertic; i++) 
					temp_step_z += pow(z_base, i);

				temp_step_z = (p_model1->zRange - temp_last_brick) / temp_step_z; // nadmiar zalezny od pow(z_base,i)
			}
			// gdy zabraknie >= 0.4 szerokosci ostatniego oczka do dolnej granicy dodatnie po Z dodatkowego oczka i rozrzuci niedomiaru 
			else
			{
				num_vertic = i; // bez -1, dodanie oczka
				temp_last_brick = temp_step_z; // przez co zmiana polozenia granicy ostatniego oczka
				
				temp_step_z = 0;
				for(i=1; i<=num_vertic; i++)  
					temp_step_z += pow(z_base, i);

				temp_step_z = (p_model1->zRange - temp_last_brick) / temp_step_z; // niedomiar zalezny od pow(z_base,i)
			}

		//******************************************************************************* implementacja siatki                        
		l = num_horiz;
		size = (num_vertic + 1)*(2*num_horiz + 1) + 2*(num_vertic + 1);  

		p_net = (block*)malloc(size * sizeof(block)); // alokacja siatki
	
		if(p_net==NULL) 
		{ 
			g_print("\n 1 Przydzielenie pamięci dla siatki nie było możliwe"); 
			getchar(); 
		} 		

		p_net[l].x1 = (p_model1->xRange / 2.0) - (temp_mid / 2.0);                    
		p_net[l].x2 = p_net[l].x1 + temp_mid; 
		p_net[l].z1 = 0;
		p_net[l].z2 = temp_first;
		p_net[l].rhob = p_model1->density;
		p_net[l].res = p_model1->resistivity;
		p_net[l].inc = p_model1->inclination;
		p_net[l].ints = p_model1->intensity;
		p_net[l].angm = p_model1->angle_mag;
		p_net[l].suscept = p_model1->susceptibility;

		if(x_base != 1)
		{
			for(j=1; j<=num_vertic+1; j++) 
			{
				for(i=1; i<=num_horiz; i++)
				{				
					// strona lewa  		
					p_net[l-i].x2 = p_net[l-i+1].x1;                  
					p_net[l-i].x1 = p_net[l-i].x2 - ceil( (temp_step_x + temp_mid) * pow(x_base, i));
					if(i == num_horiz) // dociagniecie do granicy lewej 
						p_net[l-i].x1 = 0;								
					p_net[l-i].z1 = p_net[l].z1;
					p_net[l-i].z2 = p_net[l].z2;
					p_net[l-i].rhob = p_model1->density;
					p_net[l-i].res = p_model1->resistivity;
					p_net[l-i].inc = p_model1->inclination;
					p_net[l-i].ints = p_model1->intensity;
					p_net[l-i].angm = p_model1->angle_mag;
					p_net[l-i].suscept = p_model1->susceptibility;
					p_net[l-i].flag = -1;
				
					// strona prawa
					p_net[l+i].x1 = p_net[l+i-1].x2;                   
					p_net[l+i].x2 = p_net[l+i].x1 + ceil ( (temp_step_x+temp_mid) * pow(x_base,i));
					if(i == num_horiz) // dociagniecie do prawej granicy				  
						p_net[l+i].x2 = p_model1->xRange;				     
					p_net[l+i].z1 = p_net[l].z1;
					p_net[l+i].z2 = p_net[l].z2;
					p_net[l+i].rhob = p_model1->density;
					p_net[l+i].res = p_model1->resistivity;
					p_net[l+i].inc = p_model1->inclination;
					p_net[l+i].ints = p_model1->intensity;
					p_net[l+i].angm = p_model1->angle_mag;
					p_net[l+i].suscept = p_model1->susceptibility;
					p_net[l+i].flag = -1;
				}
				
				l += 2*num_horiz + 1; // przygotowanie dla następnego rzedu
							
				p_net[l].x1 = (p_model1->xRange / 2.0) - (temp_mid / 2.0);                   
				p_net[l].x2 = p_net[l].x1 + temp_mid;
				p_net[l].z1 = p_net[l-2*num_horiz].z2;
				p_net[l].z2 = p_net[l].z1 + ceil ( (temp_first+temp_step_z) * pow(z_base,j) );

				if(j == (num_vertic)) // dociagniecie do granicy ostatniego oczka  
					p_net[l].z2 = p_model1->zRange;
			
				p_net[l].rhob = p_model1->density;
				p_net[l].res = p_model1->resistivity;
				p_net[l].inc = p_model1->inclination;
				p_net[l].ints = p_model1->intensity;
				p_net[l].angm = p_model1->angle_mag;
				p_net[l].suscept = p_model1->susceptibility;
				p_net[l].flag = -1;
			}
						
			p_model1->NumNetZ = j - 1; 
			p_model1->NumNetX = 2*num_horiz + 1; 	   
		}
		else 
		{ 		
			p_net[0].x1 = 0;                    
			p_net[0].x2 = temp_mid; 
			p_net[0].z1 = 0;
			p_net[0].z2 = temp_first;
			p_net[0].rhob = p_model1->density;
			p_net[0].res = p_model1->resistivity;
			p_net[0].inc = p_model1->inclination;
			p_net[0].ints = p_model1->intensity;
			p_net[0].angm = p_model1->angle_mag;
			p_net[0].suscept = p_model1->susceptibility;

			for(l=1, j=1; j<=num_vertic+1 ; j++)  // next row
			{
				for(i=0; p_net[l-1].x2<p_model1->xRange; i++) // next column
				{
					p_net[l].x1 = p_net[l-1].x2;
					p_net[l].z1 = p_net[l-1].z1;    
					p_net[l].x2 = p_net[l].x1 + temp_mid;
					p_net[l].z2 = p_net[l-1].z2; 	                   
					p_net[l].rhob = p_model1->density;
					p_net[l].res = p_model1->resistivity;
					p_net[l].inc = p_model1->inclination;
					p_net[l].ints = p_model1->intensity;
					p_net[l].angm = p_model1->angle_mag;
					p_net[l].suscept = p_model1->susceptibility;
					p_net[l].flag = -1;                                      
					l++;
				}	

				p_net[l].x1 = 0;
				p_net[l].z1 = p_net[l-1].z2;    
				p_net[l].x2 = temp_mid;
				p_net[l].z2 = p_net[l].z1 + ceil( (temp_first + temp_step_z) * pow(z_base, j) );	                   
				p_net[l].rhob = p_model1->density;
				p_net[l].res = p_model1->resistivity;
				p_net[l].inc = p_model1->inclination;
				p_net[l].ints = p_model1->intensity;
				p_net[l].angm = p_model1->angle_mag;
				p_net[l].suscept = p_model1->susceptibility;
				p_net[l].flag = -1;
				l++;                            
			}	

			p_model1->NumNetZ = j - 1; 
			p_model1->NumNetX = i + 1; 		
		}
			
	} // end irregular==true
	else // siatka regularna
	{   	
		if (g_slist_length(p_ListOf_net_z_range) != 0)
		{    
			temp_step_z     = ((net_z_intervals*)(g_slist_nth_data(p_ListOf_net_z_range,0)))->z_step;
			temp_last_brick = ((net_z_intervals*)(g_slist_nth_data(p_ListOf_net_z_range,0)))->depth;
			

			for(j=0,i=0,num_vertic=0; j<p_model1->zRange; j+= temp_step_z, num_vertic++ ) // obliczenie ilości wezlow
			{
				if(j >= temp_last_brick) // dopoki nie dotrze do nastepnego zakresu
				{ 
					i++;
					if ( i<g_slist_length(p_ListOf_net_z_range) )
					{
						temp_step_z     = ((net_z_intervals*)(g_slist_nth_data(p_ListOf_net_z_range,i)))->z_step;
						temp_last_brick = ((net_z_intervals*)(g_slist_nth_data(p_ListOf_net_z_range,i)))->depth;	  
					}
					else
					{
						temp_step_z	    = temp_first;
						temp_last_brick = p_model1->zRange ; 
					}
				}				  
			}
			
			size = ((p_model1->xRange)/(temp_mid)+1) * num_vertic + 2*num_vertic;  

			p_net=(block*)malloc(size * sizeof(block)); // alokacja siatki

			if (p_net==NULL) 
			{ 
				g_print("\n 2 Przydzielenie pamięci dla siatki nie było możliwe"); 
				getchar(); 
			} 

			temp_step_z     = ((net_z_intervals*)(g_slist_nth_data(p_ListOf_net_z_range,0)))->z_step;
			temp_last_brick = ((net_z_intervals*)(g_slist_nth_data(p_ListOf_net_z_range,0)))->depth;

			for(j=0, l=0; j<p_model1->zRange; j+=temp_step_z) // next row
			{
					if (j >= temp_last_brick)
					{ 
						k++;
						if ( k<g_slist_length(p_ListOf_net_z_range) )
							{
							temp_step_z     = ((net_z_intervals*)(g_slist_nth_data(p_ListOf_net_z_range,k)))->z_step;
							temp_last_brick = ((net_z_intervals*)(g_slist_nth_data(p_ListOf_net_z_range,k)))->depth;	  
							}
						else
							{
								temp_step_z	    = temp_first;
								temp_last_brick = p_model1->zRange; 
							}
					}				
				for (i=0; i<p_model1->xRange; i+=temp_mid) // next column
				{
					p_net[l].x1=i;
					p_net[l].z1=j;    
					p_net[l].x2=p_net[l].x1+temp_mid;
					p_net[l].z2=p_net[l].z1+temp_step_z; 	                   
					p_net[l].rhob= p_model1->density;
					p_net[l].res = p_model1->resistivity;
					p_net[l].inc = p_model1->inclination;
					p_net[l].ints = p_model1->intensity;
					p_net[l].angm = p_model1->angle_mag;
					p_net[l].suscept = p_model1->susceptibility;
					p_net[l].flag=-1;                                       
					l++;
				}
			}
		
			p_model1->NumNetZ = (int) (num_vertic);
			p_model1->NumNetX = (int) (i/temp_mid); 
		}		
		else
		{
			size = ((p_model1->xRange)/(temp_mid)+1) * ((p_model1->zRange)/(temp_first)+1)+ 2*((p_model1->zRange)/(temp_first));  
			p_net=(block*)malloc(size * sizeof(block)); // alokacja siatki

			if (p_net==NULL) 
			{ 
				g_print("\n 3 Przydzielenie pamięci dla siatki nie było możliwe"); 
				getchar(); 
			} 
		
			for(j=0, l=0; j<p_model1->zRange; j+=temp_first) // next row
			{
				for (i=0; i<p_model1->xRange; i+=temp_mid) // next column
				{
					p_net[l].x1=i;
					p_net[l].z1=j;    
					p_net[l].x2=p_net[l].x1+temp_mid;
					p_net[l].z2=p_net[l].z1+temp_first; 	                   
					p_net[l].rhob= p_model1->density;
					p_net[l].res = p_model1->resistivity;
					p_net[l].inc = p_model1->inclination;
					p_net[l].ints = p_model1->intensity;
					p_net[l].angm = p_model1->angle_mag;
					p_net[l].suscept = p_model1->susceptibility;
					p_net[l].flag=-1;                                       
					l++;
				}
			}
			p_model1->NumNetZ = (int) (j/temp_first);
			p_model1->NumNetX = (int) (i/temp_mid); 
		}
	}

	size = p_model1->NumNetZ*p_model1->NumNetX;
	extend = -InfiNet*p_model1->xRange;
	  
	l = size - 1; 
	  
	for(k=0; k<size; k+=p_model1->NumNetX) // extend grid, left side 
	{
	    l++;
      	p_net[l].x1 = extend; 
		p_net[l].z1 = p_net[k].z1; 
		p_net[l].x2 = 0;
		p_net[l].z2 = p_net[k].z2;
		p_net[l].flag = -1;		
	}
	
	extend = InfiNet*p_model1->xRange;
	for(k=p_model1->NumNetX-1; k<=size; k+=p_model1->NumNetX) // extend grid, right side
	{
		l++;
        p_net[l].x1 = p_net[k].x2;
		p_net[l].z1 = p_net[k].z1; 
		p_net[l].x2 = extend+p_model1->xRange; 
		p_net[l].z2 = p_net[k].z2;
		p_net[l].flag = -1;		
	}	

	return p_net;	
}

void read_line(char *p_line, FILE *p_plik)
{	
	int i, j = 0;
	char line[120];
	char *res;

	res = fgets(line, sizeof line, p_plik); 

	if(res == NULL && !feof(p_plik))
	{
		puts(line);
		g_print("Blad odczytu linii\n");
	}
			
	if((strstr(line, ":") - line) > 0) // wartosci sa po dwukropku
	{											
		j = 0;
		for(i=(strstr(line, ":")-line)+1; i<strlen(line); i++) // nadpisz stringa od pozycji ":" +1
		{		
			line[j] = line[i];										
			j++;
		}		
	
		// wyczysc reszte stringa			
		for(j=j; j<strlen(line); j++)
		{
			line[j] = ' ';
		}	

		for(i=0; i<strlen(line); i++)
			p_line[i] = line[i];			
	}  

	return;		
}

void read_line_id(char *p_line)
{
	int i, j = 0;

	char line[120];

	for(i=0; i<strlen(p_line); i++)
		line[i] = p_line[i];
				
	if((strstr(line, ":") - line) > 0) // wartosci sa po dwukropku
	{											
		j = 0;
		for(i=(strstr(line, ":")-line)+1; i<strlen(line); i++) // nadpisz stringa od pozycji ":" +1
		{		
			line[j] = line[i];										
			j++;
		}		
	
		// wyczysc reszte stringa			
		for(j=j; j<strlen(line); j++)
		{
			line[j] = ' ';
		}	

		for(i=0; i<strlen(p_line); i++)
			p_line[i] = line[i];			
	}  

	return;	
}

int sign(long a)
{    
	if(a >= 0)    		
		return 0;
	else // if(a < 0)    	
		return -1;    	
}
// End functions from ModelRaster

double anorm(int n, double *p_d)
{	
	// RETURNS THE SQUARE OF THE EUCLIDEAN NORM OF A VECTOR
	      
	int i;
	double res = 0.0;

	for(i=0; i<n; i++)
	{            
        res += p_d[i] * p_d[i];
	}	      
        
    return res;   
}

void calculate_mag(block *p_net, model* calc_model, double **p_MG_results_temp, int p_MG_sites, double **p_MG_data)
{	
	int i, j, k;
	int sig = -1;

	int NumNet  = (calc_model->NumNetX)*(calc_model->NumNetZ); // ilosc oczek siatki
	int NumNetZ = calc_model->NumNetZ; 
	
	double x0, x1, x2, z1, z2, r1sq, r2sq, term1, term2, sum, denom, mg, r3sq, r4sq;
	double incc, term3, term4, theta1, theta2, theta3, theta4, Fe, r1, r2, r3, r4;
	
	Fe=58000;
	incc = 75*PI/180;
	
	if(p_MG_results_temp == NULL) 
	{ 
		g_print("\n Magneto: Przydzielenie pamieci dla wynikow nie bylo mozliwe"); 
		getchar(); 
	} 
	
	for(j=0; j<p_MG_sites; j++) // obliczanie efektu dla pojedynczego punktu pomiarowego
	{
		x0 = p_MG_data[j][0];				
		mg = 0;

		for(i=0; i<NumNet+2*NumNetZ; i++) // petla po calej siatce  
		{ 			
			sum = 0;			

			for(k=0; k<4; k++) // liczenie efektu z jednego oczka po jego krawedziach
			{	  
				if(i<NumNet)
				{
					switch(k)
					{ 
					   case 0:  // krawedz lewa
						  x1 = p_net[i].x1 - x0;      
						  z1 = p_net[i].z2;    
						  x2 = p_net[i].x1 - x0;
						  z2 = p_net[i].z1;
						  break;
				  
					   case 1:  // krawedz dolna 
						  x1 = p_net[i].x2 - x0; 
						  x2 = p_net[i].x1 - x0;
						  z2 = p_net[i].z2;	  
						  break;
				
					   case 2:  // krawedz prawa
						  z1 = p_net[i].z1;
						  x2 = p_net[i].x2 - x0; 
						  break;

					   case 3:  //krawedz gorna
						  x1 = p_net[i].x1 - x0;
						  z1 = p_net[i].z1;
						  z2 = p_net[i].z1;
						  break;
					}
				}
				else // efekt z granic
				{
					if(i < (NumNet+2*NumNetZ))
						sig = 1; 

					switch(k)
					{ 
					   case 0:  // krawedz lewa
						  x1 = (p_net[i].x1 - x0)*sig;      
						  z1 = p_net[i].z2;    
						  x2 = (p_net[i].x1 - x0)*sig;
						  z2 = p_net[i].z1;
						  break;
				  
					   case 1:  // krawedz dolna 
						  x1 = (p_net[i].x2 - x0)*sig; 
						  x2 = (p_net[i].x1 - x0)*sig;
						  z2 = p_net[i].z2;	  
						  break;
				
					   case 2:  // krawedz prawa
						  z1 = p_net[i].z1;
						  x2 = (p_net[i].x2 - x0)*sig; 
						  break;

					   case 3:  // krawedz gorna
						  x1 = (p_net[i].x1 - x0)*sig;
						  z1 = p_net[i].z1;
						  z2 = p_net[i].z1;
						  break;
					}
				}				
			
				if(x1==0)
		    		x1 = 1e-12;
				if(x2==0)
					x2 = 1e-12;
					
				r1sq = pow(x1,2) + pow(z1,2);
				r2sq = pow(x1,2) + pow(z2,2);
				r3sq = pow(x2,2) + pow(z2,2);
				r4sq = pow(x2,2) + pow(z1,2);
				denom = z2 - z1;
				if(denom == 0)
					denom = 1e-12;
		
				r1 = sqrt(r1sq);
				r2 = sqrt(r2sq);
				r3 = sqrt(r3sq);
				r4 = sqrt(r4sq);
		
				theta1 = atan2(z1,x2);
				theta2 = atan2(z2,x2);
				theta3 = atan2(z2,x1);
				theta4 = atan2(z1,x1);
				
				term1 = sin(2*incc);
				term2 = log(r2*r3) - log(r4*r1);
				term3 = cos(2*incc);
				term4 = theta1 - theta2 - theta3 + theta4;
				
				sum = sum + (term1*term2 + (term3*term4));			
			}  							
			mg = mg + Fe*sum*p_net[i].suscept;	// efekt koncowy z oczka		
		}
						
		p_MG_results_temp[j][0] = x0;
		p_MG_results_temp[p_MG_sites-1-j][1] = mg; // wynik							
	}		
}

void get_input_files(char *model_file, char *g_file, char *mg_file)
{
	int i, j;
	i = 0;		     
	j = 0;

	FILE *plik;
	plik = fopen("input", "r");
	char line[120] = "";	
	char *res;

	res = fgets(line, sizeof line, plik); 
	if(res == NULL && !feof(plik))
	{
		puts(line);
		g_print("Blad odczytu linii\n");
	}
	
	for(i=(strstr(line, " ")-line)+1; i<strlen(line); i++) // nadpisz stringa od pozycji " "
	{		
		line[j] = line[i];										
		j++;
	}	
	strcat(model_file, line);
	model_file[strlen(model_file)-1] = 0;

	// wyczysc reszte stringa			
	for(j=j; j<strlen(line); j++)
	{
		line[j] = ' ';
	}	

	res = fgets(line, sizeof line, plik);  

	if(res == NULL && !feof(plik))
	{
		puts(line);
		g_print("Blad odczytu linii\n");
	}
	
	for(i=(strstr(line, " ")-line)+1; i<strlen(line); i++) // nadpisz stringa od pozycji " "
	{		
		line[j] = line[i];										
		j++;
	}	

	strcat(g_file, line);
	g_file[strlen(g_file)-1] = 0;

	// wyczysc reszte stringa			
	for(j=j; j<strlen(line); j++)
	{
		line[j] = ' ';
	}	

	res = fgets(line, sizeof line, plik);  

	if(res == NULL && !feof(plik))
	{
		puts(line);
		g_print("Blad odczytu linii\n");
	}
	
	for(i=(strstr(line, " ")-line)+1; i<strlen(line); i++) // nadpisz stringa od pozycji " "
	{		
		line[j] = line[i];										
		j++;
	}	

	strcat(mg_file, line);
	mg_file[strlen(mg_file)-1] = 0;

	// wyczysc reszte stringa			
	for(j=j; j<strlen(line); j++)
	{
		line[j] = ' ';
	}	

	fclose(plik);

	return;
}

double Magneto(model* calc_model, block *p_net, int p_MG_sites, double **p_MG_data, double **p_MG_results)
{			
	//FILE *results;
	
	int i = 0;	
	double misfit = 0.0;   		 // RMSE
	double *dwk4MG; 	 	     // tablica przechowujaca roznice miedzy spreparowanymi danymi polowymi a rzeczywistymi (dla pozniejszego RMSE) // KM
	double **MG_results_temp;
		
	struct timeval tic, toc;
	gettimeofday(&tic, NULL);
	//suscept = calc_model->susceptibility;
	dwk4MG = (double*)calloc(p_MG_sites, sizeof(double)); // tego nie wrzucamy do alokacji, bo ta tablica jest używana przy każdym wywołaniu funkcji		
	
	MG_results_temp = (double**)malloc((p_MG_sites)*sizeof(double*)); // alokacja pamieci na wyniki do MainForma
	for(i=0; i<p_MG_sites; i++) 
		MG_results_temp[i] = (double*)malloc(2*sizeof(double));

	calculate_mag(p_net, calc_model, MG_results_temp, p_MG_sites, p_MG_data); // fun. liczaca ef. magnetyczny 	
	
	if(areFieldMGData) // misfit ma sens tylko przy wczytaniu danych polowych
	{
		for(i=0; i<p_MG_sites; i++) 
		{                    
			dwk4MG[i] = p_MG_data[i][1] - MG_results_temp[i][1]; // kiedys musimy uwzglednic bledy pomiarowe, dzielac przez wektor bledu         
		} 
		
		misfit = anorm(p_MG_sites, dwk4MG); 
	}
	else 
		misfit = 1.0;

	p_MG_results = (double**)malloc((p_MG_sites)*sizeof(double*)); // alokacja pamieci na wyniki do MainForma
	for(i=0; i<p_MG_sites; i++) 
	{
		p_MG_results[i] = (double*)malloc(2*sizeof(double));
		p_MG_results[i][0] = MG_results_temp[i][0];
		p_MG_results[i][1] = MG_results_temp[i][1]; 
	}		
		
	for(i=0; i<p_MG_sites; i++)
		g_free(MG_results_temp[i]);
	g_free(MG_results_temp);
	
	g_free(dwk4MG);	
	gettimeofday(&toc, NULL);

	return misfit;	
}

int count_sites(char *p_file)
{
	int sites = 0;
	int res = 0;	
	FILE *data;
	data = fopen(p_file, "r");

	res = fscanf(data, "%d", &sites);

	if(res < 1) // nie wczytal ilosci stacji
		return -1;

	fclose(data);

	return sites;
}

void open_MG_data_auto(char *g_file, double **p_MG_data)
{
	double **MG_data_temp;
	int p_MG_sites = 0;

	FILE *data_MG;
	data_MG = fopen(g_file, "r");

	int i = 0;  	
	int j = 0;	
	int res = 0;									
		
	res = fscanf(data_MG, "%d", &p_MG_sites); 

	if(res < 1) // nie wczytal ilosci stacji
		return;
								
	MG_data_temp = (double**)malloc(p_MG_sites*sizeof(double*)); // alokacja pamieci na dane polowe
	for(j=0; j<p_MG_sites; j++)
		MG_data_temp[j]=(double*)malloc(2*sizeof(double));
			
	if(MG_data_temp==NULL) 
	{ 
		g_print("\n Przydzielenie pamięci dla danych polowych nie było możliwe"); 
		getchar(); 
	}	
	
	for(i=0; i<p_MG_sites; i++) // wczytanie danych polowych
	{
		res = fscanf(data_MG, "%lf %lf", &MG_data_temp[i][0], &MG_data_temp[i][1]); 

		if(res < 1) // nie wczytal danych dla stacji nr i
			return;
	}
		
	fclose(data_MG);
	
	for(i=0; i<p_MG_sites; i++)
	{
		p_MG_data[i][0] = MG_data_temp[i][0];
		p_MG_data[i][1] = MG_data_temp[i][1];
	}		
		
	for(i=0; i<p_MG_sites; i++)
        g_free(MG_data_temp[i]);
    g_free(MG_data_temp);
    MG_data_temp = NULL; 
				
	areFieldMGData = 1;	
	return;	
}

double tab_max(double **tab, int nelems)
{
	int i;
	double max = tab[0][1];

	for(i=0; i<nelems; i++)
	{
		if(tab[i][1] > max)
			max = tab[i][1];
	}

	return max;
}

double tab_min(double **tab, int nelems)
{
	int i;
	double min = tab[0][1];

	for(i=0; i<nelems; i++)
	{
		if(tab[i][1] < min)
			min = tab[i][1];	
	}

	return min;
}

block *read_model_auto(model *p_model1, char *model_file, block *p_net, GSList *p_ListOf_net_z_range)
{ 
	FILE *plik;
	plik = fopen(model_file, "r");

	char line[120];
	int i, j;
	int debugInfo = 0; // flaga determinujaca wyswietlanie w konsoli wartosci wczytanych z pliku	
	int LayerId;
	int BoundaryId;
	int VertexId;
	
	InfiNet = 2000;

	vertex *vertex1;
	boundary *boundary1;
	layer *layer1;	

	char *res;

	// SCZYTAJ Z PLIKU PARAMETRY MODELU	
	res = fgets(line, sizeof line, plik); 
	if(res == NULL && !feof(plik))
	{
		puts(line);
		g_print("Blad odczytu linii\n");
	}
	 	
	j = 0;
	for(i=(strstr(line, " ")-line)+1; i<strlen(line); i++) // nadpisz stringa od pozycji " "
	{		
		line[j] = line[i];										
		j++;
	}														
	for(i=0; i<strlen(line); i++)
		p_model1->name[i] = line[i];	
	
	read_line(line, plik);
	p_model1->id = atoi(line);
	read_line(line, plik);
	p_model1->xRange = atoi(line);
	read_line(line, plik);
	p_model1->zRange = atoi(line);
	read_line(line, plik);
	p_model1->step = atoi(line);
	read_line(line, plik);
	p_model1->nodes = atoi(line);	
	read_line(line, plik);
	p_model1->NumNetX = atoi(line);	
	read_line(line, plik);
	p_model1->NumNetZ = atoi(line);	
	read_line(line, plik);
	p_model1->thick1Lay = atoi(line);	
	read_line(line, plik);
	p_model1->widMidNod = atoi(line);
	read_line(line, plik);
	p_model1->regular = atoi(line);	
	read_line(line, plik);
	p_model1->density = atof(line);	
	read_line(line, plik);
	p_model1->resistivity = atof(line);	
	read_line(line, plik);
	p_model1->intensity = atof(line);
	read_line(line, plik);	
	p_model1->declination = atof(line);
	read_line(line, plik);
	p_model1->inclination = atof(line);
	read_line(line, plik);
	p_model1->angle_mag = atof(line);
	read_line(line, plik);
	p_model1->susceptibility = atof(line);
	read_line(line, plik);
	p_model1->ratioX = atof(line);
	read_line(line, plik);
	p_model1->ratioZ = atof(line);
	read_line(line, plik);
	// (p_model1->fillColour)->r = atof(line);
	read_line(line, plik);
	// (p_model1->fillColour)->g = atof(line);
	read_line(line, plik);	
	//(p_model1->fillColour)->b = atof(line);
	read_line(line, plik);
	//nLayers = atoi(line);	

	while(1) // sczytuj linie az do skonczenia pliku
	{
		res = fgets(line, sizeof line, plik); 
		if(res == NULL && !feof(plik))
		{
			puts(line);
			g_print("Blad odczytu linii\n");
		}
					
		if(feof(plik) != 0)	
			break;
		
		// SCZYTAJ Z PLIKU PARAMETRY WSZYSTKICH WARSTW										
		if((strstr(line, "Layer") - line) > 0) // jesli wszedl do sekcji nowej warstwy (zaczyna sie od '***')
		{							
			// zaalokuj warstwe i sczytaj reszte parametrow
			layer1 = NULL;
			layer1 = malloc(sizeof(layer)); 			
			layer1->listOfBoundaries = NULL;	
			layer1->listOfBoundariesId = NULL;			
							
			read_line(line, plik);							
			layer1->id = atoi(line);		
			read_line(line, plik);
			layer1->meterialId = atoi(line);			
			read_line(line, plik);
			layer1->used = atoi(line);					
			read_line(line, plik);		
			layer1->density = atof(line);					
			read_line(line, plik);					
			layer1->dens_range[0] = atof(line);					
			read_line(line, plik);		
			layer1->dens_range[1] = atof(line);					
			read_line(line, plik);	
			layer1->resistivity = atof(line);					
			read_line(line, plik);
			layer1->res_range[0] = atof(line);					
			read_line(line, plik);	
			layer1->res_range[1] = atof(line);					
			read_line(line, plik);
			layer1->intensity = atof(line);					
			read_line(line, plik);
			layer1->declination = atof(line);					
			read_line(line, plik);
			layer1->inclination = atof(line);					
			read_line(line, plik);
			layer1->angle_mag= atof(line);					
			read_line(line, plik);								
			layer1->susceptibility = atof(line);					
			read_line(line, plik);	
			layer1->susceptibility_range[0] = atof(line);					
			read_line(line, plik);
			layer1->susceptibility_range[1] = atof(line);					
			read_line(line, plik);															
			// (layer1->fillColour)->r = atof(line);		      
			read_line(line, plik);			
			// (layer1->fillColour)->g = atof(line);				       
			read_line(line, plik);
			// (layer1->fillColour)->b = atof(line);									

			p_model1->listOfLayersId =  g_slist_append(p_model1->listOfLayersId, GINT_TO_POINTER(layer1->id)); // dodaj warstwe do modelu	po id		
			p_model1->listOfLayers   =  g_slist_append(p_model1->listOfLayers, layer1);		 // dodaj warstwe do modelu		

			// g_free(layer1->fillColour);
			//g_free(layer1);
		}		
			
		// SCZYTAJ Z PLIKU PARAMETRY WSZYSTKICH GRANIC											
		if((strstr(line, "boundary") - line) > 0) // jesli wszedl do sekcji nowej granicy (zaczyna sie od '**')
		{						
			// nr granicy
			j = 0;
			for(i=(strstr(line, ".")-line)+1; i<strlen(line); i++) // nadpisz stringa od pozycji "." +1
			{		
				line[j] = line[i];										
				j++;
			}		
			
			// wyczysc reszte stringa			
			for(j=j; j<strlen(line); j++)
				line[j] = ' ';						
			
			// zaalokuj granice i sczytaj jej parametry
			boundary1 = NULL;
			boundary1 = malloc(sizeof(boundary)); 
			boundary1->listOfVertex = NULL;	
			boundary1->listOfVerticesId = NULL;			
			
			read_line(line, plik);
			boundary1->id = atoi(line);
			read_line(line, plik);					
			boundary1->fixed = atoi(line);
			read_line(line, plik);
			boundary1->used = atoi(line);							
			p_model1->listOfBoundaries   = g_slist_append(p_model1->listOfBoundaries, boundary1);	   // dodaj granice do modelu
		}
		
		// SCZYTAJ Z PLIKU PARAMETRY WSZYSTKICH WIERZCHOLKOW					
		if((strstr(line, "vertex") - line) > 0) // jesli wszedl do sekcji nowego wierzcholka (zaczyna sie od '*')
		{						
			// nr wierzcholka
			j = 0;
			for(i=(strstr(line, ".")-line)+1; i<strlen(line); i++) // nadpisz stringa od pozycji "." +1
			{		
				line[j] = line[i];										
				j++;
			}		
			
			// wyczysc reszte stringa			
			for(j=j; j<strlen(line); j++)
				line[j] = ' ';							
							
			// zaalokuj wierzcholek i sczytaj jego parametry	
			vertex1 = NULL;
			vertex1 = malloc(sizeof(vertex));	
			
			read_line(line, plik);
			vertex1->id = atoi(line);					
			read_line(line, plik);
			vertex1->x = atof(line);					
			read_line(line, plik);
			vertex1->z = atof(line);					
			read_line(line, plik);
			vertex1->fixed_x = atoi(line);					
			read_line(line, plik);
			vertex1->fixed_z = atoi(line);					
			read_line(line, plik);	
			vertex1->used = atoi(line);					
			read_line(line, plik);
			vertex1->xrange[0] = atof(line);					
			read_line(line, plik);	
			vertex1->xrange[1] = atof(line);					
			read_line(line, plik);	
			vertex1->zrange[0] = atof(line);					
			read_line(line, plik);	
			vertex1->zrange[1] = atof(line);											

			p_model1->listOfVertex = g_slist_append(p_model1->listOfVertex, vertex1); // dodaj wierzcholek do modelu	
		}

		layer1 = NULL;
		boundary1 = NULL;

		// SZCZYTAJ IDKI I UZUPELNIJ NIMI LISTY								
		if((strstr(line, "HERE") - line) > 0) // jesli wszedl do sekcji IDkow
		{
			if(debugInfo) puts(line);
			
			while(1)
			{
				res = fgets(line, sizeof line, plik); 
				if(res == NULL && !feof(plik))
				{
					puts(line);
					g_print("Blad odczytu linii\n");
				}
																	
				if( (feof(plik) != 0) || (strstr(line,"Mesh thicken parameters:")) )	
					break;						

				if((strstr(line, "Layer") - line) > 0)
				{
					read_line_id(line);																						
					LayerId = atoi(line);
					layer1 =((layer*)(g_list_id_get_data(p_model1,LayerId,3)));				
				}					
					
				if((strstr(line, "boundary") - line) > 0) 
				{
					read_line_id(line);																	
					BoundaryId = atoi(line);					  					
					layer1->listOfBoundariesId = g_slist_append(layer1->listOfBoundariesId, GINT_TO_POINTER (BoundaryId)); // dodaj granice do warstwy po id	
										
					// jesli istniala poprzednia granica do ktorej byly przypisane wierzcholki, zaznacz jako uzyta
					if(boundary1 != NULL)
					{
						if(g_slist_length(boundary1->listOfVerticesId) != 0)
							boundary1->used = 1;
					}

					boundary1 =((boundary*)(g_list_id_get_data(p_model1, BoundaryId,2)));
				}
				
				if((strstr(line, "vertex") - line) > 0)
				{
					read_line_id(line);
					VertexId = atoi(line);					
																							
					if(boundary1->used == 0)
						boundary1->listOfVerticesId = g_slist_append(boundary1->listOfVerticesId, GINT_TO_POINTER (VertexId)); // dodaj wierzcholek do granicy po id
				}			
			}
		}	
			
		if ( strstr(line,"Mesh thicken parameters:") )
		{
			net_z_intervals *temp;	
			temp = NULL;
			temp = malloc(sizeof(net_z_intervals));		
			temp ->z_step = 0; 
			temp ->depth  = 0; 
			
			GSList *mesh_params;
			mesh_params = NULL;
			mesh_params = g_slist_append(mesh_params,temp);
			
				
			while (1)
			{
				temp = NULL;
		
				res = fgets(line, sizeof line, plik); 
				if(res == NULL && !feof(plik))
				{
					puts(line);
					g_print("Blad odczytu linii\n");
				}
													   
				if(feof(plik) != 0) 
					break;										
				temp ->z_step  = atoi(line);
				
				res = fgets(line, sizeof line, plik); 
				if(res == NULL && !feof(plik))
				{
					puts(line);
					g_print("Blad odczytu linii\n");
				}
													   
				if(feof(plik) != 0) 
				{
					g_print("missing depth in mesh thicken parameters"); 
					break;	
				} 				
				else
				{
					temp ->depth  = atoi(line);
					mesh_params   = g_slist_append(mesh_params,temp);
				}
			}
			g_free(temp);
			thicken_mesh(p_model1, mesh_params, p_net, p_ListOf_net_z_range);	
		}
	}	
			
	// wyczyszczenie flag uzycia granic			
	for(i=0; i<g_slist_length(p_model1->listOfLayers); i++) // for layers of model
	{  		
		LayerId = (intptr_t)(g_slist_nth(p_model1->listOfLayersId, i)->data);   
		
		for(j=0; j<g_slist_length( ((layer*)(g_list_id_get_data(p_model1,LayerId,3))) -> listOfBoundariesId); j++) // for all boundaries
		{
			BoundaryId = (intptr_t)(g_slist_nth( (((layer*)(g_list_id_get_data(p_model1,LayerId,3)))->listOfBoundariesId ), j)->data);
			boundary1 =((boundary*)(g_list_id_get_data(p_model1, BoundaryId,2)));		
			boundary1->used = 0;				
		}
	}			
	isRegular = p_model1->regular;				
	p_net = net_maker(p_model1, isRegular, p_net, p_ListOf_net_z_range);	
	areFwdSolvers = 0;						
	fclose(plik);								

	return p_net;	
}

void run_fwd_solvers(model *p_model1, block *p_net, int p_MG_sites, double *p_MGmisfit, double **p_MG_data, double **p_MG_results)
{		
	if(fill_mesh_with_model(p_model1, p_net) == 0) // fill_mesh_with_model zwraca 0 gdy jest ok
	{			
		if(MGIsIn) *p_MGmisfit = Magneto(p_model1, p_net, p_MG_sites, p_MG_data, p_MG_results);	
        areFwdSolvers = 1;				
	}
	else
	{
		if(MGIsIn) *p_MGmisfit = DBL_MAX; // 1.0;	
        areFwdSolvers = 0; //?
	}
    if(*p_MGmisfit < DBL_MAX)
	{
		g_print("MG Misfit = %lf ", *p_MGmisfit);

		// model params	
		int count = 0;
    	int VertexId = 0;
    	vertex *vertextemp;	

		// od 5go do 10go wierzcholka
		for(VertexId = 5; VertexId < 11; VertexId++)
		{
			vertextemp = NULL;        
			vertextemp = (vertex*)(g_list_id_get_data(p_model1, VertexId,1));

			count += 2;

	        g_print("%lf, %lf, ", vertextemp->x, vertextemp->z);
		}
			
		layer *layer1, *layer2;
		int LayerId = 0;

		LayerId = (intptr_t)(g_slist_nth_data(p_model1->listOfLayersId,0)); // dla 1 warstwy
		layer1 = (layer*)(g_list_id_get_data(p_model1,LayerId,3));
		g_print("%lf, %lf, ",  layer1->density, layer1->susceptibility);

		LayerId = (intptr_t)(g_slist_nth_data(p_model1->listOfLayersId,1)); // dla 2 warstwy
		layer2 = (layer*)(g_list_id_get_data(p_model1,LayerId,3));
		g_print("%lf, %lf\n",  layer2->density, layer2->susceptibility);
	}	

	return;			
}  

void run_rasterization(model *p_model1, block *p_net, GSList *p_ListOf_net_z_range)
{	
	if(p_net == NULL)
		p_net = net_maker(p_model1, isRegular, p_net, p_ListOf_net_z_range);  // Create mesh 	
	
	fill_mesh_with_model(p_model1, p_net); 									           

	return;	   
}  

void thicken_mesh(model *p_model1, GSList *p_entThMeshParams, block *p_net, GSList *p_ListOf_net_z_range)
{
	if(p_entThMeshParams != NULL)
	{
		net_z_intervals *temp;
		temp = NULL;
		temp = malloc(sizeof(net_z_intervals));
		
		int i = 0;
		
		p_ListOf_net_z_range = NULL;
				
		for(i=1; i < g_slist_length(p_entThMeshParams); i++)
		{	
			temp = NULL;			
			temp -> z_step =  ((net_z_intervals*)(g_slist_nth_data(p_entThMeshParams,i))) -> z_step;
			temp -> depth  =  ((net_z_intervals*)(g_slist_nth_data(p_entThMeshParams,i))) -> depth ; 					
			p_ListOf_net_z_range = g_slist_append (p_ListOf_net_z_range, temp);
		}
		g_free(temp);
	}	
	
	g_free(p_net); 
	p_net = net_maker(p_model1, isRegular, p_net, p_ListOf_net_z_range);

	return;
}

void Magneto_new(double *p_model_params, double *misfit)
{   
	// = c(30.0, 198.0, 30.0, 110.0, 30.0, 20.0, 110.0, 20.0, 110.0, 110.0, 110.0, 198.0, 2100.0, 0.00005, 3000.0, 0.05)

	GSList *ListOf_net_z_range = NULL;

	int i = 0;
	int MG_sites = 0;  				// il. pkt. pomiarowych gravi
	char model_file[120] = "";
	char g_file[120] = "";
	char mg_file[120] = "";	

	double MGmisfit	= DBL_MAX; 
	double **MG_data; 	  // site localization [i][0] and field data [i][1] for Gravi - field data ok
	double **MG_results;   // site localization [i][0] and model response [i][1] for Gravi

	block *net1;
	net1 = NULL;
	
	isShiftCalculated = 0;
	areFieldMGData    = 0; // Flag: are Gravi field data read?
	areFwdSolvers     = 0;	// Flag: are the fwd solvers done?
	isShiftCalculated = 0;
	isRegular         = 1;	// Flag: is mesh regular? 0-no 1-yes

	model *model1;
	model1 = NULL;
	model1 = malloc(sizeof(model)); 
	//model1->fillColour = (colour*)(malloc(sizeof(colour)));	
	model1->listOfPolygons 	 = NULL; 
	model1->listOfLayers 	 = NULL; 
	model1->listOfLayersId   = NULL;
	model1->listOfBoundaries = NULL; 
	model1->listOfVertex 	 = NULL; 

	get_input_files(model_file, g_file, mg_file);	
	net1 = read_model_auto(model1, model_file, net1, ListOf_net_z_range);	
	MG_sites = count_sites(mg_file);
	alloc_2D_array(&MG_data, MG_sites, 2);
	alloc_2D_array(&MG_results, MG_sites, 2);

	open_MG_data_auto(mg_file, MG_data);
	
	// overwrite model params from given in R vector model_params
	// e.g. for prostokat_19II  x = c(30.0, 198.0, 30.0, 110.0, 110.0, 110.0, 110.0, 198.0, 30.0, 20.0, 110.0, 20.0, 2100.0, 0.00005, 3000.0, 0.05)
	
	// wierzcholki warstwy 2
	int count = 0;
	int VertexId = 0;
	vertex *vertextemp;	

	// od 5go do 10go wierzcholka
	for(VertexId = 5; VertexId < 11; VertexId++)
	{
		vertextemp = NULL;        
		//vertextemp = (vertex*)(g_slist_nth_data(model1->listOfVertex, VertexId));    
		vertextemp = (vertex*)(g_list_id_get_data(model1,VertexId,1));
		vertextemp->x = p_model_params[count];  
		vertextemp->z = p_model_params[count+1];
		count += 2;
	}

	model1->density        = p_model_params[12]; 
	model1->susceptibility = p_model_params[13]; 
	
	layer *layer_a, *layer_b;

	int LayerId = 0;

	LayerId = (intptr_t)(g_slist_nth_data(model1->listOfLayersId,0)); // dla 1 warstwy
	layer_a = (layer*)(g_list_id_get_data(model1,LayerId,3));
	layer_a->density = p_model_params[12];
	layer_a->susceptibility = p_model_params[13];

	LayerId = (intptr_t)(g_slist_nth_data(model1->listOfLayersId,1)); // dla 2 warstwy
	layer_b = (layer*)(g_list_id_get_data(model1,LayerId,3));
	layer_b->density = p_model_params[14];
	layer_b->susceptibility = p_model_params[15];
	
	run_rasterization(model1, net1, ListOf_net_z_range);
	run_fwd_solvers(model1, net1, MG_sites, &MGmisfit, MG_data, MG_results);

	*misfit = MGmisfit;

	if(areFieldMGData)
	{
		for(i=0; i<MG_sites; i++)
			g_free(MG_data[i]);
		g_free(MG_data);
	}

	if(areFwdSolvers)
	{
		for(i=0; i<MG_sites; i++)
			g_free(MG_results[i]);
		g_free(MG_results);      
	} 

	g_free(net1);

	for(i = 1; i <= g_slist_length(model1->listOfVertex); i++)			
		g_free((vertex*)(g_list_id_get_data(model1,i,1)));

	for(i = 1; i <= g_slist_length(model1->listOfBoundaries); i++)			
		g_free((boundary*)(g_list_id_get_data(model1,i,2)));

	for(i = 1; i <= g_slist_length(model1->listOfLayers); i++)			
		g_free((layer*)(g_list_id_get_data(model1,i,3)));

	g_free(model1);

	muntrace();

    return;
}