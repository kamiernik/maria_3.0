typedef struct {
	
    double r;
    double g;
    double b;
    
} colour;

typedef struct {
	 
	int x1;  
	int z1; 
	int x2;
	int z2;
	double rhob;
	double res;
	double inc; // inclination
	double angm; // angle of magnetisation
	double ints; // magnetisation
	double suscept; // susceptibility
	colour *fillColour;
	int flag; 
	
} block; 

typedef struct {

	int x1;
	int z1;
	int x2;
	int z2;

	double rhob; // not necessary

	int ifBelongs; // 1 - if point belongs to dajka, 0 if not
	
} modelmap;

typedef struct {
	
	int id;
	double x;
	double z;	
	int fixed_x; // 0-editable, 1-blocked
	int fixed_z; // 0-editable, 1-blocked
 	int used;
	double xrange[2];
	double zrange[2];
 	
} vertex;

typedef struct {

	int id;
	vertex v1;
	vertex v2;

} segment;

typedef struct {
	
	GSList *listOfVertex;
	GSList *listOfVerticesId;
	double density;
	double resistivity;	
	double intensity;
	double declination;
	double inclination;
	double magnetisation; 
	double angle_mag;
	double susceptibility; 
	colour *fillColour;
	
} polygon;

typedef struct{
	
	int id;
	GSList *listOfPolygons;	
	GSList *listOfLayers;	
	GSList *listOfLayersId;	
	GSList *listOfBoundaries;
	GSList *listOfVertex;
	int xRange;
	int zRange;
	int step; 
	int nodes;  // if regular mesh nodes = (NumNetX/step) * (NumNetZ/step)
	int NumNetX;
	int NumNetZ;
	int thick1Lay;	
	int widMidNod;	
	int regular;	
	double density;
	double resistivity;
	double intensity;
	double declination;
	double inclination;	
	double magnetisation; 
	double angle_mag;
	double susceptibility;
	double ratioX;		
	double ratioZ;			
	char name[50];
	colour *fillColour;

} model;

typedef struct {
	
	int id;
	GSList *listOfBoundaries;
	GSList *listOfBoundariesId;
	double density;
	double resistivity;
	double intensity;
	double declination;
	double inclination;
	double angle_mag; 
	double susceptibility; 	
	int meterialId;
	int used;
	double dens_range[2];
	double res_range[2];
	double intensity_range[2];
	double angmag_range[2];
	double susceptibility_range[2];
	colour *fillColour;

} layer;

typedef struct {
	
	int id;
	GSList *listOfVertex;
	GSList *listOfVerticesId;
	int fixed; // 0-editable, 1-blocked
	int used;
	
} boundary;

typedef struct {
	
	int id;
	char *name;
	double density;
	double resistivity;
	double inclination;
	double intensity;
	double angle_mag; 	
	double susceptibility;

} material;

typedef struct {
	
	int LayerId;
	int NearestVertexId[2];
	int VertexId;
		
} marginal_points;

typedef struct {
	
	int depth;
	int z_step;
		
} net_z_intervals;

typedef struct {
	
	int maxiter;    // maximum amount of iterations
	int nparticles; // no. of particles
	double eps; 	// accuracy - if more than 50% of particles stick fast arount local minimum (difference = eps) and have to be catapulted		
	int nthreads;	// no. of threads for parallel pso

} psoParamsS;

// structure holding pso Params and entries for them
typedef struct {
	
	psoParamsS *psoParams;
	GSList *entPSOParams;   // list of entries
	
} psoVarsS;
