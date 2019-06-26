#ifndef __RHO__
#define __RHO__

struct header {
	float x_min, x_max, y_min, y_max, z_min, z_max;
	int nx, ny, nz;
};

struct model {
	struct header header;
	float data[];
};

struct rho_xyz_cell{
	int i, j, k;               // cell
	float rho;                 // density
	float distance_next_cell;  // distance to the next nearest cell
};


void write(const char * path);
struct model * read(const char * path);
void rho_display(struct model * model);

/*
	gives rho(x, y, z)
*/
float rho_model(float x, float y, float z);

// return rho > 0 if in model & -1 if not
struct rho_xyz_cell * rho_xyz_get(struct model * model, 
	float x, float y, float z, 
	float ux, float uy, float uz
);

/*
	min distance from a point (x, y, z) to a voxel decribes by its center (x0, y0, z0) and size in a direction u
*/
float distance_next_voxel(float x, float y, float z, // position
	float x0, float y0, float z0,              // center of the voxel
	float sizex, float sizey, float sizez, 	   // size of the voxel
	float ux, float uy, float uz);             // direction

#endif
