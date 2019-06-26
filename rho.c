#include "rho.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <math.h>

void write(const char * path)
{
	struct header header = {
		.x_min = -1., .x_max = 1., .nx = 3,
		.y_min = -2., .y_max = 2., .ny = 5,
		.z_min = -3., .z_max = 3., .nz = 7 };
	float * data = NULL;
	FILE * fid;

	size_t data_size = header.nx * header.ny * header.nz * sizeof(*data);
	data = malloc(data_size);
	// init data with 0
	memset(data, 0., data_size);
	
	// define x, y, z pixel on each axis
	float x_pixel, y_pixel, z_pixel;
	x_pixel = (header.x_max - header.x_min) / header.nx;
	y_pixel = (header.y_max - header.y_min) / header.ny;
	z_pixel = (header.z_max - header.z_min) / header.nz;
	float x, y, z;
	
	for(int i = 0; i < header.nx; i++) {
		x = header.x_min + x_pixel * (i + 0.5);
		for(int j = 0; j < header.ny; j++) {
			y = header.y_min + y_pixel * (j + 0.5);
			for(int k = 0; k < header.nz; k++) {
				z = header.z_min + z_pixel * (k + 0.5);
				// get rho(x, y, z)
				data[header.nz * (i * header.ny + j) + k] = rho_model(x, y, z);
			}
		}
	}
	
	fid = fopen(path, "wb");
	fwrite(&header, sizeof(header), 1, fid);
	fwrite(data, sizeof(*data), header.nx * header.ny * header.nz, fid);
	
	fclose(fid);
	
	free(data);
}


struct model * read(const char * path)
{
	struct header header;
	struct model * model;
	FILE * fid;

	fid = fopen(path, "rb");

	fread(&header, sizeof(header), 1, fid);	
	model = malloc(sizeof(*model) + header.nx * header.ny * header.nz * sizeof (*model->data));
	
	memcpy(&model->header, &header, sizeof(header));
	fread(model->data, sizeof (*model->data), header.nx * header.ny * header.nz, fid);
	
	fclose(fid);

	return model;
}

void rho_display(struct model * model) {
	const struct header * const header = &model->header;
	float x_pixel, y_pixel, z_pixel;
	float x, y, z;
	int m, n, l;

	x_pixel = (header->x_max - header->x_min) / header->nx;
	y_pixel = (header->y_max - header->y_min) / header->ny;
	z_pixel = (header->z_max - header->z_min) / header->nz;

	for(int i = 0; i < header->nx; i++) {
		x = header->x_min + x_pixel * (i + 0.5);
		m = (int)((x - header->x_min) / x_pixel - 0.5); // must be in [0, header->nx - 1]
	
		for(int j = 0; j < header->ny; j++) {
			y = header->y_min + y_pixel * (j + 0.5);
			n = (int)((y - header->y_min) / y_pixel - 0.5); // must be in [0, header->ny - 1]
	
			for(int k = 0; k < header->nz; k++) {
				z = header->z_min + z_pixel * (k + 0.5);
				l = (int)((z - header->z_min) / z_pixel - 0.5); // must be in [0, header->nz - 1]
				printf("[%d, %d, %d] -> rho(%.10f, %.10f, %.10f) = %.10f\n", m, n, l, x, y, z, model->data[header->nz * (m * header->ny + n) + l]);
			}
		}
	}
}




struct rho_xyz_cell * rho_xyz_get(struct model * model, 
	float x, float y, float z, 
	float ux, float uy, float uz) {

	const struct header * const header = &model->header; // alias from model->header to header
	int i, j, k;
	float x_pixel, y_pixel, z_pixel;

	x_pixel = (header->x_max - header->x_min) / header->nx;
	y_pixel = (header->y_max - header->y_min) / header->ny;
	z_pixel = (header->z_max - header->z_min) / header->nz;

	i = (int)((x - header->x_min) / x_pixel - 0.5); // must be in [0, header->nx - 1]
	j = (int)((y - header->y_min) / y_pixel - 0.5); // must be in [0, header->ny - 1]
	k = (int)((z - header->z_min) / z_pixel - 0.5); // must be in [0, header->nz - 1]

	struct rho_xyz_cell * rho_xyz_cell;
	rho_xyz_cell = NULL;
	rho_xyz_cell = malloc(sizeof(*rho_xyz_cell));

	// save indices (i, j, k)
	rho_xyz_cell->i = i;
	rho_xyz_cell->j = j;
	rho_xyz_cell->k = k;

	// in the model
	if(i >= 0 && i < header->nx && j >= 0 && j < header->ny && k >= 0 && k < header->nz) {
		rho_xyz_cell->rho = model->data[header->nz * (i * header->ny + j) + k];

		// get the distance with the next nearest cell taking into account the entire model
		rho_xyz_cell->distance_next_cell = -1.; // for avoiding errors

		// get the center of that voxel
		float voxel_cx, voxel_cy, voxel_cz;
		voxel_cx = header->x_min + x_pixel * (i + 0.5);
		voxel_cy = header->y_min + y_pixel * (j + 0.5);
		voxel_cz = header->z_min + z_pixel * (k + 0.5);
		// get min distance
		rho_xyz_cell->distance_next_cell = distance_next_voxel(x, y, z, voxel_cx, 
			voxel_cy, voxel_cz, x_pixel, y_pixel, z_pixel, ux, uy, uz);
		
	}
	else {
		rho_xyz_cell->rho = -1.;
		// get the distance with the next nearest cell taking into account the entire model
		rho_xyz_cell->distance_next_cell = -1.; // for avoiding errors
	}
	
	memcpy(&rho_xyz_cell, &rho_xyz_cell, sizeof(*rho_xyz_cell));
	return rho_xyz_cell;
}

float distance_next_voxel(float x, float y, float z, // position
	float x0, float y0, float z0,                    // center of the voxel
	float sizex, float sizey, float sizez, 	         // size of the voxel
	float ux, float uy, float uz) {                  // direction
	
	float distance = INFINITY;

	float X, Y, Z;
	float t1, t2, t3;
	t1 = t2 = t3 = INFINITY;

	// plans in zoz directions
	if(uz != 0) {
		t3 = (z0 + 0.5 * sizez - z) / uz;
		if(t3 < 0) {
			t3 = (z0 - 0.5 * sizez - z) / uz;
		}
	}
	
	if(t3 < 0) {
		t3 = INFINITY;
	}

	if(distance > t3)
		distance = t3;

	// plans in yoy directions
	if(uy != 0) {
		t2 = (y0 + 0.5 * sizey - y) / uy;
		if(t2 < 0) {
			t2 = (y0 - 0.5 * sizey - y) / uy;
		}
	}

	if(t2 < 0) {
		t2 = INFINITY;
	}

	if(distance > t2)
		distance = t2;

	// plans in xox directions
	if(ux != 0) {
		t1 = (x0 + 0.5 * sizex - x) / ux;
		if(t1 < 0) {
			t1 = (x0 - 0.5 * sizex - x) / ux;
		}
	}

	if(t1 < 0) {
		t1 = INFINITY;
	}

	if(distance > t1)
		distance = t1;
	
	distance *= sqrt(ux * ux + uy * uy + uz * uz);
	
	return distance;
}


float rho_model(float x, float y, float z) {
	float rho_xyz;
	// define rho(x, y, z) with a special known model
	rho_xyz = 2 + (x + y + z) / (x * x + y * y + z * z + x + y + z); // 0.;

	if(rho_xyz < 0)
		rho_xyz *= -1;

	return rho_xyz;
}
