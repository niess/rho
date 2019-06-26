#include "rho.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


int main(int argc, char * argv[])
{
	const char * path = (argv[1] != NULL) ? argv[1] : "geometry.b";

	write(path);
	struct model * model = read(path);

	// display to verify
	rho_display(model);

	float x, y, z, ux, uy, uz;	

	x = 0.6666667461;
	y = 0.8000000715;
	z = -1.7142857313;
	ux = 0.;
	uy = 0.;
	uz = 1.;

	struct rho_xyz_cell * rho_xyz_cell = rho_xyz_get(model, x, y, z, ux, uy, uz);

	printf("============================================================================\n");
	printf("rho(%.10f, %.10f, %.10f) = %.10f\n", x, y, z, rho_xyz_cell->rho);
	printf("----------------------------------------------------------------------------\n");
	printf("distance with the next nearest cell for u(%.1f, %.1f, %.1f) : %.10f\n", ux, uy, uz, rho_xyz_cell->distance_next_cell);
	printf("============================================================================\n");



	free(rho_xyz_cell);

	free(model);

	exit(EXIT_SUCCESS);
}
