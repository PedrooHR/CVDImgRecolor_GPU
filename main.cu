#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "grid.h"
#include "solver.h"
#include "dataset.h"
#include "utils.h"

#define INVERTER 1
#define CENTERWHITE 1

// Programa Principal
int main(int argc, char **argv) {
	char * imgpath = "data/testeimg1.jpg";
	/* Read input from the arguments */
	if (argc > 1)
		imgpath = argv[1];

	/* Init data structures */
	Grid *grid = new Grid(PROTANOPE, T_GRID);
	Dataset *dataset = new Dataset(imgpath);
	Solver *solver = new Solver(grid, dataset);
	printf("Data Structures Initialized\n");

	/* Run algorithm */
	int max_epochs = solver->numEpochs;
	while (solver->CurrentEpoch < max_epochs) {
		solver->constructSysMatrix();
		solver->solveLS();
		solver->CurrentEpoch++;
		printf("Done Epoch %d\n", solver->CurrentEpoch);
	}
	solver->centerWhite();
	printf("Elastic map white centering\n");
	solver->invert();
	printf("Inverting process done\n");
	solver->projectPoints();
	printf("Points projected on the CVD plane\n");
	solver->drawRecolored("output.jpg");
	printf("Image Saved as output.jpg\n");
	exit(0);

}

