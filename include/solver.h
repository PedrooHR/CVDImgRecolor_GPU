/*
 * solver.h
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <vector>
#include <iostream>

#include "grid.h"
#include "dataset.h"
#include "utils.h"

#define T_GRID 100
#define GPU_ENABLE 1

#define MAX_VALUE 1e10

class Solver {
public:
	Solver();
	Solver(Grid *A, Dataset *B);

	Grid * grid;
	Dataset * dataset;
	int numNodes;
	bool inverte;
	int numEpochs;
	float MSE = MAX_VALUE;
	int GridDimension;
	int CurrentEpoch;
	float Rmodule;
	float Emodule;
	int zeroNode;

	//Projeção final
	std::vector<std::vector<float> > OriginalMap;
	std::vector<float> zMax;
	std::vector<float> zMin;

	std::vector<std::vector<int> > Taxons;
	std::vector<float> EpochM;
	std::vector<std::vector<float> > SysMatrix;
	std::vector<std::vector<float> > ResultVector;
	std::vector<std::vector<float> > Solution;
	std::vector<std::vector<float> > NodesSolution;
	void drawRecolored(const char * filepath);

	void resetSysMatrix();
	void constructSysMatrix();
	void calcTaxons();
	void solveLS();
	void ElMap();
	void centerWhite();
	void invert();
	void projectPoints();
};

#endif /* SOLVER_H_ */
