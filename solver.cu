/*
 * solver.cpp
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#include <thrust/device_ptr.h>
#include "solver.h"
#include "grid.h"
#include "dataset.h"

#define Min2(x,y) x < y ? x : y
#define Min3(x,y,z) Min2(x,Min2(y,z))
#define Max2(x,y) x > y ? x : y
#define Max3(x,y,z) Max2(x,Max2(y,z))

Solver::Solver(Grid *A, Dataset *B){
	inverte = false;
	grid = A;
	dataset = B;
	numEpochs = 7;
	zeroNode = -1;
	Emodule = 0.1;
	Rmodule = 50;
	CurrentEpoch = 0;
	GridDimension = 2;
	numNodes = A->GRAPH_SIZE;
	SysMatrix.resize(numNodes);
	Taxons.resize(numNodes);
	ResultVector.resize(numNodes);
	Solution.resize(numNodes);
	NodesSolution.resize(dataset->Datasize);
	for (int i = 0; i < numNodes; i++){
		SysMatrix[i].resize(numNodes);
		Solution[i].resize(2);
		ResultVector[i].resize(2);
	}
	float epo = 0.05;
	for (int i = 0; i < numEpochs; i++){
		EpochM.push_back(epo);
		epo /= (3 - i*0.04);
	}
}

void Solver::resetSysMatrix(){
	// Clear stage
	SysMatrix.clear();
	Taxons.clear();
	ResultVector.clear();
	Solution.clear();

	//Resize stage
	SysMatrix.resize(numNodes);
	Taxons.resize(numNodes);
	ResultVector.resize(numNodes);
	Solution.resize(numNodes);
	for (int i = 0; i < numNodes; i++){
		SysMatrix[i].resize(numNodes);
		Solution[i].resize(2);
		ResultVector[i].resize(2);
	}
}

void Solver::constructSysMatrix(){
	resetSysMatrix();

	Emodule = EpochM[CurrentEpoch]*pow(grid->numEdges, (2 - GridDimension)/GridDimension)*10;
	Rmodule = EpochM[CurrentEpoch]*pow(grid->numRibs, (2 - GridDimension)/GridDimension)*100;

	calcTaxons();

 	// Funcional
	for (int i = 0; i < grid->numRibs; i++){
		int N0 = grid->Ribs[i][0];
		int N1 = grid->Ribs[i][1];
		int N2 = grid->Ribs[i][2];

		SysMatrix[N0][N0] += Rmodule * 4.0;
		SysMatrix[N0][N1] -= Rmodule * 2.0;
		SysMatrix[N0][N2] -= Rmodule * 2.0;
		SysMatrix[N1][N0] -= Rmodule * 2.0;
		SysMatrix[N2][N0] -= Rmodule * 2.0;
		SysMatrix[N1][N1] += Rmodule;
		SysMatrix[N1][N2] += Rmodule;
		SysMatrix[N2][N1] += Rmodule;
		SysMatrix[N2][N2] += Rmodule;
	}

	for (int i = 0; i < grid->numEdges; i++){
		int N0 = grid->Edges[i][0];
		int N1 = grid->Edges[i][1];

		SysMatrix[N0][N0] += Emodule;
		SysMatrix[N1][N1] += Emodule;
		SysMatrix[N0][N1] -= Emodule;
		SysMatrix[N1][N0] -= Emodule;
	}

	for (int i = 0; i < numNodes; i++){
		SysMatrix[i][i] += Taxons[i].size() / (1.0*dataset->Datasize);
	}

	for (int i = 0; i < numNodes; i++){
		ResultVector[i][0] = 0;
		ResultVector[i][1] = 0;
		for (int j = 0; j < Taxons[i].size(); j++){
			ResultVector[i][0] += dataset->Datapoints[Taxons[i][j]].Lab[1];
			ResultVector[i][1] += dataset->Datapoints[Taxons[i][j]].Lab[2];
		}
		ResultVector[i][0] /= (1.0*dataset->Datasize);
		ResultVector[i][1] /= (1.0*dataset->Datasize);
	}
}

void Solver::solveLS(){
	//escalonador
	for (int k = 0; k < numNodes; k++) {
		for (int j = k + 1; j < numNodes; j++){
			float multi = SysMatrix[j][k] / SysMatrix[k][k];
			for (int i = k; i < numNodes; i++)
				SysMatrix[j][i] -= (SysMatrix[k][i] * multi);

			ResultVector[j][0] -= ResultVector[k][0] * multi;
			ResultVector[j][1] -= ResultVector[k][1] * multi;
		}
	}
	//resolve
	for (int k = numNodes-1; k >= 0; k--) {
		double somax = 0;
		double somay = 0;
		for (int j = k + 1; j < numNodes; j++){
			somax += SysMatrix[k][j] * Solution[j][0];
			somay += SysMatrix[k][j] * Solution[j][1];
		}
		Solution[k][0] = (ResultVector[k][0] - somax) / SysMatrix[k][k];
		Solution[k][1] = (ResultVector[k][1] - somay) / SysMatrix[k][k];
	}

	for (int i = 0; i < numNodes; i++){
			grid->Graph[i].Position[0] = Solution[i][0];
			grid->Graph[i].Position[1] = Solution[i][1];
	}
}

__global__ void calcTaxonI (Pixel *Datapoints, Node *Graph, int size){
	int i = threadIdx.x + (blockDim.x * blockIdx.x);
	if (i < size) {
		float MinDist = (float) MAX_VALUE;
		int MinNodeRef;
		for (int j = 0; j < T_GRID; j++) {
			float localDist = (Datapoints[i].a - Graph[j].a)
					* (Datapoints[i].a - Graph[j].a);
			localDist += (Datapoints[i].b - Graph[j].b)
					* (Datapoints[i].b - Graph[j].b);

			if (localDist < MinDist) {
				MinDist = localDist;
				MinNodeRef = j;
			}
		}
		Datapoints[i].closest_node = MinNodeRef;
	}
}

void Solver::calcTaxons(){
	for (int i = 0; i < numNodes; i++)
		Taxons[i].clear();

	int div = (dataset->Datasize / 1024) + 1;
	calcTaxonI <<<div, 1024>>> (dataset->Datapoints, grid->Graph, dataset->Datasize);

	cudaDeviceSynchronize();

	for (int i = 0; i < dataset->Datasize; i++){
		Taxons[dataset->Datapoints[i].closest_node].push_back(i);
	}
}

void Solver::invert(){
	float sumDistAB = 0, sumDistBA = 0;
		for(int i = 0; i < Taxons.size(); i++){ //para cada nodo
			for (int j = 0; j < Taxons[i].size(); j++){

				float Deuto[3][3] = {{0.367322, 0.860646, -0.227968},
				                    {0.280085, 0.672501, 0.047413},
				                    {-0.011820, 0.042940, 0.968881}};
				int current_point = Taxons[i][j];

				std::vector<float> algproj = {grid->Graph[i].CVDposition[0], grid->Graph[i].CVDposition[1], dataset->Datapoints[current_point].Lab[0]};
				std::vector<float> alternativa = {grid->Graph[numNodes - i -1 ].CVDposition[0], grid->Graph[numNodes - i -1 ].CVDposition[1], dataset->Datapoints[current_point].Lab[0]};
				std::vector<float> original = {dataset->Datapoints[current_point].Lab[1], dataset->Datapoints[current_point].Lab[2], dataset->Datapoints[current_point].Lab[0]};

				std::vector <float> RgbProjection;
				RgbProjection.resize(3);
				RgbProjection[0] = std::min(std::max(Deuto[0][0]*original[0] + Deuto[0][1] * original[1] + Deuto[0][2]*original[2],0.0f),1.0f) * 255.0;
				RgbProjection[1] = std::min(std::max(Deuto[1][0]*original[0] + Deuto[1][1] * original[1] + Deuto[1][2]*original[2],0.0f),1.0f) * 255.0;
				RgbProjection[2] = std::min(std::max(Deuto[2][0]*original[0] + Deuto[2][1] * original[1] + Deuto[2][2]*original[2],0.0f),1.0f) * 255.0;

				sumDistAB += sqED_2D(algproj[0], algproj[1], RgbProjection[0], RgbProjection[1]);
				sumDistBA += sqED_2D(alternativa[0], alternativa[1], RgbProjection[0], RgbProjection[1]);
			}
		}
		inverte = false;

		if (sumDistAB > sumDistBA)
			inverte = true;
}

void Solver::centerWhite(){
	float min_distance = (float) MAX_VALUE;
	int zeroNode = numNodes / 2;
	for (int i = 0; i < numNodes; i++){
		for (int j = 0; j < Taxons[i].size(); j++){
			int point = Taxons[i][j];
			float localZeroDist = sqED_2D(dataset->Datapoints[point].a, dataset->Datapoints[point].b, 0, 0);
			if (localZeroDist < min_distance){
				zeroNode = i;
				min_distance = localZeroDist;
			}
		}
	}

	int size_positive = numNodes - zeroNode;
	int size_negative = zeroNode;
	float x_step = fabs((size_positive > size_negative) ? (grid->A[0] / (1.0 * size_positive)) : (-grid->C[0] / (1.0 * size_negative)));

	OriginalMap.resize(numNodes);

	float x_start = 0;
	for (int i = zeroNode; i < numNodes; i++){
		std::vector <float> t = {x_start, grid->miAB * x_start };
		OriginalMap[i] = t;
		x_start += x_step;
	}

	x_start = 0 - x_step;
	for (int i = zeroNode - 1; i >= 0; i--){
		std::vector <float> t = {x_start, grid->miBC * x_start};
		OriginalMap[i] = t;
		x_start -= x_step;
	}
}

void Solver::projectPoints(){
	NodesSolution.clear();
	NodesSolution.resize(dataset->Datasize);
	for(int i = 0; i < Taxons.size(); i++){ //para cada nodo
		for (int j = 0; j < Taxons[i].size(); j++){
			int current_point = Taxons[i][j];
			std::vector<float> algproj, original;

			algproj = {grid->Graph[i].CVDposition[0], grid->Graph[i].CVDposition[1], dataset->Datapoints[current_point].Lab[0]};

			if(OriginalMap.size() > 0){
				algproj = {OriginalMap[i][0], OriginalMap[i][1], dataset->Datapoints[current_point].Lab[0]};
			}

			NodesSolution[current_point].push_back(algproj[0]);
			NodesSolution[current_point].push_back((inverte) ? -algproj[1] : algproj[1]);
			NodesSolution[current_point].push_back(algproj[2]);
		}
	}
}

void Solver::drawRecolored(const char * filepath){
	cimg_library::CImg<float> image(dataset->width, dataset->height, 1, 3, 0);
	for(int i = 0; i < NodesSolution.size(); i++){
		int x = i%dataset->width;
		int y = i/dataset->width;
		std::vector <float> RGBcolor = getRGBColor(NodesSolution[i][2], NodesSolution[i][0], NodesSolution[i][1]);
		for (int j = 0; j < 3; j++){
			RGBcolor[j] = (RGBcolor[j] > 1.0) ? 1.0 : RGBcolor[j];
			RGBcolor[j] = (RGBcolor[j] < 0.0) ? 0.0 : RGBcolor[j];
		}
		image(x,y,0) = RGBcolor[0]*255.0;
		image(x,y,1) = RGBcolor[1]*255.0;
		image(x,y,2) = RGBcolor[2]*255.0;
	}
	image.save_jpeg(filepath, 100);
}
