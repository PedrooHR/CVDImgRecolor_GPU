/*
 * grid.cpp
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#include "grid.h"
#include <iostream>
#include "utils.h"

Grid::Grid(int type, int graph_size) {
	GRAPH_SIZE = graph_size;
	numEdges = GRAPH_SIZE - 1;
	numRibs = GRAPH_SIZE - 2;
	cudaMallocHost((void**) &Graph, GRAPH_SIZE * sizeof(Node));
	// do graphnodes
	if (type == PROTANOPE) {
		//PROTANOPE LIMITS
		A = {8.648425, -73.086372, 56.664734};
		C = {-14.907598, 86.293831, 89.536812};

		miAB = A[1] / A[0];
		miBC = C[1] / C[0];

		float x_step = fabs((A[0] - C[0]) / (GRAPH_SIZE * 1.0));
		float x_start = C[0];

		for(int i = 0; i < GRAPH_SIZE; i++) {
			Graph[i].id = i;
			Graph[i].Position = (float*) malloc (2*sizeof(float));
			Graph[i].Position[0] = (x_start <= 0) ? miBC * x_start : miAB * x_start;
			Graph[i].Position[1] = -x_start;
			Graph[i].CVDposition = (float*) malloc (2*sizeof(float));
			Graph[i].CVDposition[0] = x_start;
			Graph[i].CVDposition[1] = (x_start <= 0) ? miBC * x_start : miAB * x_start;
			Graph[i].Weight = 1.0f;
			Graph[i].a = (x_start <= 0) ? miBC * x_start : miAB * x_start;
			Graph[i].b = -x_start;
			x_start += x_step;
		}
	}

	// do edges
	for (int i = 0; i < GRAPH_SIZE - 1; i++) {
		std::vector<int> t = { i, i + 1 };
		Edges.push_back(t);
	}

	// do ribs
	for (int i = 0; i < GRAPH_SIZE - 2; i++) {
		std::vector<int> t = { i + 1, i, i + 2 };
		Ribs.push_back(t);
	}
}

