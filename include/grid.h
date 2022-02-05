/*
 * grid.h
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#ifndef GRID_H_
#define GRID_H_

#include <vector>

#define PROTANOPE 	1
#define DEUTERANOPE 0
#define TRITANOPE   0

typedef struct Node {
	int id;
	float Weight;
	float *CVDposition;
	float *Position;
	float a;
	float b;
} Node;

class Grid {
public:
	// LIMITES DOS PLANOS
	std::vector<float> A;
	std::vector<float> C;

	int GRAPH_SIZE;
	int numRibs;
	int numEdges;

	float miAB;
	float miBC;

	std::vector<std::vector<int> > Edges;
	std::vector<std::vector<int> > Ribs;
	Node *Graph;

	Grid(int type, int graph_size);
};

#endif /* GRID_H_ */
