#ifndef UTILS_H_
#define UTILS_H_

#include <math.h>
#include <stdlib.h>
#include <vector>

#define euler 	2.71828182846
#define SIGMA 	1.0
#define MI 	  	0.0

//Operações entre espaço de cores
std::vector<float> getLabColor(unsigned int sR, unsigned int sG,
		unsigned int sB);
std::vector<float> getRGBColor(float L, float a, float b);

//Operações com vetor
float sqED_2D(float x1, float y1, float x2, float y2);

#endif

