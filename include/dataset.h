/*
 * dataset.h
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#ifndef DATASET_H_
#define DATASET_H_

#define cimg_display 0
#define cimg_use_jpeg

#include <vector>
#include <CImg.h>

typedef struct Pixel {
	float* RGB;
	float* Lab;
	float R;
	float G;
	float B;
	float L;
	float a;
	float b;
	float PL;
	float Pa;
	float Pb;
	int id;
	int closest_node;
} Pixel;

class Dataset {
public:
	int width;
	int height;

	Pixel * Datapoints;
	int Datasize;

	Dataset(const char*imgPath);
};

#endif /* DATASET_H_ */
