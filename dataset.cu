/*
 * dataset.cpp
 *
 *  Created on: 17/03/2018
 *      Author: phr
 */

#include "dataset.h"
#include "utils.h"

Dataset::Dataset(const char *imgPath){
	int local_id = 0;
	cimg_library::CImg<unsigned int> image(imgPath);
	width = image.width();
	height = image.height();
	int size = width * height;

	cudaMallocHost((void**)&Datapoints, size * sizeof(Pixel));
	cimg_forXY(image,x,y) {
		Datapoints[local_id].id = local_id;
		Datapoints[local_id].closest_node = -1;
		//create RGB
		Datapoints[local_id].RGB = (float*) malloc(3* sizeof(float));
		Datapoints[local_id].R = Datapoints[local_id].RGB[0] = image(x,y,0)/255.0f;
		Datapoints[local_id].G = Datapoints[local_id].RGB[1] = image(x,y,1)/255.0f;
		Datapoints[local_id].B = Datapoints[local_id].RGB[2] = image(x,y,2)/255.0f;
		//create Lab
		std::vector <float> lab_color = getLabColor(image(x,y,0), image(x,y,1), image(x,y,2));
		Datapoints[local_id].Lab = (float*) malloc(3* sizeof(float));
		Datapoints[local_id].L = Datapoints[local_id].Lab[0] = lab_color[0];
		Datapoints[local_id].a = Datapoints[local_id].Lab[1] = lab_color[1];
		Datapoints[local_id].b = Datapoints[local_id].Lab[2] = lab_color[2];

		local_id++;
	};
	Datasize = size;
}
