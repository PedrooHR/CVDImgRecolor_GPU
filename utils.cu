
#include "utils.h"
#include <stdio.h>


float Proto[3][3] = {{0.152286, 1.052583, -0.204868},
                    {0.114503, 0.786281, 0.099216},
                    {-0.003882, -0.048116, 1.051998}};

float Deuto[3][3] = {{0.367322, 0.860646, -0.227968},
                    {0.280085, 0.672501, 0.047413},
                    {-0.011820, 0.042940, 0.968881}};

float Trito[3][3] = {{1.255528, -0.076749, -0.178779},
                    {-0.078411, 0.930809, 0.147602},
                    {0.004733, 0.691367	, 0.303900}};


#define maximo(a,b) a > b ? a : b
#define minimo(a,b) a < b ? a : b

float NDcoeff = (1.0 / sqrtf(2*M_PI*SIGMA*SIGMA));
float NDmax   = NDcoeff * pow(euler, -((MI*MI) / (2*SIGMA*SIGMA)));

std::vector <float> getRGBColor(float L, float a, float b){
	float X_ref =  94.811 / 100.0; // D65 Observer - Daylight, sRGB, Adobe-RGB
	float Y_ref = 100.000 / 100.0;
	float Z_ref = 107.304 / 100.0;

	float fy = (L + 16.0) / 116.0,
		  fx = (a / 500.0) + fy,
		  fz = fy - (b / 200.0),
		  xr, yr, zr,
		  X, Y, Z,
		  r, g, B;

	xr = (fx * fx * fx > 0.008856) ? fx * fx * fx : ((116.0 * fx) - 16.0) / 903.3;
	yr = (L > (903.3 * 0.008856)) ? fy * fy * fy : L / 903.3;
	zr = (fz * fz * fz > 0.008856) ? fz * fz * fz : ((116.0 * fz) - 16.0) / 903.3;

	X = xr * X_ref;
	Y = yr * Y_ref;
	Z = zr * Z_ref;

	r = X *  3.2406 + Y * -1.5372 + Z * -0.4986;
	g = X * -0.9689 + Y *  1.8758 + Z *  0.0415;
	B = X *  0.0557 + Y * -0.2040 + Z *  1.0570;

	std::vector <float> RGB = {r, g, B};

	return RGB;
}

std::vector <float> getLabColor(unsigned int sR, unsigned int sG, unsigned int sB)
{
	float X_ref =  94.811; // D65 Observer - Daylight, sRGB, Adobe-RGB
	float Y_ref = 100.000;
	float Z_ref = 107.304;

	float r = sR / 255.0,
		  g = sG / 255.0,
		  b = sB / 255.0,
		  x,y,z,
	      fx, fy, fz;

	//RGB to XYZ
	x = (r * 0.4124 + g * 0.3576 + b * 0.1805) / (X_ref/100.0);
	y = (r * 0.2126 + g * 0.7152 + b * 0.0722) / (Y_ref/100.0);
	z = (r * 0.0193 + g * 0.1192 + b * 0.9505) / (Z_ref/100.0);

	//XYZ to Lab
	fx = (x > 0.008856) ? pow(x, 1/3.0) : (7.787 * x) + 16.0/116.0;
	fy = (y > 0.008856) ? pow(y, 1/3.0) : (7.787 * y) + 16.0/116.0;
	fz = (z > 0.008856) ? pow(z, 1/3.0) : (7.787 * z) + 16.0/116.0;

	std::vector <float> v = {(116 * fy) - 16,
							 500 * (fx - fy),
							 200 * (fy - fz)};

	return v;
}

float sqED_2D(float x1, float y1, float x2, float y2){
	float dist = (x1 - x2) * (x1 - x2);
	dist += (y1 - y2) * (y1 - y2);
	return dist;
}


