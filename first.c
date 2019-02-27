#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tga.h"
#include "model.h"

int sign (int x){
if (x>0) {
	return 1; 
}
if (x<0) {
	return -1;
}
return 0;
}

void swap(int *x, int *y){
int z = *x;
*x = *y;
*y = z;
}

void line(tgaImage *image, int x0,int y0,int x1, int y1){
tgaColor color = tgaRGB(255, 255, 255);
int steep=0;

if (abs(y1-y0)>abs(x1-x0)) {
	steep=1;
	swap(&x0,&y0);
	swap(&x1,&y1);
} 
else {
	steep=0;
} 

if (x0>x1) {
	swap(&x0,&x1);
	swap(&y0,&y1);
}

int dx=x1-x0;
int dy=y1-y0;
int de=2*abs(dy);
int e=0;

int y=y0;
for (int x=x0; x<x1; ++x){
	if (steep==1) {
		tgaSetPixel(image, y, x, color);
	} 
	else {
		tgaSetPixel(image, x, y, color);	
	}
	e=e+de;
	if (e>dx) {
		y=y+sign(dy);
		e=e-2*dx;
	}
}
}

void meshgrid(tgaImage *image, Model *model) {
double max_x = -1000;
double max_y = -1000;
double min_x = 1000;
double min_y = 1000;

for (int i = 0; i < model->nvert; i++){
	double x = (model->vertices[i][0]);
	double y = (model->vertices[i][1]);

	if (x > max_x) {
		max_x = x;
	}
	if (y > max_y) {
		max_y = y;
	}
	if (x < min_x) {
		min_x = x;
	}
	if (y < min_y) {
		min_y = y;
	}
}

// Проверка нахождения минимумов и максимумов
if ((max_x == 0) && (min_x == 0) && (max_y == 0) && (min_y == 0)) {
	line(image, 0, 400, 800, 400);
} 
else {
	double a, b;
	double add_y = 0;
	double add_x = 0;

	if ((max_x - min_x) > (max_y - min_y)) {
		a = max_x;
		b = min_x;
		add_y = ((max_y+min_y) - b - a)/(a - b); //(1-((max_y - min_y) - b - a)/(a - b)); 
	} 
	else {
		a = max_y;
		b = min_y;
		add_x = -((max_x+min_x) - b - a)/(a - b);  
	}

	for (int i = 0; i < model->nface; i++) {
		int screen_coords[3][2];
	
		for (int j = 0; j < 3; j++) {
			Vec3 *v = &(model->vertices[model->faces[i][3*j]]);

			screen_coords[j][0] = (((2*(*v)[0] - b - a)/(a - b)+ 1 + add_x) * image->width / 2);   //x
			screen_coords[j][1] = (1 - (2*(*v)[1] - b - a)/(a - b) + add_y) * image->height / 2;  //y
		}

		for (int j = 0; j < 3; ++j){
			line(image, screen_coords[j][0], screen_coords[j][1], screen_coords[(j+1)%3][0], screen_coords[(j+1)%3][1]);
		}
	}
}
}

int main(int argc, char const *argv[], char const *argv_1[]) {

// Загружаем модель из файла
Model *model = loadFromObj(argv[1]);

// Создаем область для рисования
int height = 1020;
int width = 1020;

tgaImage * image = tgaNewImage(height, width, RGB);
if (image == NULL) {
    fprintf(stderr, "%s\n", "Failed to create image");
    return EXIT_FAILURE;
}

// Рисуем сетчатую модель
meshgrid(image, model);

if (-1 == tgaSaveToFile(image, argv_1[1])) {
    perror("tgaSateToFile");
    return EXIT_FAILURE;
}

tgaFreeImage(image);
freeModel(model);

return EXIT_SUCCESS;
}
