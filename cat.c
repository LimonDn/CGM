#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "tga.h"
#include "model.h"

const int height = 1020;
const int width = 1020;
const int deph = 255;
const int phi = 30;


void swap (int *x, int *y){
int z = *x;
*x = *y;
*y = z;
}

void swapf (double *x, double *y){
double z = *x;
*x = *y;
*y = z;
}

void transform (int n, int m, double *A, double *x, double *y){
  for (int i=0; i<n; ++i){
    y[i] = 0;
    for (int j=0; j<m; ++j){
      y[i]+=A[i*m+j]*x[j];
    }
  }
} 

void triangle (Model *model, tgaImage *image, int x0, int y0, int z0, int x1, int y1, int z1, int x2, int y2, int z2, 
                       Vec3 uv0, Vec3 uv1, Vec3 uv2, int zbuf[image->width][image->height], tgaColor color,float I){
  
  int a, b;
  int mas[2][3] = { {x0,x1,x2}, {y0,y1,y2}};
  double uv[3][3] = {{uv0[0],uv0[1],uv0[2]}, {uv1[0],uv1[1],uv1[2]}, {uv2[0],uv2[1],uv2[2]}};
  int z, za, zb;

  for (int i = 0; i < 3; ++i){
    for (int j = 1; j < 3; ++j) {
      if (mas[1][j-1]>mas[1][j]) {
        swap(&mas[1][j-1], &mas[1][j]);
        swap(&mas[0][j-1], &mas[0][j]);
        for (int k=0; k<3; ++k){
          swapf(&uv[j-1][k], &uv[j][k]);
        }
        }
      }
  }

  if (((x0 == x1) && ( x0 == x2))||((y0 == y1) && ( y0 == y2))) {
    if ((y0 == y1) && ( y0 == y2)) {
      for (int x = mas[0][0]; x < mas[0][2]; ++x){
        Vec3 U;
        z = (x-mas[0][0])*(z2-z0)/(mas[0][2]-mas[0][0])+z0;
        for (int k=0; k<3; ++k){
          U[k] = (x-mas[0][0])*(uv[2][k]-uv[0][k])/(mas[0][2]-mas[0][0])+uv[0][k];
        }
        if (z > zbuf[x][y0]) {
          zbuf[x][y0] = z;
          tgaColor color1;
          color1 = getDiffuseColor(model, U);
          tgaSetPixel(image, x, y0, tgaRGB(fabs(I)*Red(color1), fabs(I)*Green(color1), fabs(I)*Blue(color1)));
        }
      }
    } else {
      for (int y = mas[1][0]; y < mas[1][2]; ++y){
        Vec3 U;
        z = (y-mas[1][0])*(z2-z0)/(mas[1][2]-mas[1][0])+z0;
        for (int k=0; k<3; ++k){
          U[k] = (y-mas[1][0])*(uv[2][k]-uv[0][k])/(mas[1][2]-mas[1][0])+uv[0][k];
        }
        if (z > zbuf[x0][y]) {
          zbuf[x0][y] = z;
          tgaColor color1;
          color1 = getDiffuseColor(model, U);
          tgaSetPixel(image, x0, y, tgaRGB(fabs(I)*Red(color1), fabs(I)*Green(color1), fabs(I)*Blue(color1)));
        }
      }
    }
  } 
  else {

  for (int y = mas[1][0]; y < mas[1][2]; ++y) {
    Vec3 uvA, uvB;

    if ((y > mas[1][1]) || (mas[1][0]==mas[1][1])) {
      a = (y-mas[1][1])*(mas[0][2]-mas[0][1])/(mas[1][2]-mas[1][1])+mas[0][1];
      za = (y-mas[1][1])*(z2-z1)/(mas[1][2]-mas[1][1])+z1;
      for (int k = 0; k < 3; ++k){
        uvA[k] = uv[1][k] + (uv[2][k]-uv[1][k])*(y-mas[1][1])/(mas[1][2]-mas[1][1]);
      }
    } else {
      a = (y-mas[1][0])*(mas[0][1]-mas[0][0])/(mas[1][1]-mas[1][0])+mas[0][0];
      za = (y-mas[1][0])*(z1-z0)/(mas[1][1]-mas[1][0])+z0;
      for (int k = 0; k < 3; ++k){
        uvA[k] = uv[0][k] + (uv[1][k]-uv[0][k])*(y-mas[1][0])/(mas[1][1]-mas[1][0]);
      }
    }

    for (int k = 0; k < 3; ++k){
      uvB[k] = uv[0][k] + (uv[2][k]-uv[0][k])*(y-mas[1][0])/(mas[1][2]-mas[1][0]);
    }
    b = (y-mas[1][0])*(mas[0][2]-mas[0][0])/(mas[1][2]-mas[1][0])+mas[0][0];
    zb = (y-mas[1][0])*(z2-z0)/(mas[1][2]-mas[1][0])+z0;

    if (b < a) {
      swap(&a,&b);
      for (int k=0; k<3; ++k){
        swapf(&uvA[k], &uvB[k]);
      }
    }

    for (int x = a; x < b; ++x){
      Vec3 U;
      if ((b-a)!=0){
        z = (x-a)*(zb-za)/(b-a)+za;
        for (int k=0; k<3; ++k){
          U[k] = ((x-a)*(uvB[k]-uvA[k])/(b-a)+uvA[k]);
        }
        if (z > zbuf[x][y]) {
          zbuf[x][y] = z;
          tgaColor color1;
          color1 = getDiffuseColor(model, U);
          tgaSetPixel(image, x, y, tgaRGB(fabs(I)*Red(color1), fabs(I)*Green(color1), fabs(I)*Blue(color1)));
        } 
      }
    }
  }

}
}

void meshgrid(tgaImage *image, Model *model) {
double max_x = DBL_MIN;
double max_y = DBL_MIN;
double min_x = DBL_MAX;
double min_y = DBL_MAX;
double A[16]= {cos(phi*3.14/180),0,-sin(phi*3.14/180),0, 0,1,0,0, sin(phi*3.14/180),0,cos(phi*3.14/180),0, 0,0,0,1};

loadDiffuseMap(model, "cat_d.tga");

int zbuf[image->width][image->height];
for (int i = 0; i < image->width; ++i){
  for (int j = 0; j < image->height; ++j) {
    zbuf[i][j] = -1000000;;
  }
}

for (int i = 0; i < model->nvert; i++){
  double rx[4] = {(model->vertices[i][0]),(model->vertices[i][1]),(model->vertices[i][2]),1};
  double coord[4] ={};
  transform(4,4,&A, &rx, &coord);

  if (coord[0] > max_x) {
    max_x = coord[0];
  }
  if (coord[1] > max_y) {
    max_y = coord[1];
  }
  if (coord[0] < min_x) {
    min_x = coord[0];
  }
  if (coord[1] < min_y) {
    min_y = coord[1];
  }
}

// Проверка нахождения минимумов и максимумов
if ((max_x == 0) && (min_x == 0) && (max_y == 0) && (min_y == 0)) {
  printf("ERROR");
} 
else {
  double a, b;
  double add_y = 0;
  double add_x = 0;

  if ((max_x - min_x) > (max_y - min_y)){
    a = max_x;
    b = min_x;
    if ((a-b)!=0) {
      add_y = ((max_y+min_y) - b - a)/(a - b); 
    } else {
      add_y = ((max_y+min_y) - b - a);
    }
  } 
  else {
    
    a = max_y;
    b = min_y;
    if ((a-b)!=0) {
      add_x = -((max_x + min_x) - b - a)/(a - b);  
    } else {
      add_x = -((max_x + min_x) - b - a);
    }
  }

  for (int i = 0; i < model->nface; i++) {
    int screen_coords[3][3];
    double p[3][4] = {};
    Vec3 *uv[3];

    for (int j = 0; j < 3; j++) {
      Vec3 *v = &(model->vertices[model->faces[i][3*j]]);
      uv[j] = (getDiffuseUV(model,i,j));
      double rx[4] = {(*v)[0],(*v)[1],(*v)[2],1};

      transform(4,4,&A, &rx, &p[j]);

      screen_coords[j][0] = (((2*p[j][0] - b - a)/(a - b)+ 1 + add_x) * image->width / 2);   //x
      screen_coords[j][1] = (1 - (2*p[j][1] - b - a)/(a - b) + add_y) * image->height / 2;  //y
      screen_coords[j][2] = (p[j][2]+1) * deph / 2;
    }
    double a[3], b[3];
    for (int i = 0; i<3; i++){
      a[i] = p[1][i]-p[0][i];
      b[i] = p[2][i]-p[0][i];
    }

    double n[3];
    int light [3] = {0,0,-1};

    n[0] = a[1]*b[2]-a[2]*b[1];
    n[1] = -(a[0]*b[2]-a[2]*b[0]);
    n[2] = a[0]*b[1]-a[1]*b[0];

    double norm = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);    
    double I = (light[0]*n[0]+light[1]*n[1]+light[2]*n[2])/norm;
    
    if (I < 0) {
      triangle(model, image, screen_coords[0][0], screen_coords[0][1], screen_coords[0][2], 
            screen_coords[1][0], screen_coords[1][1], screen_coords[1][2], 
              screen_coords[2][0], screen_coords[2][1], screen_coords[2][2], uv[0], uv[1], uv[2],
                  zbuf, tgaRGB(fabs(I)*255, fabs(I)*255, fabs(I)*255), I);
    }
  }
}
}

int main(int argc, char const *argv[]) {
Model *model = loadFromObj(argv[1]);
tgaImage * image = tgaNewImage(height, width, RGB);

meshgrid(image, model);

if (-1 == tgaSaveToFile(image, "128.tga")) {
    perror("tgaSateToFile");
    return EXIT_FAILURE;
}

tgaFreeImage(image);
freeModel(model);

return EXIT_SUCCESS;
}