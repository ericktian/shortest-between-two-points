#include <iostream>
#include <stdio.h>
#include <string>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <float.h>
#include <stdlib.h>
#include <time.h>

#define pi 3.14159265358979323846264338

using namespace std;

class rgb {
	private:
		int r,g,b;
	public:
		rgb()		{ r=1;g=1;b=1; }
		rgb(int Rin, int Gin, int Bin)	{ r=Rin;g=Gin;b=Bin; }
		void set_values(int,int,int);
		int getR();
		int getG();
		int getB();
};
void rgb::set_values (int Rin, int Gin, int Bin) {
	r = Rin;
	g = Gin;
	b = Bin;
}
int rgb::getR(){
	return r;
}
int rgb::getG(){
	return g;
}
int rgb::getB(){
	return b;
}
int randVal(int min, int max){
	return min + (rand() % static_cast<int>(max - min + 1));
}

int main(void){
	int numpoints = 30;
	int endnumpoints = numpoints;
	int xvals[numpoints];
	int yvals[numpoints];
	int finalxvals[endnumpoints];
	int finalyvals[endnumpoints];
	int size = 800;
	rgb image[size][size];
	srand(0);//time(NULL));

	// make matrix of image pixel values
	printf("P3 800 800 1\n");

	
	// make random points	
	for (int j = 0; j < numpoints; j++){
		int randx = randVal(0,size);
		int randy = randVal(0,size);

		xvals[j] = randx;
		yvals[j] = randy;
		finalxvals[j] = randx;
		finalyvals[j] = randy;
		image[xvals[j]][yvals[j]].set_values(0,0,0);
	}


	float min = float(sqrt(2*(size*size)));
        int minind = 0;
        int minind2 = 1;

        for(int i = 0; i<numpoints; i++){
                for(int j = 0; j<numpoints; j++){
                    if(i!=j){// and xvals[i]!=xvals[j] and yvals[i]!=yvals[j]){
                        float curdist = sqrt(pow(float(xvals[i]-xvals[j]),2)+pow(float(yvals[i]-yvals[j]),2));
						// if(curdist==0){
						// 	printf("%i %i %i %i %i %i\n",i,j,xvals[i],xvals[j],yvals[i],yvals[j]);
						// }
                        if(curdist<min){
                            min = curdist;
                            minind = i;
                            minind2 = j;
                        }
                    }
                }
        }

	//DRAWLINE
	int x1 = xvals[minind];
	int x2 = xvals[minind2];
	int y1 = yvals[minind];
	int y2 = yvals[minind2];

	double trials = 1000.0;
	double xslope = (x2-x1)/trials;
	double yslope = (y2-y1)/trials;
	for (int iter = 0; iter < 1000; iter++){
		image[(int)(x1+xslope*iter)][(int)(y1+yslope*iter)].set_values(0,0,0);
	}
	//DRAWLINE

	

	// print to ppm file using ./a.out > assignment1.ppm
	for (int i = 0; i < size; i++) {
		for (int i2=0;i2<size;i2++){
			printf("%d %d %d ", image[i][i2].getR(), image[i][i2].getG(), image[i][i2].getB());
		}
	}
	return 0;
}
