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
	printf("\npoints\tbf\t");
	int numpoints;
    int size = 10000;
	for(numpoints = 100; numpoints<=1500; numpoints+=100){
		//int numpoints = 25;
		int endnumpoints = numpoints;
		int xvals[numpoints];
		int yvals[numpoints];
		srand(0);//time(NULL));

		// make matrix of image pixel values
		//printf("P3 800 800 1\n");

	
		// make random points	
		for (int j = 0; j < numpoints; j++){
			int randx = randVal(0,size);
			int randy = randVal(0,size);

			xvals[j] = randx;
			yvals[j] = randy;
		}
		// for(int j= 0; j<numpoints; j++){
		// 	printf("%i ",xvals[j]);
		// }


/*		// record all points of convex hull
		for (int a = 0; a < numpoints; a++){
			for (int b = 0; b < numpoints; b++){
				for (int c = 0; c < numpoints; c++){
					 for (int d = 0; d < numpoints; d++){
						if (d!=a && d!=b && d!=c){
							float alpha = (float)((yvals[b] - yvals[c])*(xvals[d] - xvals[c]) + (xvals[c] - xvals[b])*(yvals[d] - yvals[c])) / ((yvals[b] - yvals[c])*(xvals[a] - xvals[c]) + (xvals[c] - xvals[b])*(yvals[a] - yvals[c]));
							float beta = (float)((yvals[c] - yvals[a])*(xvals[d] - xvals[c]) + (xvals[a] - xvals[c])*(yvals[d] - yvals[c])) / ((yvals[b] - yvals[c])*(xvals[a] - xvals[c]) + (xvals[c] - xvals[b])*(yvals[a] - yvals[c]));
							float gamma = 1.0f - alpha - beta;
						
							if(alpha>0 && beta>0 && gamma >0){
								//check if point in final
								bool isIn = false;
								int ind = 0;
								for (int j = 0; j<endnumpoints; j++){
									if(finalxvals[j]==xvals[d] && finalyvals[j]==yvals[d]){
										isIn = true;
										ind = j;
								}}

								//update finalvals
								if (isIn){
									endnumpoints--;
									for (int e = ind; e < endnumpoints; e++){
										finalxvals[e] = finalxvals[e+1];
										finalyvals[e] = finalyvals[e+1];
									}
								}
							}
						}
					}
				}
			}
		}

		int hull[2][endnumpoints];
		for (int l = 0; l<endnumpoints; l++){
			hull[0][l] = finalxvals[l];
			hull[1][l] = finalyvals[l];
		}

		//sort hull by x value
		int it1, it2, min_ind;
		for (it1 = 0; it1<endnumpoints-1; it1++){
			min_ind = it1;
			for (it2 = it1+1; it2 < endnumpoints; it2++){
				if (hull[0][it2] < hull[0][min_ind]){
					min_ind = it2;
				}
			}
			//swap
			int temp = hull[0][min_ind];
			hull[0][min_ind] = hull[0][it1];
			hull[0][it1] = temp;

			int temp2 = hull[1][min_ind];
			hull[1][min_ind] = hull[1][it1];
			hull[1][it1] = temp2;
		}

		//connect top part of hull
		int pastind = 0;
		int starty = hull[1][0];
		for (int o = 0; o<endnumpoints; o++){
			if (hull[1][o]>starty){
				//DRAWLINE
				int x1 = hull[0][pastind];
				int x2 = hull[0][o];
				int y1 = hull[1][pastind];
				int y2 = hull[1][o];

				double trials = 1000.0;
				double xslope = (x2-x1)/trials;
				double yslope = (y2-y1)/trials;
				for (int iter = 0; iter < 1000; iter++){
					image[(int)(x1+xslope*iter)][(int)(y1+yslope*iter)].set_values(0,0,0);
				}
				//DRAWLINE
				pastind = o;
			}
		}
		//connect ends
		//DRAWLINE
		int x1 = hull[0][pastind];
		int x2 = hull[0][endnumpoints-1];
		int y1 = hull[1][pastind];
		int y2 = hull[1][endnumpoints-1];

		double trials = 1000.0;
		double xslope = (x2-x1)/trials;
		double yslope = (y2-y1)/trials;
		for (int iter = 0; iter < 1000; iter++){
			image[(int)(x1+xslope*iter)][(int)(y1+yslope*iter)].set_values(0,0,0);
		}
		//DRAWLINE

		//sort hull by x value BACKWARDS
		int max_ind;
		for (it1 = 0; it1<endnumpoints-1; it1++){
			max_ind = it1;
			for (it2 = it1+1; it2 < endnumpoints; it2++){
				if (hull[0][it2] > hull[0][max_ind]){
					max_ind = it2;
				}
			}
			//swap
			int temp = hull[0][max_ind];
			hull[0][max_ind] = hull[0][it1];
			hull[0][it1] = temp;
		
			int temp2 = hull[1][max_ind];
			hull[1][max_ind] = hull[1][it1];
			hull[1][it1] = temp2;
		}

		//connect top part of hull
		int pastind2 = 0;
		starty = hull[1][0];
		for (int o = 0; o<endnumpoints; o++){
			if (hull[1][o]<starty){
				//DRAWLINE
				int x1 = hull[0][pastind2];
				int x2 = hull[0][o];
				int y1 = hull[1][pastind2];
				int y2 = hull[1][o];

				double trials = 1000.0;
				double xslope = (x2-x1)/trials;
				double yslope = (y2-y1)/trials;
				for (int iter = 0; iter < 1000; iter++){
					image[(int)(x1+xslope*iter)][(int)(y1+yslope*iter)].set_values(0,0,0);
				}
				//DRAWLINE
				pastind2 = o;
			}
		}
		//connect ends
		//DRAWLINE
		int Bx1 = hull[0][pastind2];
		int Bx2 = hull[0][endnumpoints-1];
		int By1 = hull[1][pastind2];
		int By2 = hull[1][endnumpoints-1];

		double Btrials = 1000.0;
		double Bxslope = (Bx2-Bx1)/Btrials;
		double Byslope = (By2-By1)/Btrials;
		for (int iter = 0; iter < 1000; iter++){
			image[(int)(Bx1+Bxslope*iter)][(int)(By1+Byslope*iter)].set_values(0,0,0);
		}
		//DRAWLINE
	

		// print to ppm file using ./a.out > assignment1.ppm
		for (int i = 0; i < size; i++) {
			for (int i2=0;i2<size;i2++){
				;//printf("%d %d %d ", image[i][i2].getR(), image[i][i2].getG(), image[i][i2].getB());
			}
		}
*/		int start_s=clock();/////time it
		/*Point P[numpoints];
		for (int i = 0; i<numpoints; i++){
			P[i].x = xvals[i];
			P[i].y = yvals[i];
		}
		int n = sizeof(P) / sizeof(P[0]); 
		printf("The smallest distance is %f \n", closest(P, n)); */
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

        printf("Distance: %f \n", min);

		int stop_s=clock();/////time it
		printf("%i\t%f\t",numpoints,(stop_s-start_s)/double(CLOCKS_PER_SEC));/////

	}
	printf("\n\n");
	return 0;
}
