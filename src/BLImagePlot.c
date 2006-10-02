#include <math.h>
#include <R.h>
#include "beadarray.h"

void BLImagePlot(int *nProbes, double *fground, double *yvalues, int *ygrid, double *result, int *nrow){

    int i,j;
	double counter, total;
    for(i = 0; i < (*nrow); i++){
	counter = 0;
	total = 0;

        for(j = 0; j < (*nProbes); j++){

            if((yvalues[j] >= ygrid[i]) & (yvalues[j] < ygrid[i+1])){
				if(fground[j] > 0){
                    total = total + (log(fground[j])/log(2.00));
	    			counter = counter + 1;
				}
			}
		}
	//Error checking. Total can become NaN if there was a negative value in fground.
	//Hopefully this shouldn't ever print out.
		if(isnan(total) == 1){
		    Rprintf("Total = NaN\n");
		}
		result[i] = total/counter;
	}
}
