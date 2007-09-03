#include "beadarray.h"
#include <stdlib.h>
#include <math.h>
#include <R.h>


int binarySearch(int key, int *arr, int left, int right){
	int mid;
	
	if (right <= left) {
		return -1;
	}
	
	mid = floor((right-left)/2)+left;
	if(key > arr[mid]){
		return binarySearch(key, arr, mid + 1, right);
	}
	else if (key < arr[mid]){
		return binarySearch(key, arr, left, mid - 1);
	}
	else {
	 	return mid;
	}
}

double sd(double *arr, int length, double mean){

    int i;
	double sum, temp;

	sum = 0;

	for(i = 0; i < length; i++){
	    temp = arr[i] - mean;
		sum += pow(temp, 2);
		}

	temp = sum/(length-1);
	return(sqrt(temp));
}

void quicksortDouble(double arr[], int low, int high) {
 int i = low;
 int j = high;
 double y = 0;
 /* compare value */
 double z = arr[(low + high) / 2];

 /* partition */
 do {
  /* find member above ... */
  while(arr[i] < z) i++;

  /* find element below ... */
  while(arr[j] > z) j--;

  if(i <= j) {
   /* swap two elements */
   y = arr[i];
   arr[i] = arr[j]; 
   arr[j] = y;
   i++; 
   j--;
  }
 } while(i <= j);

 /* recurse */
 if(low < j) 
  quicksortDouble(arr, low, j);

 if(i < high) 
  quicksortDouble(arr, i, high); 
}

double median(double *arr, int length){
	int i;
	double j;
	double *tempArr;

	tempArr = (double *)malloc(sizeof(double) * length);

	for(i = 0; i < length; i++){
		  tempArr[i] = arr[i];
	}

	quicksortDouble(tempArr, 0, (length-1));
	i = floor(length/2);

	if(i == (length/2.0)){ //if its of even length take the average of the middle two
		j = (tempArr[i-1]+tempArr[i])/2;
	}
	else { //if odd length
		j = tempArr[i];
	}

	free(tempArr);
	return j;
}

double mad(double *arr, int length){
	int i;
	double temp, med;
	double *arr2;

	med = median(arr, length);

	arr2 = (double *)malloc(sizeof(double) * length);

	for(i = 0;i < length; i++){
		  temp = (arr[i] - med);
		  arr2[i] = fabs(temp);
	}
	temp = 1.4826*(median(arr2, length));
	free(arr2);
	return(temp);
}

int* getProbeIndices(int *probeList, int probeID, int *start, int numBeads){
	int i , j, ind;
	int *indices;
	indices = (int *)malloc(sizeof(int) * 2);
	/* find the location of the one probe of the correct ID */
	ind = binarySearch(probeID, probeList, 0, numBeads);
	/*check if that probe ID was found */
	if(ind != -1){
	
		i = j = ind;
		/* go left until you stop finding that probeID and record the index of the last one */
		while((i >= 0) && (probeList[i] == probeID)){
	 		i--;
        	}
		indices[0] = i+1;

		/* go right and do the same */
		while((j < numBeads) && (probeList[j] == probeID)){
	 		j++;
 		}
		indices[1] = j-1;
	}
	else {
		Rprintf("ProbeID %d not found\n", probeID);
		indices[0] = 1;
		indices[1] = 0;
	}

	return(indices);
}

beadStatusStruct* findBeadStatus(double *intensities, int *probeList, int probeID, int numBeads, int *count, int *start, double *nmads){
	beadStatusStruct *status;
	int *outlierInds, *validInds, *indices;
	double *inten;
	double m, ma;
	int i,j,k, nsize1, nsize2;

	//get the start and end indices for this probeID
	indices = getProbeIndices(probeList, probeID, start, numBeads);

	*start = (indices[1]+1);  //set the next start position
	*count = (indices[1] - indices[0] + 1);  //calculate how many beads there are of this type
	inten = (double *)malloc(sizeof(double) * (*count));

	i = indices[0];
	j = indices[1];
	k = 0;
	
	while(i < (j+1)){
		inten[k] = intensities[i];
		i++;
		k++;
	}
	m = median(inten, (*count));
	ma = mad(inten, (*count));

	status = (beadStatusStruct *)malloc(sizeof(beadStatusStruct));
	validInds = (int *)malloc(sizeof(int));
	outlierInds = (int *)malloc(sizeof(int));

	i = j = k = nsize1 = nsize2 = 0;	

 
	/* find the indices of valid probes */
	while(k < (*count)){
	    if(ma==0) {
		validInds[i] = (indices[0]+k);
		i++;
		nsize1 = (sizeof(int) * (i+1));
		validInds = (int *)realloc(validInds, nsize1);	
	    }
            else {
		if((inten[k] < (m + *nmads*ma)) && (inten[k] > (m - *nmads*ma))){
//					 Rprintf("%d - valid\n",k);
//					 Rprintf("i: %d, Index: %d\n", i, (indices[0]+k));
			validInds[i] = (indices[0]+k);
			i++;
			nsize1 = (sizeof(int) * (i+1));
			validInds = (int *)realloc(validInds, nsize1);	
		}
		else{
//			 Rprintf("%d - outlier\n", k);
			outlierInds[j] = (indices[0]+k);
			j++;
			nsize2 = (sizeof(int) * (j+1));
			outlierInds = (int *)realloc(outlierInds, nsize2);
		}
	    }
		k++;
	}

	validInds[i] = -1;
	outlierInds[j] = -1;

	status->validInds = validInds;
	status->outlierInds = outlierInds;

	free(inten);
	free(indices);
	return(status);
}

void findAllOutliers(double *finten, int *binaryStatus, int *probeList, int *probeIDs, int *numProbes, int *numBeads, int *nextStart, double *nmads){

	beadStatusStruct *status; 
	int *valids, *count;
	int k, i, probeID, temp;

	count = (int *)malloc(sizeof(int));
	*count = 0; //initialise counter to record how many beads there are of this type

	for(k = 0; k < (*numProbes); k++){
		  
		  probeID = probeIDs[k];
		  status = findBeadStatus(finten, probeList, probeID, *numBeads, count, nextStart, nmads);
		  valids = status->validInds;
//		  outliers = status->outlierInds;

		  i = 0;
		  while(valids[i] != -1){
		  		temp = valids[i];
		  		binaryStatus[temp] = 1;
				i++;
				}
		}
	free(count);
	free(status->validInds);
	free(status->outlierInds);
	free(status);
}

void createBeadSummary(double *finten, double *binten, int *probeList, int *probeIDs, int *numProbes, int *numBeads, double *foreground, 
	 						  double *background, double *stdev, int *numValid, int *numOutlier, int *nextStart, double *nmads,  int *method, double *trimfrac){
	 
	 beadStatusStruct *status;
	 int *valids, *count, i, j, k, l, ind, probeID, *indices, *start, low, high;
	 double *validInten, *trimInten; 
	 double fground, bground;

         start = (int *)malloc(sizeof(int));
	 count = (int *)malloc(sizeof(int));
	 *count = 0; //initialise counter to record how many beads there are of this type
         *start = 0;

         if(*method==1) { // Illumina default method
  	   for(k = 0; k < (*numProbes); k++){
	      probeID = probeIDs[k];

              status = findBeadStatus(finten, probeList, probeID, *numBeads, count, nextStart, nmads);

	      valids = status->validInds;

  	      i=0;
              fground = bground = 0.0;
              while(valids[i] != -1){
        	  ind = valids[i];
	 	  fground += finten[ind];
		  bground += binten[ind];
		  i++;
	      }

  	      foreground[k] = (fground/i);
	      background[k] = (bground/i);

  	      validInten = (double *)malloc(sizeof(double) * i);

	      i = 0;
	      while(valids[i] != -1){
	 	  ind = valids[i];
	 	  validInten[i] = finten[ind];
		  i++;
	      }
	      stdev[k] = sd(validInten, i, foreground[k]);
	      numValid[k] = i;
	      numOutlier[k] = ((*count)-i);

	      free(validInten);
//	      free(valids);
	      free(status->validInds);
	      free(status->outlierInds);
	      free(status);
           }
	 }

	 if(*method==2) { // regular mean and sd  	   
           for(k = 0; k < (*numProbes); k++){
              probeID = probeIDs[k];
	      indices = getProbeIndices(probeList, probeID, start, *numBeads);
              i = indices[0];
	      j = indices[1];
              *start = (indices[1]+1); 

              validInten = (double *)malloc(sizeof(double) * (j-i+1));	      
              fground = 0.0;
              l = 0;

	      while(i < (j+1)) {
	 	  fground += finten[i];
                  validInten[l] = finten[i];
 		  i++;
                  l++;
	      }
              foreground[k] = (fground/l);
	      background[k] = 0;

              stdev[k] = sd(validInten, l, foreground[k]);
              numValid[k] = l;
              numOutlier[k] = ((*count)-l);
              free(validInten);
              free(indices);
	    }
	 }
	 if(*method==3) { // trimmed mean and sd
           for(k = 0; k < (*numProbes); k++){
              probeID = probeIDs[k];
              low = 0, high=0;
  	      indices = getProbeIndices(probeList, probeID, start, *numBeads);              
              i = indices[0];
	      j = indices[1];
              *start = (indices[1]+1);
	      validInten = (double *)malloc(sizeof(double) * (j-i+1));

              l = 0;
	      while(i < (j+1)){
                  validInten[l] = finten[i];
		  i++;
		  l++;
	      }
              quicksortDouble(validInten, 0, l-1);
              low = round(*trimfrac*(l-1)); /* floor ceil */
              high = round((1-*trimfrac)*(l-1)); /* floor ceil */
              trimInten = (double *)malloc(sizeof(double) * (high-low+1));
              fground = 0.0;
              l = 0;
              for(i = low; i <= high; i++){
	 	fground += validInten[i];
                trimInten[l] = validInten[i];
                l++;
	      }
              foreground[k] = (fground/l);
	      background[k] = 0;
              stdev[k] = sd(trimInten, l, foreground[k]);
              numValid[k] = l;
              numOutlier[k] = ((*count)-l);
              free(validInten);
              free(trimInten);
              free(indices);
	   }
	 }
         if(*method==4) { // winsorized mean and sd           
           for(k = 0; k < (*numProbes); k++){
              probeID = probeIDs[k];
              low = 0, high=0; 
	      indices = getProbeIndices(probeList, probeID, start, *numBeads);
              i = indices[0];
	      j = indices[1];
              *start = (indices[1]+1);
	      l = 0;
	      validInten = (double *)malloc(sizeof(double) * (j-i+1));

	      while(i < (j+1)){
                  validInten[l] = finten[i];
		  i++;
		  l++;
	      }
              quicksortDouble(validInten, 0, l-1);
              low = round(*trimfrac*(l-1)); /* floor ceil */
              high = round((1-*trimfrac)*(l-1)); /* floor ceil */
              i = 0;
              while(i<low) {
	          validInten[i] = validInten[low];
                  i++;
              }
              i = high+1;
              while(i<l) {
	          validInten[i] = validInten[high];
                  i++;
	      }
              fground = 0.0;
              for(i = 0; i < l; i++){
	 	 fground += validInten[i];
	      }
              foreground[k] = fground/l;
	      background[k] = 0;
              stdev[k] = sd(validInten, l, foreground[k]);
              numValid[k] = l;
              numOutlier[k] = ((*count)-l);
              free(validInten);
              free(indices);
	    }
	 }
         if(*method==5) { // median and mad           
           for(k = 0; k < (*numProbes); k++){
              probeID = probeIDs[k];
	      indices = getProbeIndices(probeList, probeID, start, *numBeads);
              i = indices[0];
	      j = indices[1];
              *start = (indices[1]+1); 

              validInten = (double *)malloc(sizeof(double) * (j-i+1));	      
              fground = 0.0;
              l = 0;
	      while(i < (j+1)) {
                  validInten[l] = finten[i];
 		  i++;
                  l++;
	      }
              foreground[k] = median(validInten, l);
	      background[k] = 0;

              stdev[k] = mad(validInten, l);
              numValid[k] = l;
              numOutlier[k] = 0;
              free(validInten);
              free(indices);
	    }
	 }
   free(count);
   free(start);
}
