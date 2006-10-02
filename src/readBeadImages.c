#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include "beadarray.h"

#define B0 256

void quicksort(int arr[], int low, int high) {
 int i = low;
 int j = high;
 int y = 0;
 /* compare value */
 int z = arr[(low + high) / 2];

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
  quicksort(arr, low, j);

 if(i < high) 
  quicksort(arr, i, high); 
}

       
int number(int l, char c[])
{
   int n, b, i;

   n = 0;
   b = 1;
   for (i = 0; i < l; i++) {
      n += b*((unsigned int)c[i]%B0);
      b *= B0;
   }
   return(n);
}

void sharpen(int **pixels, int ImageWidth, int ImageHeight){
     
     int **sharpened;
     int i, j, k;
     
     sharpened = malloc(sizeof(int *) * ImageWidth);
     for(i = 0; i < ImageWidth; i++){
           sharpened[i] = malloc(sizeof(int) * ImageHeight);
           }

     for(j = 1; j < ImageHeight-1; j++){
           for(k = 1; k < ImageWidth-1; k++){
                 sharpened[k][j] = pixels[k][j]-(0.5*(pixels[k][j+1]+pixels[k][j-1]+pixels[k-1][j]+pixels[k+1][j]-(4*pixels[k][j])));
           }
     }

     for (i = 0; i < ImageWidth; i++){
         for (j = 0; j < ImageHeight; j++) {
             pixels[i][j] = sharpened[i][j];    
         }
     }
     
     /*free the memory */
     for (i = 0; i < ImageWidth; i++){
             free(sharpened[i]);   
         }
     free(sharpened);
     
}



void asf(int **pixels, int ImageWidth, int ImageHeight)
#define SE 4 /* max half-size of the square structuring element */
{
    int i, j, l, m, n, ii, jj, s;
    int **x, **w;
     
    x = malloc(sizeof(int *) * ImageWidth);
    w = malloc(sizeof(int *) * ImageWidth);
     
    for(i = 0;i < ImageWidth; i++){
          x[i] = malloc(sizeof(int) * ImageHeight);
          w[i] = malloc(sizeof(int) * ImageHeight);
          }
    
       /* now filter pepper noise */
   Rprintf("Filtering pepper noise...\n");
   for (j = 1; j < ImageHeight-1; j++)
      for (i = 1; i < ImageWidth-1; i++) {
         m = pixels[i][j];
         n = pixels[i][j-1];
         pixels[i][j] = n;
         for (ii = i-1; ii <= i+1; ii++)
            for (jj = j-1; jj <= j+1; jj++)
               if (pixels[ii][jj] < n) n = pixels[ii][jj]; 
         pixels[i][j] = (m > n) ? m : n; /* max(m, n) */
      }

   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth; i++) {
         w[i][j] = pixels[i][j];
         x[i][j] = w[i][j];
      }

   /* now do the first step with a 2x2 structuring element */
   Rprintf("Opening and closing with a 2x2 element...\n");

   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth-1; i++) 
         x[i][j] = (w[i][j] > w[i+1][j]) ? w[i][j] : w[i+1][j];
   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth-1; i++) 
         w[i][j] = (x[i][j] < x[i+1][j]) ? x[i][j] : x[i+1][j];
   for (i = 0; i < ImageWidth; i++) 
      for (j = 0; j < ImageHeight-1; j++)
         x[i][j] = (w[i][j] > w[i][j+1]) ? w[i][j] : w[i][j+1];
   for (i = 0; i < ImageWidth; i++) 
      for (j = 0; j < ImageHeight-1; j++)
         w[i][j] = (x[i][j] < x[i][j+1]) ? x[i][j] : x[i][j+1];

   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth-1; i++) 
         x[i][j] = (w[i][j] < w[i+1][j]) ? w[i][j] : w[i+1][j];
   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth-1; i++) 
         w[i][j] = (x[i][j] > x[i+1][j]) ? x[i][j] : x[i+1][j];
   for (i = 0; i < ImageWidth; i++) 
      for (j = 0; j < ImageHeight-1; j++)
         x[i][j] = (w[i][j] < w[i][j+1]) ? w[i][j] : w[i][j+1];
   for (i = 0; i < ImageWidth; i++) 
      for (j = 0; j < ImageHeight-1; j++)
         w[i][j] = (x[i][j] > x[i][j+1]) ? x[i][j] : x[i][j+1];

   /* now do the alternating sequential filtering */
   for (s = 1; s <= SE; s++) {
      Rprintf("Calculating ASF with s = %i...\n", s);

      for (j = 0; j < ImageHeight; j++)
         for (i = s; i < ImageWidth-s; i++) 
            for (l = -s+1; l <= s; l++)
               x[i][j] = (w[i-s][j] > w[i+l][j]) ? w[i-s][j] : w[i+l][j];
      for (j = 0; j < ImageHeight; j++)
         for (i = s; i < ImageWidth-s; i++)
            for (l = -s+1; l <= s; l++)
               w[i][j] = (x[i-s][j] < x[i+l][j]) ? x[i-s][j] : x[i+l][j];
      for (i = 0; i < ImageWidth; i++) 
         for (j = s; j < ImageHeight-s; j++)
            for (l = -s+1; l <= s; l++)
               x[i][j] = (w[i][j-s] > w[i][j+l]) ? w[i][j-s] : w[i][j+l];
      for (i = 0; i < ImageWidth; i++) 
         for (j = s; j < ImageHeight-s; j++)
            for (l = -s+1; l <= s; l++)
               w[i][j] = (x[i][j-s] < x[i][j+l]) ? x[i][j-s] : x[i][j+l];

      for (j = 0; j < ImageHeight; j++)
         for (i = s; i < ImageWidth-s; i++) 
            for (l = -s+1; l <= s; l++)
               x[i][j] = (w[i-s][j] < w[i+l][j]) ? w[i-s][j] : w[i+l][j];
      for (j = 0; j < ImageHeight; j++)
         for (i = s; i < ImageWidth-s; i++)
            for (l = -s+1; l <= s; l++)
               w[i][j] = (x[i-s][j] > x[i+l][j]) ? x[i-s][j] : x[i+l][j];
      for (i = 0; i < ImageWidth; i++) 
         for (j = s; j < ImageHeight-s; j++)
            for (l = -s+1; l <= s; l++)
               x[i][j] = (w[i][j-s] < w[i][j+l]) ? w[i][j-s] : w[i][j+l];
      for (i = 0; i < ImageWidth; i++) 
         for (j = s; j < ImageHeight-s; j++)
            for (l = -s+1; l <= s; l++)
               w[i][j] = (x[i][j-s] > x[i][j+l]) ? x[i][j-s] : x[i][j+l];
   }

   for (i = 0; i < ImageWidth; i++){
      for (j = 0; j < ImageHeight; j++) {
         pixels[i][j] -= w[i][j]; /* y now contains the background */
         if (pixels[i][j] < 0) pixels[i][j] = 0; /* w now contains the background-subtracted image */
      }
   }
   
   for(i = 0; i < ImageWidth; i++){
         free(x[i]);
         free(w[i]);
   }
   free(x);
   free(w);
}

void asfFaster(int **pixels, int ImageWidth, int ImageHeight)
#define SE 4 /* max half-size of the square structuring element */
{
    int i, j, l, m, n, ii, jj, s;
    int **w;
	int *xstrip, *ystrip;	
     
//    x = malloc(sizeof(int *) * ImageWidth);
    w = malloc(sizeof(int *) * ImageWidth);

	if ((xstrip = (int *)malloc(sizeof(int)*ImageWidth) ) == NULL ){
        printf("\nError, memory not allocated.\n");
        exit(1);
    }
	if ((ystrip = (int *)malloc(sizeof(int)*ImageHeight) ) == NULL ){
        printf("\nError, memory not allocated.\n");
        exit(1);
    }     
    
	for(i = 0;i < ImageWidth; i++){
//          x[i] = malloc(sizeof(int) * ImageHeight);
          w[i] = malloc(sizeof(int) * ImageHeight);
          }
    
       /* now filter pepper noise */
   Rprintf("Filtering pepper noise...\n");
   for (j = 1; j < ImageHeight-1; j++)
      for (i = 1; i < ImageWidth-1; i++) {
         m = pixels[i][j];
         n = pixels[i][j-1];
         pixels[i][j] = n;
         for (ii = i-1; ii <= i+1; ii++)
            for (jj = j-1; jj <= j+1; jj++)
               if (pixels[ii][jj] < n) n = pixels[ii][jj]; 
         pixels[i][j] = (m > n) ? m : n; /* max(m, n) */
      }

   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth; i++) {
         w[i][j] = pixels[i][j];
//         x[i][j] = w[i][j];
      }

   /* now do the first step with a 2x2 structuring element */
   Rprintf("Opening and closing with a 2x2 element...\n");

   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth-1; i++) 
         w[i][j] = (w[i][j] > w[i+1][j]) ? w[i][j] : w[i+1][j];
   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth-1; i++) 
         w[i][j] = (w[i][j] < w[i+1][j]) ? w[i][j] : w[i+1][j];
   for (i = 0; i < ImageWidth; i++) 
      for (j = 0; j < ImageHeight-1; j++)
         w[i][j] = (w[i][j] > w[i][j+1]) ? w[i][j] : w[i][j+1];
   for (i = 0; i < ImageWidth; i++) 
      for (j = 0; j < ImageHeight-1; j++)
         w[i][j] = (w[i][j] < w[i][j+1]) ? w[i][j] : w[i][j+1];

   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth-1; i++) 
         w[i][j] = (w[i][j] < w[i+1][j]) ? w[i][j] : w[i+1][j];
   for (j = 0; j < ImageHeight; j++)
      for (i = 0; i < ImageWidth-1; i++) 
         w[i][j] = (w[i][j] > w[i+1][j]) ? w[i][j] : w[i+1][j];
   for (i = 0; i < ImageWidth; i++) 
      for (j = 0; j < ImageHeight-1; j++)
         w[i][j] = (w[i][j] < w[i][j+1]) ? w[i][j] : w[i][j+1];
   for (i = 0; i < ImageWidth; i++) 
      for (j = 0; j < ImageHeight-1; j++)
         w[i][j] = (w[i][j] > w[i][j+1]) ? w[i][j] : w[i][j+1];

   /* now do the alternating sequential filtering */
//   for (s = 1; s <= SE; s++) {
	 s = SE;
      Rprintf("Calculating ASF with s = %i...\n", s);

      for (j = 0; j < ImageHeight; j++){
         for (i = s; i < ImageWidth-s; i++) 
            for (l = -s+1; l <= s; l++)
				 xstrip[i] = (w[i-s][j] > w[i+l][j]) ? w[i-s][j] : w[i+l][j];
			for (i = s; i < ImageWidth-s; i++)
         		w[i][j] = xstrip[i];
		}

      for (j = 0; j < ImageHeight; j++){
         for (i = s; i < ImageWidth-s; i++)
            for (l = -s+1; l <= s; l++)
				xstrip[i] = (w[i-s][j] < w[i+l][j]) ? w[i-s][j] : w[i+l][j];
			for (i = s; i < ImageWidth-s; i++)
         		w[i][j] = xstrip[i];
		}

      for (i = 0; i < ImageWidth; i++){
         for (j = s; j < ImageHeight-s; j++)
            for (l = -s+1; l <= s; l++)
				ystrip[j] = (w[i][j-s] > w[i][j+l]) ? w[i][j-s] : w[i][j+l];
			for (j = s; j < ImageHeight-s; j++)
			    w[i][j] = ystrip[j];
		}


      for (i = 0; i < ImageWidth; i++) {
         for (j = s; j < ImageHeight-s; j++)
            for (l = -s+1; l <= s; l++)
               ystrip[j] = (w[i][j-s] < w[i][j+l]) ? w[i][j-s] : w[i][j+l];
			for (j = s; j < ImageHeight-s; j++)
			   w[i][j] = ystrip[j];
		}

      for (j = 0; j < ImageHeight; j++){
         for (i = s; i < ImageWidth-s; i++) 
            for (l = -s+1; l <= s; l++)
               xstrip[i] = (w[i-s][j] < w[i+l][j]) ? w[i-s][j] : w[i+l][j];
      		for (i = s; i < ImageWidth-s; i++)
         		w[i][j] = xstrip[i];
   		}

      for (j = 0; j < ImageHeight; j++){
         for (i = s; i < ImageWidth-s; i++)
            for (l = -s+1; l <= s; l++)
               xstrip[i] = (w[i-s][j] > w[i+l][j]) ? w[i-s][j] : w[i+l][j];
      		for (i = s; i < ImageWidth-s; i++)
         		w[i][j] = xstrip[i];
   		}

      for (i = 0; i < ImageWidth; i++) {
         for (j = s; j < ImageHeight-s; j++)
            for (l = -s+1; l <= s; l++)
				ystrip[j] = (w[i][j-s] < w[i][j+l]) ? w[i][j-s] : w[i][j+l];
			for (j = s; j < ImageHeight-s; j++)
         		w[i][j] = ystrip[j];
		}

      for (i = 0; i < ImageWidth; i++) {
         for (j = s; j < ImageHeight-s; j++)
            for (l = -s+1; l <= s; l++)
				ystrip[j] = (w[i][j-s] > w[i][j+l]) ? w[i][j-s] : w[i][j+l];
      		for (j = s; j < ImageHeight-s; j++)
         		w[i][j] = ystrip[j];
   		}
   //}

   for (i = 0; i < ImageWidth; i++){
      for (j = 0; j < ImageHeight; j++) {
         pixels[i][j] -= w[i][j]; /* w now contains the background */
         if (pixels[i][j] < 0) pixels[i][j] = 0; /* pixels now contains the background-subtracted image */
      }
   }
   
   for(i = 0; i < ImageWidth; i++){
//         free(x[i]);
         free(w[i]);
   }
//   free(x);
	 free(xstrip);
	 free(ystrip);
   free(w);
}

void getPixelIntensities(int **pixels, FILE *fp, int ImageWidth, int ImageHeight, int ImageSize, int BeginData){

      char *c;
      char r[1];
      int l, i, j;    
           
      c = malloc(ImageSize);      
           
      for (i = 1; i <= BeginData; i++) {
          fread(r, sizeof(r), 1, fp);
          }
    
      fread(c, 1, ImageSize, fp);

      /* convert tif binary data into image intensities */
      l = 0;
      for (j = 0; j < ImageHeight; j++){
          for (i = 0; i < ImageWidth; i++) {
              pixels[i][j] = number(2, c+l); 
              l += 2;
              }
      }
      free(c);
}

double matrixMean(int **matrix, int xStart, int yStart){
       int i, j;
       double result = 0;
     
       for(i = 0; i < 3; i++){
             for(j = 0; j < 3; j++){
                   result += matrix[xStart+i][yStart+j];
             }
       }

       return((result/9));
}

void calculateBackground(int **pixels, double *xs, double *ys, int numBeads, int ImageWidth, int ImageHeight, double *background, int n){

     double dist[4];
     double xc, yc, opt1, opt2, opt3, opt4;
     int i, j, k, count, newcoord[2], cX[4], cY[4], M[289];
     int n2 = floor(n/2);
     int temp = 0;
     
     for(i = 0; i < numBeads; i++){
           xc = xs[i] - floor(xs[i]);
           yc = ys[i] - floor(ys[i]);

		   opt1 = pow(xc, 2.0);
		   opt2 = pow(yc, 2.0);
		   opt3 = pow((xc - 1), 2.0);
		   opt4 = pow((yc - 1), 2.0);

		   dist[0] = opt1 + opt2;
		   dist[1] = opt3 + opt2;
		   dist[2] = opt1 + opt4;
		   dist[3] = opt3 + opt4;
/*           
           dist[0] = pow(xc, 2.0) + pow(yc, 2.0);
           dist[1] = pow((xc - 1), 2.0) + pow(yc, 2.0);
           dist[2] = pow(xc, 2.0) + pow((yc - 1), 2.0);
           dist[3] = pow((xc - 1), 2.0) + pow((yc - 1), 2.0);
*/           
           cX[0] = cX[2] = cY[0] = cY[1] = 0;
           cX[1] = cX[3] = cY[2] = cY[3] = 1;
           
           /* find the closest pixel to given co-ords */
           for(j = 0; j < 4; j++){
                 if(dist[j] < dist[temp]){
                            temp = j;
                 }
           }
           
           newcoord[0] = floor(xs[i]) + cX[temp];
           newcoord[1] = floor(ys[i]) + cY[temp];
     
            if (((newcoord[0] - n2) < 0) || ((newcoord[0] + n2) > ImageWidth) || 
            ((newcoord[1] - n2) < 0) || ((newcoord[1] + n2) > ImageHeight)) {
            /* crap average to fill in a blank if a bead is too close to the edge of the image */
  //          Rprintf("crap background\n");
                         background[i] = 700;
            }
            else {
                  count = 0;
                  for(j = 0; j < n; j++){
                        for(k = 0; k < n; k++){
                              M[count+k] = pixels[(newcoord[0] - n2) + j][(newcoord[1] - n2) + k];
                        }
                        count = count+n;
                  }
                  
           // qsort(M, 289, sizeof(int), (void *)comp_nums);
		   quicksort(M, 0, 289);

            background[i] = (M[0] + M[1] + M[2] + M[3] + M[4])/5;
            }
     }
} 

void HIPForeground(int **pixels, double *xs, double *ys, int numBeads, int ImageWidth, int ImageHeight, double *foreground){
      
      int x2, y2, i, j, k, count;
      int M[25];
      
      for(i=0; i < numBeads; i++){       
          x2 = floor(xs[i]);
          y2 = floor(ys[i]);
               
          if((x2 < 3) || (x2 > (ImageWidth) - 3) || (y2 < 3) || (y2 > (ImageHeight - 3))){
              foreground[i] = 0;
              //Rprintf("Bead %d is too close to the edge of the image to be evaluated and has been ignored.\n", i);
          }
          else{
               count = 0;
               for(j = 0; j < 5; j++){
                        for(k = 0; k < 5; k++){
                              M[count+k] = pixels[(x2 - 2) + j][(y2 - 2) + k];
                        }
                        count = count+5;
                  }
//          qsort(M, 25, sizeof(int), (void *)comp_nums);
		  quicksort(M, 0, 25);          
          foreground[i] = M[24];
        }
    }
}
                 

void IlluminaForeground(int **pixels, double *xs, double *ys, int numBeads, int ImageWidth, int ImageHeight, double *foreground){

//     int x2[numBeads], y2[numBeads];
     int *x2, *y2;
     int i;
     double av[4], w[4];//, L[numBeads];
     double xc, yc;
     
     x2 = malloc(sizeof(int) * numBeads);
     y2 = malloc(sizeof(int) * numBeads);

     for(i = 0; i < numBeads; i++){
           x2[i] = floor(xs[i]);
           y2[i] = floor(ys[i]);
           }          

     for(i = 0; i < numBeads; i++){
           if((x2[i] < 3) || (x2[i] > (ImageWidth) - 3) || (y2[i] < 3) || (y2[i] > (ImageHeight - 3))){
                 foreground[i] = 0;
                 //Rprintf("Bead %d is too close to the edge of the image to be evaluated and has been ignored.\n", i);
           }
           else {
                av[0] = matrixMean(pixels, (x2[i] - 1), (y2[i] - 1));
                av[1] = matrixMean(pixels, (x2[i] - 1), y2[i]);
                av[2] = matrixMean(pixels, x2[i], y2[i]);
                av[3] = matrixMean(pixels, x2[i], (y2[i] - 1));

                xc = xs[i] - floor(xs[i]);
                yc = ys[i] - floor(ys[i]);

                w[0] = ((1 - xc) * (1 - yc));
                w[1] = ((1 - xc) * yc);
                w[2] = (xc * yc); 
                w[3] = (xc * (1 - yc));
                foreground[i] = ((w[0] * av[0]) + (w[1] * av[1]) + (w[2] * av[2]) + (w[3] * av[3]));
           }
     }
     free(x2);
     free(y2);
}

void startEndPos(int *ProbeIDs, int numBeads, int *starts, int *ends){
	int i, j;

	j = 0;
	for(i = 1; i < numBeads; i++){
		  starts[0] = 1;
		  if(ProbeIDs[i] != ProbeIDs[i-1]){
		  	ends[j] = i;
			starts[j+1] = i+1;
			j++;
		}
	}
	ends[j] = i;
}


void readBeadImage(char **tif, double *xs, double *ys, int *numBeads, double *foreground, double *background, int *n, int *manip, int *fground){
     
      FILE *fp;
      int i, BeginData, len, ImageWidth, ImageHeight, ImageSize, StripOffset;
      long tag, type, length, offset, records;
      char c2[2], c4[4];    
      int **pixels;

      if ( (fp = fopen(tif[0], "rb")) == NULL)
      {
         Rprintf("Error opening file %s", tif[0]);
         exit(0);
      }

      if (fread(c2, 1, 2, fp) == 2) {
 //        Rprintf("Mode: %u %u\n", c2[0], c2[1]);
         }
         
      if (fread(c2, 1, 2, fp) == 2) {
//         Rprintf("Version: %i\n", number(2, c2));
         }
         
      if (fread(c4, 1, 4, fp) == 4) {
         offset = number(4, c4);
//         Rprintf("Offset: %i\n", offset);
         }
         
      if (fseek(fp, offset, SEEK_SET) != 0) {
         Rprintf("Error in fseek()\n");
         exit(0);
         }
         
   if (fread(c2, 1, 2, fp) == 2) {
      records = number(2, c2);
 //     Rprintf("Records: %i\n", records);
   }
   for (i = 1; i <= records; i++) {
      fread(c2, 1, 2, fp);
      tag = number(2, c2);
      fread(c2, 1, 2, fp);
      type = number(2, c2);
      fread(c4, 1, 4, fp);
      length = number(4, c4);
      fread(c4, 1, 4, fp);
      offset = number(4, c4);
  //    Rprintf("%i %i %i %i\n", tag, type, length, offset);
      switch (tag) {
         case (256) :
            ImageWidth = offset;
         case (257) :
            ImageHeight = offset;
         case (273) : 
            StripOffset = offset;
      }
   }
   if (fseek(fp, StripOffset, SEEK_SET) != 0) {
      Rprintf("Error in fseek()\n");
      exit(0);
   }
   if (fread(c4, 1, 4, fp) == 4) {
      BeginData = number(4, c4);
   }
   ImageSize = ImageWidth*ImageHeight*2;
   if (fseek(fp, 0, SEEK_END) != 0) {
      Rprintf("Error in fseek()\n");
      exit(0);
   }
   len = ftell(fp);
   if ((BeginData) + (ImageSize) > len) {
      BeginData = StripOffset; /* it's not a tif */
   }

   rewind(fp);
   
   pixels = malloc(sizeof(int *) * ImageWidth);

      for(i = 0;i < ImageWidth; i++){
           pixels[i] = malloc(sizeof(int) * ImageHeight);
           }
           
   Rprintf("Reading pixels of %s\n", tif[0]); 
   getPixelIntensities(pixels, fp, ImageWidth, ImageHeight, ImageSize, BeginData);
/*
   if((*manip) != 2){
   			   Rprintf("Calculating background\n");
               calculateBackground(pixels, xs, ys, *numBeads, ImageWidth, ImageHeight, background, *n);
               }
*/
   switch(*manip){
                  case 1: 
				  	   Rprintf("Calculating background\n");
               		   calculateBackground(pixels, xs, ys, *numBeads, ImageWidth, ImageHeight, background, *n);
               		   Rprintf("Sharpening Image\n");
                       sharpen(pixels, ImageWidth, ImageHeight);
                       break;
                  case 2:
                       Rprintf("Morphological Background\n");
                       asf(pixels, ImageWidth, ImageHeight);
                       break;
				  case 3:
				  	   Rprintf("ASF Faster\n");
					   asfFaster(pixels, ImageWidth, ImageHeight);
					   break;	   				
                  default:
				  	   Rprintf("Calculating background\n");
               		   calculateBackground(pixels, xs, ys, *numBeads, ImageWidth, ImageHeight, background, *n);	  
                       break;
                  }

   Rprintf("Calculating foregound\n");  
   switch(*fground){
                    case 0:
						IlluminaForeground(pixels, xs, ys, *numBeads, ImageWidth, ImageHeight, foreground);
						break;
					case 1:
						HIPForeground(pixels, xs, ys, *numBeads, ImageWidth, ImageHeight, foreground);
						break;
					default:
						break;
					}  

//   startEndPos(ProbeIDs, *numBeads, starts, ends);
   
   for(i = 0; i < ImageWidth; i++){
         free(pixels[i]);
   }
   free(pixels);
   
}
