#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <errno.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <time.h>
#include "implementierung.h"
void grayscale_naive(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    //initialise the R G and B values
    float* result = (float*) tmp;
    uint8_t R;
    uint8_t G;
    uint8_t B;
    for (size_t i=0;i < height * width * 3;i+=3) {
        //loading the colours' variables 
        R= img[i];
        G= img[i+1];
        B= img[i+2];
        //calculate the grayscale D value round it and store it to its corresponding place
        result[i/3] = round(a * R + b * G + c * B);
    }
}

// lookup grayscale
// Function to generate a precomputed table for a given coefficient
void gentable(float coeff, float* table) {
    for (size_t i = 0; i < 256; ++i) {
        table[i] = i * coeff;
    }
}

void grayscale_lookup(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    // Generate tables for coefficients
    float tableA[256];
    float tableB[256];
    float tableC[256];
    gentable(a, tableA);
    gentable(b, tableB);
    gentable(c, tableC);
    float* result =(float*)(tmp);
    //initialise the colour values which are also indexes to the prelookup tables
     uint8_t R;
     uint8_t G;
     uint8_t B;
    for (size_t i = 0; i < height * width*3; i+=3) {
        //load R G and B 
         R = img[i];
         G = img[i+ 1];
         B = img[i+ 2];
        // Lookup precomputed values from tables
        float D = round(tableA[R] + tableB[G] + tableC[B]);       
        result[i/3] = D;
    }
}

// simd grayscale
void grayscale_simd(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    float* result =(float*)(tmp);
    // Prepare SIMD constants for coefficients
    __m128 coefA = _mm_set1_ps(a);
    __m128 coefB = _mm_set1_ps(b);
    __m128 coefC = _mm_set1_ps(c);
    //initialisation to the colours' vectors
    __m128 pixelR;
    __m128 pixelG;
    __m128 pixelB;
    //initialisation to the variables of the results of the multiplication
    __m128 Ra;
    __m128 Gb;
    __m128 Bc;
    size_t i;
    // Process four pixels at a time
    for (i = 0; i < (height * width * 3)-((height*width*3 )%12); i += 12) {
        // Load pixel values as floats using SIMD instruction _mm_set_ps for R, G, B channels
         pixelR = _mm_set_ps((float)img[i+9], (float)img[i+6], (float)img[i+3], (float)img[i]);
         pixelG = _mm_set_ps((float)img[i+10], (float)img[i+7], (float)img[i+4], (float)img[i+1]);
         pixelB = _mm_set_ps((float)img[i+11], (float)img[i+8], (float)img[i+5], (float)img[i+2]);

        // Multiply each color channel with the corresponding coefficient
         Ra = _mm_mul_ps(pixelR, coefA);
         Gb = _mm_mul_ps(pixelG, coefB);
         Bc = _mm_mul_ps(pixelB, coefC);

        // Sum the results to get grayscale values
__m128 D = _mm_add_ps(_mm_add_ps(Ra, Gb), Bc);
__m128 rounding =_mm_set1_ps(0.5f);
  __m128 D1 = _mm_add_ps(D, rounding);
        // Round the floating-point values to the nearest integer
__m128 D_rounded = _mm_floor_ps(D1);
_mm_storeu_ps(result + i/3, D_rounded);
     }

    // Handle any remaining pixels
    for (  ; i < height * width * 3; i += 3) {
     
        result[i/3] =round(  a *img[i] + b * img[i+1] + c * img[i+2]);
    }
}


// naive interpolation calculation 
void interpolation_calculation_naive(size_t width, size_t height, size_t factor, float* inputArray, uint8_t* result) {
    // Calculate scaled width and height
    size_t scaledW = width * factor;  
    size_t scaledH = height * factor; 

    // Inverse factor for interpolation
    float invFactor = 1.0f / factor;

    // Initialize variables outside the loops
    float gy, oneMinusGy, gx, oneMinusGx;
    size_t yIdx, yNextIdx, xIdx, xNextIdx;
    float Q00, Qs0, Q0s, Qss;

    // Iterate over each pixel in the scaled image
    for (size_t yScaled = 0; yScaled < scaledH; yScaled++) {
        //precalculate  intermediate values
        gy = (float)(yScaled % factor) * invFactor;
        oneMinusGy = 1.0f - gy;
        //calculate indexes for the given 
        yIdx = yScaled / factor;
        yNextIdx = (yIdx + 1) % height;

        for (size_t xScaled = 0; xScaled < scaledW; xScaled++) {
            gx = (float)(xScaled % factor) * invFactor;
            oneMinusGx = 1.0f - gx;
            xIdx = xScaled / factor;
            xNextIdx = (xIdx + 1) % width;

            // Retrieve pixel values from the four surrounding corners
            Q00 = inputArray[yIdx * width + xIdx];
            Qs0 = inputArray[yIdx * width + xNextIdx];
            Q0s = inputArray[yNextIdx * width + xIdx];
            Qss = inputArray[yNextIdx * width + xNextIdx];

            // Interpolate and store the result
            result[yScaled * scaledW + xScaled] = (uint8_t)((oneMinusGy * (oneMinusGx * Q00 + gx * Qs0) + gy * (oneMinusGx * Q0s + gx * Qss)) + 0.5f);
        }
    }
}

// algorithmically optimized interpolation calculation
void interpolation_calculation_algorithm_optimized(size_t width, size_t height, size_t factor,float* inputArray, uint8_t* result,size_t startWidth, size_t startHeight) {
   float Q00, Qs0, Q0s, Qss; // Variables to store pixel values for interpolation
   float val1, val2; // Intermediate values for bilinear interpolation

   // Iterate over height
   for (size_t h = startHeight; h < height; h++) {
       // Iterate over width
       for (size_t b = startWidth; b < width; b++) {
           // Get the pixel values for the current and adjacent pixels
           Q00 = inputArray[b + h * width]; // Current pixel 
           Qs0 = inputArray[(b + 1) % width + h * width]; // Right neighbor
           Q0s = inputArray[b + ((h + 1) % height) * width]; // Bottom neighbor
           Qss = inputArray[(b + 1) % width + ((h + 1) % height) * width]; // Diagonal neighbor

           // Perform bilinear interpolation for each subpixel within the current sector
           for (size_t f1 = 0; f1 < factor; f1++) {
               // Calculate intermediate values
               val1 = ((float)(factor - f1) / factor) * Q00 + ((float)f1 / factor) * Qs0;
               val2 = ((float)(factor - f1) / factor) * Q0s + ((float)f1 / factor) * Qss;

               // Calculate the final interpolated values and assign them to the result
               for (size_t f = 0; f < factor; f++) {
                   result[(h * factor + f) * width * factor + (b * factor + f1)] =
                       ((float)(factor - f) / factor) * val1 +
                       ((float)f / factor) * val2 + 0.5f; // Adding 0.5f for rounding
               }
           }
       }
   }
}

// simd interpolation calculation
void interpolation_calculation_simd( size_t width, size_t height, size_t factor,float* inputArray, uint8_t* result) {
    size_t scaledWidth = width * factor;//calculate new width
    float invFactor = 1.0f / factor;  // Inverse of scaling factor for coefficient calculation

    // Allocate memory for precomputed interpolation coefficients a b c d
    float* coeffTable = (float*)malloc(sizeof(float) * factor * factor * 4); 

    float gx,gy;//initialise intermediate values
    size_t index,x,y,sectorBaseX,coeffIndex,resultIndex;//initlialise indexes
    __m128 Q00,Qs0,Q0s,Qss,a,b,c,d,interpolated,coeff;//initialise intermediate values and final interpolated value
    for (size_t i = 0; i < factor; i++) {
        for (size_t j = 0; j < factor; j++) {
             gx = j * invFactor;
             gy = i * invFactor;
             index = (i * factor + j) * 4;
            // Precompute bilinear interpolation coefficients for each sub-pixel position
            coeffTable[index] = (1 - gx) * (1 - gy);     // Coefficient a
            coeffTable[index + 1] = gx * (1 - gy);       // Coefficient b
            coeffTable[index + 2] = (1 - gx) * gy;       // Coefficient c
            coeffTable[index + 3] = gx * gy;             // Coefficient d
        }
    }

    //  SIMD processing loop over the image
    for (y = 0; y < height - 1; y++) {
        for ( x = 0; x < (width - 1)-((width-1)%4); x += 4) { 
            // Load pixel values for Q00, Q0s, Qs0, Qss using SIMD
             Q00 = _mm_loadu_ps(&inputArray[y * width + x]);
             Q0s = _mm_loadu_ps(&inputArray[(y + 1) * width + x]);
             Qs0 = _mm_loadu_ps(&inputArray[y * width + x + 1]);
             Qss = _mm_loadu_ps(&inputArray[(y + 1) * width + x + 1]);

            // Iterate over each sub-pixel position within the pixel
            for (size_t i = 0; i < factor; i++) {
                for (size_t j = 0; j < factor; j++) {
                     coeffIndex = (i * factor + j) * 4;
                     coeff = _mm_loadu_ps(&coeffTable[coeffIndex]);

                    // Load the precomputed coefficients into SIMD 128bit registers
                     a = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(0, 0, 0, 0));
                     b = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(1, 1, 1, 1));
                     c = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(2, 2, 2, 2));
                     d = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(3, 3, 3, 3));

                    // Perform bilinear interpolation using SIMD
                     interpolated = _mm_add_ps(_mm_add_ps(_mm_mul_ps(a, Q00), _mm_mul_ps(b, Qs0)), _mm_add_ps(_mm_mul_ps(c, Q0s), _mm_mul_ps(d, Qss)));

                    // Store the interpolated results into the result buffer
                    for (size_t k = 0; k < 4; k++) {
                         sectorBaseX = (x + k) * factor;
                         resultIndex = (y * factor + i) * scaledWidth + sectorBaseX + j;
                        result[resultIndex] = (interpolated[k] + 0.5f); // Convert to uint8_t with rounding
                    }
                }
            }
        }
    }

    // Process the edge cases by calling 
    interpolation_calculation_algorithm_optimized( width, height, factor,inputArray, result,(width - 1)-((width-1)%4), 0); // Rightmost column
    interpolation_calculation_algorithm_optimized( width, height, factor,inputArray, result, 0, height-1); // Bottom row

    free(coeffTable); // Free the allocated memory for coefficients
}


// handle scale factor 1
void scale1_handling(uint8_t* tmp, uint8_t* result, size_t width, size_t height) {
    for (size_t i = 0; i < width * height; i++) {
             result[i] = tmp[i];
    }
}

// naive interpolate
void interpolate_naive(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    grayscale_naive(img, width, height, a, b, c, tmp);
     if (scale_factor == 1) {
        scale1_handling(tmp, result, width, height);
        return;
    }
    interpolation_calculation_naive(width, height, scale_factor,(float*) tmp, result);
}

// algorithmically optimized interpolate
void interpolate_algorithm_optimized(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    grayscale_lookup(img, width, height, a, b, c, tmp);
    if (scale_factor == 1) {
        scale1_handling( tmp, result, width, height);
        return;
    }
    interpolation_calculation_algorithm_optimized(width, height, scale_factor,(float*) tmp, result,0,0);
}

// simd interpolate
void interpolate_simd(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
  //  grayscale_simd(img, width, height, a, b, c, tmp);
  grayscale_simd(img, width, height, a, b, c, tmp);
    if (scale_factor == 1) {
        scale1_handling( tmp, result, width, height);
        return;
    }
    interpolation_calculation_simd(width, height, scale_factor,(float*) tmp, result);
}

// interpolate (standardversion) (the best interpolation is interpolate_simd)
void interpolate(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    interpolate_simd(img, width, height, a, b, c, scale_factor, tmp, result);
}