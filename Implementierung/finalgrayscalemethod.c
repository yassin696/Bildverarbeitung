#include <stdbool.h>
#include "stdio.h"
#include "inttypes.h"
#include <emmintrin.h> 
#include <stdint.h>
#include <time.h>
#include <math.h>

// Function to generate a precomputed table for a given coefficient
void gentable(float coeff, float* table) {
    for (size_t i = 0; i < 256; ++i) {
        table[i] = i * coeff;
    }
}

void grayscalePrecomputed(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* result) {
    printf("grayscale with precomputed tables\n");
float* result =(float*)(result);
    // Generate tables for coefficients
    float tableA[256];
    float tableB[256];
    float tableC[256];
    gentable(a, tableA);
    gentable(b, tableB);
    gentable(c, tableC);

//initialise the colour values which are also indexes to the prelookup tables
 uint8_t R;
 uint8_t G;
 uint8_t B;
    for (size_t i = 0; i < height * width*3; i+=3) {
        //load R G and B 
         R = img[i];
         G = img[i + 1];
         B = img[i + 2];
        // Lookup precomputed values from tables
        float D = (tableA[R] + tableB[G] + tableC[B]);       
        result[i/3] = D;
    }
}
//naiv grayscale implementation
void grayscale1(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    float* tmp =(float*)(tmp);
    //initialise the R G and B values
    uint8_t R;
    uint8_t G;
    uint8_t B;
    float D;
    for (size_t i=0;i < height * width * 3;i+=3) {
        //loading the colours' variables 
        R= img[i];
        G= img[i+1];
        B= img[i+2];
       
        //calculate the grayscale D value round it and store it to its corresponding place
        tmp[i/3] =  (a *R + b * G + c * B);
}
}
void grayscale(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* result) {
    float* result =(float*)(result);
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
        __m128 D1 = _mm_add_ps(_mm_add_ps(Ra, Gb), Bc);
      
        // Store the SIMD computed values directly into the result array
        _mm_storeu_ps(result + i/3, D1);
     }

    // Handle any remaining pixels
    for ( i ; i < height * width * 3; i += 3) {
     
        result[i/3] =  a *img[i] + b * img[i+1] + c * img[i+2];
    }
}