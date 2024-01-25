
#include <stdbool.h>
#include "stdio.h"
#include "inttypes.h"
#include <emmintrin.h> 
#include <stdint.h>
#include <time.h>
#include <math.h>

void grayscale(const uint8_t* img, int width, int height, float a, float b, float c, uint8_t* result) {
    // Prepare SIMD constants for coefficients
    __m128 coefA = _mm_set1_ps(a);
    __m128 coefB = _mm_set1_ps(b);
    __m128 coefC = _mm_set1_ps(c);
    // Process four pixels at a time
    for (size_t i = 0; i < (height * width * 3)-((height*width*3 )%12); i += 12) {
        // Load pixel values as floats using SIMD instruction _mm_set_ps for R, G, B channels
        __m128 pixelR = _mm_set_ps((float)img[i+9], (float)img[i+6], (float)img[i+3], (float)img[i]);
        __m128 pixelG = _mm_set_ps((float)img[i+10], (float)img[i+7], (float)img[i+4], (float)img[i+1]);
        __m128 pixelB = _mm_set_ps((float)img[i+11], (float)img[i+8], (float)img[i+5], (float)img[i+2]);

        // Multiply each color channel with the corresponding coefficient
        __m128 Ra = _mm_mul_ps(pixelR, coefA);
        __m128 Gb = _mm_mul_ps(pixelG, coefB);
        __m128 Bc = _mm_mul_ps(pixelB, coefC);

        // Sum the results to get grayscale values
        __m128 grayF = _mm_add_ps(_mm_add_ps(Ra, Gb), Bc);
     //   __m128 grayF  = _mm_add_ps(grayF1,roound );
        __m128i grayI = _mm_cvtps_epi32(grayF); // Convert to 32-bit integers
        __m128i grayPacked = _mm_packus_epi16(_mm_packs_epi32(grayI, grayI), _mm_setzero_si128()); // Pack to 16-bit and then to 8-bit
        _mm_storel_epi64((__m128i *)(result + i/3), grayPacked); // Store lower 64 bits (4 packed grayscale values)



     }

    // Handle any remaining pixels
    for (size_t i ; i < height * width * 3; i += 3) {
        // Non-SIMD processing for the last few pixels
        float pixelR = (float)img[i];
        float pixelG = (float)img[i + 1];
        float pixelB = (float)img[i + 2];

        // Perform the grayscale conversion
        float D = a * pixelR + b * pixelG + c * pixelB;
        result[i / 3] = (uint8_t)(D > 255 ? 255 : (D < 0 ? 0 : D));
    }
}