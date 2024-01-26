#include <stdbool.h>
#include "stdio.h"
#include "inttypes.h"
#include <emmintrin.h> 
#include <stdint.h>
#include <time.h>
#include <math.h>

//algorithmically optimised interpolation
void interpolationopt(float* inputArray, uint8_t* result, size_t width, size_t height, size_t factor) {
   float v1_1, v1_2, v2_1, v2_2; // Variables to store pixel values for interpolation
   float val1, val2; // Intermediate values for bilinear interpolation

   // Iterate over height
   for (size_t h = 0; h < height; h++) {
       // Iterate over width
       for (size_t b = 0; b < width; b++) {
           // Get the pixel values for the current and adjacent pixels
           v1_1 = inputArray[b + h * width]; // Current pixel 
           v1_2 = inputArray[(b + 1) % width + h * width]; // Right neighbor
           v2_1 = inputArray[b + ((h + 1) % height) * width]; // Bottom neighbor
           v2_2 = inputArray[(b + 1) % width + ((h + 1) % height) * width]; // Diagonal neighbor

           // Perform bilinear interpolation for each subpixel within the current sector
           for (size_t f1 = 0; f1 < factor; f1++) {
               // Calculate intermediate values
               val1 = ((float)(factor - f1) / factor) * v1_1 + ((float)f1 / factor) * v1_2;
               val2 = ((float)(factor - f1) / factor) * v2_1 + ((float)f1 / factor) * v2_2;

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

//simd interpolation starts here
void processEdge(const float* inputArray, uint8_t* result, size_t width, size_t height, size_t factor, size_t startWidth, size_t startHeight) {
   float v1_1, v1_2, v2_1, v2_2; // Variables to store pixel values for interpolation
   float val1, val2; // Intermediate values for bilinear interpolation
   //iterate over height
    for (size_t h = startHeight; h < height; h++) {
        //iterate over width
        for (size_t b = startWidth; b < width; b++) {
            // Get the pixel values for the current and adjacent pixels
             v1_1 = inputArray[b + h * width]; //current pixel
             v1_2 = inputArray[(b + 1) % width + h * width]; //right neighbor
             v2_1 = inputArray[b + ((h + 1) % height) * width]; // Bottom neighbor
             v2_2 = inputArray[(b + 1) % width + ((h + 1) % height) * width]; // Diagonal neighbor

            // Perform bilinear interpolation for each subpixel within the current sector
            for (size_t f1 = 0; f1 < factor; f1++) {
                // Calculate intermediate values
                 val1 = ((float)(factor - f1) / factor) * v1_1 + ((float)f1 / factor) * v1_2;
                 val2 = ((float)(factor - f1) / factor) * v2_1 + ((float)f1 / factor) * v2_2;
                
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

void simdInterpolate(const float* inputArray, uint8_t* result, size_t width, size_t height, size_t factor) {
    size_t scaledWidth = width * factor;
    float invFactor = 1.0f / factor;  // Inverse of scaling factor for coefficient calculation

    // Allocate memory for precomputed interpolation coefficients a b c d
    float* coeffTable = (float*)malloc(sizeof(float) * factor * factor * 4); 
    for (size_t i = 0; i < factor; i++) {
        for (size_t j = 0; j < factor; j++) {
            float gx = j * invFactor;
            float gy = i * invFactor;
            size_t index = (i * factor + j) * 4;
            // Precompute bilinear interpolation coefficients for each sub-pixel position
            coeffTable[index] = (1 - gx) * (1 - gy);     // Coefficient a
            coeffTable[index + 1] = gx * (1 - gy);       // Coefficient b
            coeffTable[index + 2] = (1 - gx) * gy;       // Coefficient c
            coeffTable[index + 3] = gx * gy;             // Coefficient d
        }
    }

    //  SIMD processing loop over the image
    for (size_t y = 0; y < height - 1; y++) {
        for (size_t x = 0; x < width - 1; x += 4) { 
            // Load pixel values for Q11, Q12, Q21, Q22 using SIMD
            __m128 Q11 = _mm_loadu_ps(&inputArray[y * width + x]);
            __m128 Q21 = _mm_loadu_ps(&inputArray[(y + 1) * width + x]);
            __m128 Q12 = _mm_loadu_ps(&inputArray[y * width + x + 1]);
            __m128 Q22 = _mm_loadu_ps(&inputArray[(y + 1) * width + x + 1]);

            // Iterate over each sub-pixel position within the pixel
            for (size_t i = 0; i < factor; i++) {
                for (size_t j = 0; j < factor; j++) {
                    size_t coeffIndex = (i * factor + j) * 4;
                    __m128 coeff = _mm_loadu_ps(&coeffTable[coeffIndex]);

                    // Load the precomputed coefficients into SIMD 128bit registers
                    __m128 a = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(0, 0, 0, 0));
                    __m128 b = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(1, 1, 1, 1));
                    __m128 c = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(2, 2, 2, 2));
                    __m128 d = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(3, 3, 3, 3));

                    // Perform bilinear interpolation using SIMD
                    __m128 interpolated = _mm_add_ps(_mm_add_ps(_mm_mul_ps(a, Q11), _mm_mul_ps(b, Q12)), _mm_add_ps(_mm_mul_ps(c, Q21), _mm_mul_ps(d, Q22)));

                    // Store the interpolated results into the result buffer
                    for (size_t k = 0; k < 4; k++) {
                        size_t sectorBaseX = (x + k) * factor;
                        size_t resultIndex = (y * factor + i) * scaledWidth + sectorBaseX + j;
                        result[resultIndex] = (uint8_t)(interpolated[k] + 0.5f); // Convert to uint8_t with rounding
                    }
                }
            }
        }
    }

    // Process the edge cases by calling 
    processEdge(inputArray, result, width, height, factor, width - 1, 0); // Rightmost column
    processEdge(inputArray, result, width, height, factor, 0, height - 1); // Bottom row

    free(coeffTable); // Free the allocated memory for coefficients
}

//new naiv interpolation
void interpolate100(const float* inputArray, uint8_t* result, size_t width, size_t height, size_t factor) {
    size_t scaledW = width * factor;  // Scaled width
    size_t scaledH = height * factor; // Scaled height
    float invFactor = 1.0f / factor; // Inverse of the scaling factor

    // Iterate over each pixel in the scaled image
    for (size_t yScaled = 0; yScaled < scaledH; yScaled++) {
        //redundant intermediate calculated values
        float gy = (float)(yScaled % factor) * invFactor;
        float oneMinusGy = 1.0f - gy;
        size_t yIdx = (yScaled / factor) % height;
        size_t yNextIdx = ((yScaled + factor) / factor) % height;

        for (size_t xScaled = 0; xScaled < scaledW; xScaled++) {
                    //intermediate calculated values

            float gx = (float)(xScaled % factor) * invFactor;
            float oneMinusGx = 1.0f - gx;
            size_t xIdx = (xScaled / factor) % width;
            size_t xNextIdx = ((xScaled + factor) / factor) % width;

            // Fetch the pixel values from the four surrounding corners
            uint8_t Q11 = inputArray[yIdx * width + xIdx];
            uint8_t Q21 = inputArray[yIdx * width + xNextIdx];
            uint8_t Q12 = inputArray[yNextIdx * width + xIdx];
            uint8_t Q22 = inputArray[yNextIdx * width + xNextIdx];

            // Calculate interpolation coefficients
            float a = oneMinusGx * oneMinusGy;
            float b = gx * oneMinusGy;
            float c = oneMinusGx * gy;
            float d = gx * gy;

            // Interpolate and store the result
            result[yScaled * scaledW + xScaled] =((a * Q11 + b * Q21 + c * Q12 + d * Q22)+0.5f);
        }
    }
}