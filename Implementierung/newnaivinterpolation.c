#include <stdbool.h>
#include "stdio.h"
#include "inttypes.h"
#include <emmintrin.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

void interpolate100(const uint8_t* inputArray, uint8_t* result, int width, int height, int factor) {
    int scaledW = width * factor;  // Scaled width
    int scaledH = height * factor; // Scaled height
    float invFactor = 1.0f / factor; // Inverse of the scaling factor

    // Iterate over each pixel in the scaled image
    for (int yScaled = 0; yScaled < scaledH; yScaled++) {
        float gy = (float)(yScaled % factor) * invFactor;
        float oneMinusGy = 1.0f - gy;
        int yIdx = (yScaled / factor) % height;
        int yNextIdx = ((yScaled + factor) / factor) % height;

        for (int xScaled = 0; xScaled < scaledW; xScaled++) {
            float gx = (float)(xScaled % factor) * invFactor;
            float oneMinusGx = 1.0f - gx;
            int xIdx = (xScaled / factor) % width;
            int xNextIdx = ((xScaled + factor) / factor) % width;

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
            result[yScaled * scaledW + xScaled] = (uint8_t)(a * Q11 + b * Q21 + c * Q12 + d * Q22);
        }
    }
}
