#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <immintrin.h>
#include <emmintrin.h>

void grayscale(const uint8_t* img, int width, int height, float a, float b, float c, uint8_t* result) {
    // Prepare SIMD constants for coefficients
    __m128 coef = _mm_set_ps(0.0f, c, b, a);

    // Process pixels one at a time
    for (size_t i = 0; i < height * width * 3; i += 3) {
        // Load the pixel values as floats using SIMD
        __m128 pixelR = _mm_set_ss((float)img[i]);
        __m128 pixelG = _mm_set_ss((float)img[i + 1]);
        __m128 pixelB = _mm_set_ss((float)img[i + 2]);

        // Multiply each color channel with the corresponding coefficient
        __m128 resultFloat = _mm_add_ps(_mm_add_ps(_mm_mul_ps(pixelR, coef), _mm_mul_ps(pixelG, coef)), _mm_mul_ps(pixelB, coef));

        // Convert to grayscale value
        uint8_t gray = (uint8_t)_mm_cvtss_f32(resultFloat);

        // Store the grayscale value in the result array
        result[i / 3] = gray;
    }
}
