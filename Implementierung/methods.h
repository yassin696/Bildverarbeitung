#ifndef METHODS_H
#define METHODS_H

#include <stdint.h>
#include <stddef.h>

// Function declarations here

void grayscale_naive(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp);
void gentable(float coeff, float* table);
void grayscale_lookup(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp);
void grayscale_simd(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp);
void interpolation_calculation_naive(size_t width, size_t height, size_t factor, float* inputArray, uint8_t* result);
void interpolation_calculation_algorithm_optimized(size_t width, size_t height, size_t factor,float* inputArray, uint8_t* result,size_t startWidth, size_t startHeight);
void interpolation_calculation_simd(size_t width, size_t height, size_t factor,float* inputArray, uint8_t* result);
void scale1_handling(const float* tmp, uint8_t* result, size_t width, size_t height);
void interpolate_naive(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result);
void interpolate_algorithm_optimized(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result);
void interpolate_simd(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result);

#endif // METHODS_H
