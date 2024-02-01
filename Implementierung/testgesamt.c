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
#include "methods.h"


double get_elapsed_time(struct timespec start, struct timespec end) {
    double start_sec = (double)start.tv_sec + (double)start.tv_nsec / 1e9; // Start time in seconds
    double end_sec = (double)end.tv_sec + (double)end.tv_nsec / 1e9;       // End time in seconds
    return end_sec - start_sec; // Return the difference in seconds
}
int main() {
    size_t pixel_counts[] = {100000,500000,1000000, 2000000, 4000000, 8000000, 16000000};
    size_t num_tests = sizeof(pixel_counts) / sizeof(pixel_counts[0]);
    size_t scale_factor = 2; // Assuming a scaling factor of 2
    float a = 0.299, b = 0.587, c = 0.114;

    for (size_t test = 0; test < num_tests; test++) {
        size_t num_pixels = pixel_counts[test];
        size_t width = sqrt(pixel_counts[test]);
        size_t height = num_pixels / width;

        printf("for a number of pixels :%li\n",num_pixels);
        // Allocate memory for random image data
        uint8_t* img = (uint8_t*)malloc(width * height * 3 * sizeof(uint8_t));
        if (!img) {
            fprintf(stderr, "Memory allocation for image failed\n");
            exit(EXIT_FAILURE);
        }

        // Fill the image with random data
        for (size_t i = 0; i < width * height * 3; i++) {
            img[i] = rand() % 256;
        }

        // Allocate memory for grayscale conversion and interpolation
        float* tmp = (float*)malloc(width * height * sizeof(float));
        if (!tmp) {
            fprintf(stderr, "Memory allocation for tmp failed\n");
            free(img);
            exit(EXIT_FAILURE);
        }

        // Allocate separate result buffers for each interpolation method
        uint8_t* result_naive = (uint8_t*)malloc(width * height * scale_factor * scale_factor * sizeof(uint8_t));
        uint8_t* result_optimized = (uint8_t*)malloc(width * height * scale_factor * scale_factor * sizeof(uint8_t));
        uint8_t* result_simd = (uint8_t*)malloc(width * height * scale_factor * scale_factor * sizeof(uint8_t));
        if (!result_naive || !result_optimized || !result_simd) {
            fprintf(stderr, "Memory allocation for result buffers failed\n");
            free(img);
            free(tmp);
            free(result_naive); // It's okay to call free on NULL
            free(result_optimized);
            free(result_simd);
            exit(EXIT_FAILURE);
        }

        // Measure and print the elapsed time for naive implementation
        printf("runtime testing : \n");
        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        interpolate_naive(img, width, height, a, b, c, scale_factor, (uint8_t*)tmp, result_naive);
        clock_gettime(CLOCK_MONOTONIC, &end);
        printf("Naive method runtime: %f seconds\n", get_elapsed_time(start, end));

        //Measure and print the elapsed time for algorithmically optimised implementation
      
        clock_gettime(CLOCK_MONOTONIC, &start);
        interpolate_algorithm_optimized(img, width, height, a, b, c, scale_factor, (uint8_t*)tmp, result_optimized);
        clock_gettime(CLOCK_MONOTONIC, &end);
        printf("algorithm optimised method runtime: %f seconds\n", get_elapsed_time(start, end));
        
        //Measure and print the elapsed time for SIMD optimised implementation
       
        clock_gettime(CLOCK_MONOTONIC, &start);
        interpolate_simd(img, width, height, a, b, c, scale_factor, (uint8_t*)tmp, result_simd);
        clock_gettime(CLOCK_MONOTONIC, &end);
        printf("simd method runtime: %f seconds\n", get_elapsed_time(start, end));
        // Compare the results from the three methods
        printf("\n");
        printf("    Accuracy testing part :\n");
        size_t diffs_naive_optimized = 0, diffs_naive_simd = 0, diffs_optimized_simd = 0;
        size_t total_pixels = width * height * scale_factor * scale_factor;
        for (size_t i = 0; i < total_pixels; i++) {
    if (result_naive[i] != result_optimized[i]) {
        if (result_naive[i] != result_optimized[i] ) {
            diffs_naive_optimized++;
        }
    }

    if (result_naive[i] != result_simd[i]) {
        if (result_naive[i] != result_simd[i] ) {
            diffs_naive_simd++;
        }
    }

    if (result_optimized[i] != result_simd[i]) {
        if (result_optimized[i] != result_simd[i] ) {
            diffs_optimized_simd++;
        }
    }
}

        // Print the comparison results
        printf("    Differences between naive and optimized: %zu\n", diffs_naive_optimized);
        printf("    Differences between naive and SIMD: %zu\n", diffs_naive_simd);
        printf("    Differences between optimized and SIMD: %zu\n", diffs_optimized_simd);
        printf("\n");
      
    // Free the allocated memory for result arrays
        free(img);
        free(tmp);
        free(result_naive);
        free(result_optimized);
        free(result_simd);  }  
        // Clean up
 size_t widths = 512, heights = 512;
//uint8_t* img = (uint8_t*)malloc(width * height * 3 * sizeof(uint8_t));



 //   float a = 0.299, b = 0.587, c = 0.114;
   

    // Allocate memory for random image data
    uint8_t* imgs = (uint8_t*)malloc(widths * heights * 3 * sizeof(uint8_t));
    if (!imgs) {
        fprintf(stderr, "Memory allocation for image failed\n");
        return 1;
    }

    // Fill the image with random data
    for (size_t i = 0; i < widths * heights * 3; i++) {
        imgs[i] = rand() % 256;
    }

    for (size_t scale_factor = 1; scale_factor <= 20; scale_factor++) {
        if (scale_factor==20)scale_factor=26;
        printf("\n");
        printf("Testing for scale factor %zu:\n", scale_factor);

        // Allocate memory for grayscale conversion
        float* tmps = (float*)malloc(widths * heights * sizeof(float));
        if (!tmps) {
            fprintf(stderr, "Memory allocation for tmp failed\n");
            free(imgs);
            return 1;
        }

        // Allocate memory for the result buffers
        size_t scaled_width = widths * scale_factor;
        size_t scaled_height = heights * scale_factor;
        uint8_t* result_naive = (uint8_t*)malloc(scaled_width * scaled_height * sizeof(uint8_t));
        uint8_t* result_optimized = (uint8_t*)malloc(scaled_width * scaled_height * sizeof(uint8_t));
        uint8_t* result_simd = (uint8_t*)malloc(scaled_width * scaled_height * sizeof(uint8_t));

        if (!result_naive || !result_optimized || !result_simd) {
            fprintf(stderr, "Memory allocation for result buffers failed\n");
            free(imgs);
            free(tmps);
            free(result_naive);
            free(result_optimized);
            free(result_simd);
            return 1;
        }

        struct timespec start, end;

        // Naive method
        clock_gettime(CLOCK_MONOTONIC, &start);
        interpolate_naive(imgs, widths, heights, a, b, c, scale_factor, (uint8_t*)tmps, result_naive);
        clock_gettime(CLOCK_MONOTONIC, &end);
        printf("Naive method runtime: %f seconds\n", get_elapsed_time(start, end));

        // Optimized method
        clock_gettime(CLOCK_MONOTONIC, &start);
        interpolate_algorithm_optimized(imgs, widths, heights, a, b, c, scale_factor, (uint8_t*)tmps, result_optimized);
        clock_gettime(CLOCK_MONOTONIC, &end);
        printf("Optimized method runtime: %f seconds\n", get_elapsed_time(start, end));

        // SIMD method
        clock_gettime(CLOCK_MONOTONIC, &start);
        interpolate_simd(imgs, widths, heights, a, b, c, scale_factor, (uint8_t*)tmps, result_simd);
        clock_gettime(CLOCK_MONOTONIC, &end);
        printf("SIMD method runtime: %f seconds\n", get_elapsed_time(start, end));

        // Free allocated memory
        free(tmps);
        free(result_naive);
        free(result_optimized);
        free(result_simd);
    }

    free(imgs);
        
    

    return 0;
}