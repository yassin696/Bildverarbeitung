#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

// Function to generate a precomputed table for a given coefficient
void gentable(float coeff, int16_t* table) {
    for (int i = 0; i < 256; ++i) {
        table[i] = (int16_t)(i * coeff);
    }
}

void grayscale(const uint8_t* img, int width, int height, float a, float b, float c, uint8_t* result) {
    printf("grayscale with precomputed tables\n");

    // Generate tables for coefficients
    int16_t tableA[256];
    int16_t tableB[256];
    int16_t tableC[256];
    gentable(a, tableA);
    gentable(b, tableB);
    gentable(c, tableC);

    for (int i = 0; i < height * width; ++i) {
        uint8_t R = img[i * 3];
        uint8_t G = img[i * 3 + 1];
        uint8_t B = img[i * 3 + 2];

        // Lookup precomputed values from tables
        int16_t weightedSum = tableA[R] + tableB[G] + tableC[B];

        
        result[i] = (uint8_t)(weightedSum );
    }
}