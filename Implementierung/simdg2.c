#include <stdbool.h>
#include "stdio.h"
#include "inttypes.h"
#include <emmintrin.h> 
#include <stdint.h>
#include <time.h>
#include <math.h>

// Function to generate a precomputed table for a given coefficient
void gentable(float coeff, float* table) {
    for (int i = 0; i < 256; ++i) {
        table[i] = i * coeff;
    }
}

void grayscalePrecomputed(const uint8_t* img, int width, int height, float a, float b, float c, uint8_t* result) {
    printf("grayscale with precomputed tables\n");

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
    for (int i = 0; i < height * width*3; i+=3) {
        //load R G and B 
         R = img[i ];
         G = img[i + 1];
         B = img[i + 2];
        // Lookup precomputed values from tables
        float D = (tableA[R] + tableB[G] + tableC[B]+(float)(0.5));       
        result[i/3] = D;
    }
}
//naiv grayscale implementation
void grayscale1(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    //initialise the R G and B values
    uint8_t R;
    uint8_t G;
    uint8_t B;
    for (size_t i=0;i < height * width * 3;i+=3) {
        //loading the colours' variables 
        R= img[i];
        G= img[i+1];
        B= img[i+2];
        //calculate the grayscale D value round it and store it to its corresponding place
        tmp[i/3] =  ((a *R + b * G + c * B)+(float)0.5);
}
}
void grayscale(const uint8_t* img, int width, int height, float a, float b, float c, uint8_t* result) {
    // Prepare SIMD constants for coefficients
    __m128 coefA = _mm_set1_ps(a);
    __m128 coefB = _mm_set1_ps(b);
    __m128 coefC = _mm_set1_ps(c);
    //set a vector by 4 duplicated 0.5 that we will need to convert later
    __m128 conv  =_mm_set1_ps((float)0.5);
    //initialisation to the colours' vectors
    __m128 pixelR;
    __m128 pixelG;
    __m128 pixelB;
    //initialisation to the variables of the results of the multiplication
    __m128 Ra;
    __m128 Gb;
    __m128 Bc;
    // Process four pixels at a time
    for (size_t i = 0; i < (height * width * 3)-((height*width*3 )%12); i += 12) {
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
        // add 0.5 to each grayscale value so that the cast to uin8_t will behave like a round function
        __m128 D  = _mm_add_ps(D1,conv);

            result[i/3]=D[0];
            result[i/3 +1]=D[1];
            result[i/3 +2]=D[2];
            result[i/3 +3]=D[3];
     }

    // Handle any remaining pixels
    for (size_t i ; i < height * width * 3; i += 3) {
     
        result[i/3] =  ((a *img[i] + b * img[i+1] + c * img[i+2])+(float)0.5);
    }
}
int main() {
    // Open the input file in binary mode
    FILE* inputFile = fopen("mandjdid.ppm", "rb");
    if (!inputFile) {
        perror("Error opening input file");
        return 1;
    }

    // Read and validate the PPM header
    char magic[3];
    fscanf(inputFile, "%2s", magic);
    if (magic[0] != 'P' || magic[1] != '6') {
        fprintf(stderr, "Invalid PPM format\n");
        fclose(inputFile);
        return 1;
    }

    int width, height, maxColor;
    fscanf(inputFile, "%d %d %d", &width, &height, &maxColor);
    fgetc(inputFile); // Consume newline character

    // Allocate memory for pixel data
    uint8_t* pixels = (uint8_t*)malloc(width * height * 3 * sizeof(uint8_t));
    if (!pixels) {
        perror("Memory allocation error");
        fclose(inputFile);
        return 1;
    }

    // Read pixel data
    fread(pixels, sizeof(uint8_t), width * height * 3, inputFile);

    // Allocate memory for the grayscale results
    uint8_t* resultSIMD = (uint8_t*)_mm_malloc(width * height * sizeof(uint8_t), 16);
    uint8_t* resultNonSIMD = (uint8_t*)malloc(width * height * sizeof(uint8_t));
    uint8_t* resultPrecomputed = (uint8_t*)malloc(width * height * sizeof(uint8_t));
    if (!resultSIMD || !resultNonSIMD || !resultPrecomputed) {
        perror("Memory allocation error");
        free(pixels);
        fclose(inputFile);
        if (resultSIMD) _mm_free(resultSIMD);
        if (resultNonSIMD) free(resultNonSIMD);
        if (resultPrecomputed) free(resultPrecomputed);
        return 1;
    }

    clock_t start, end;
    double cpu_time_used;

    // Run SIMD-based grayscale conversion and time it
    start = clock();
    grayscale(pixels, width, height, 0.299, 0.587, 0.114, resultSIMD);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("SIMD grayscale conversion took %f seconds\n", cpu_time_used);

    // Run non-SIMD-based grayscale conversion and time it
    start = clock();
    grayscale1(pixels, width, height, 0.299, 0.587, 0.114, resultNonSIMD);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Non-SIMD grayscale conversion took %f seconds\n", cpu_time_used);

    // Run precomputed-table grayscale conversion and time it
    start = clock();
    grayscalePrecomputed(pixels, width, height, 0.299, 0.587, 0.114, resultPrecomputed);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Precomputed table grayscale conversion took %f seconds\n", cpu_time_used);

    // Compare SIMD and Non-SIMD results for discrepancies
    int diffCount = 0;
    for (int i = 0; i < width * height; ++i) {
        if (resultSIMD[i] != resultNonSIMD[i]) {
       //     printf("%d ",resultNonSIMD[i]-resultSIMD[i]);
            diffCount++;
        }
    }
    printf("Number of differing pixels between SIMD and Non-SIMD: %d\n", diffCount);
    float po=(diffCount*100)/(width*height);
    printf("percent of inaccuracy is %f \n",po);

    //test the difference between prelookup and naiv 
    // Compare Non-SIMD and Precomputed-Table results for discrepancies
int diffCountNonSIMDPrecomputed = 0;
for (int i = 0; i < width * height; ++i) {
    if(resultNonSIMD[i]!=resultPrecomputed[i]) diffCountNonSIMDPrecomputed++;
    int diff = (int)resultNonSIMD[i] - (int)resultPrecomputed[i];
    if (diff > 1 || diff < -1) {
        printf("Pixel %d: Non-SIMD = %d, Precomputed = %d, Diff = %d\n", i, resultNonSIMD[i], resultPrecomputed[i], diff);
        
    }
}
printf("Number of differing pixels with more than 1 difference between Non-SIMD and Precomputed: %d\n", diffCountNonSIMDPrecomputed);
float percentDiffNonSIMDPrecomputed = (diffCountNonSIMDPrecomputed * 100.0f) / (width * height);
printf("Percentage of significant differences: %f%%\n", percentDiffNonSIMDPrecomputed);
    // Write SIMD result to output file
    FILE* outputFile = fopen("ergebnis.pgm", "wb");
    if (!outputFile) {
        perror("Error opening output file for writing");
        free(pixels);
        _mm_free(resultSIMD);
        free(resultNonSIMD);
        free(resultPrecomputed);
        fclose(inputFile);
        return 1;
    }
    fprintf(outputFile, "P5\n%d %d\n255\n", width, height);
    fwrite(resultSIMD, sizeof(uint8_t), width * height, outputFile);

    // Close files and free memory
    fclose(inputFile);
    fclose(outputFile);
    free(pixels);
    _mm_free(resultSIMD);
    free(resultNonSIMD);
    free(resultPrecomputed);

    return 0;
}