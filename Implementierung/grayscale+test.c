#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

void grayscale(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* result) {
    for (size_t i = 0; i < height * width * 3; i += 3) {
        uint8_t R = img[i];
        uint8_t G = img[i + 1];
        uint8_t B = img[i + 2];
        float D = (R * a + G * b + B * c) / (a + b + c);
        result[i] = (uint8_t)(round(D));
        result[i + 1] = (uint8_t)(round(D));
        result[i + 2] = (uint8_t)(round(D));
    }
}

int main() {
    FILE* inputFile = fopen("mandrill.ppm", "rb"); // Open the file in binary mode

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

    // Allocate memory for pixel data
    uint8_t* pixels = (uint8_t*)malloc(width * height * 3 * sizeof(uint8_t));

    if (!pixels) {
        perror("Memory allocation error");
        fclose(inputFile);
        return 1;
    }

    // Consume newline character following the maxColor value
    fgetc(inputFile);

    // Read pixel data
    fread(pixels, sizeof(uint8_t), width * height * 3, inputFile);

    // Example: Print the color of the pixel at position (0, 0)
    printf("Color at (0, 0): R=%d, G=%d, B=%d\n", pixels[0], pixels[1], pixels[2]);

    // Allocate memory for the result
    uint8_t* result = (uint8_t*)malloc(width * height * 3 * sizeof(uint8_t));

    if (!result) {
        perror("Memory allocation error");
        free(pixels);
        fclose(inputFile);
        return 1;
    }

    // Call the grayscale function
    grayscale(pixels, width, height, 0.299, 0.587, 0.114, result);

    // Open the output file
    FILE* outputFile = fopen("ergebnis.ppm", "wb"); // Open the file in binary write mode

    if (!outputFile) {
        perror("Error opening output file for writing");
        free(pixels);
        free(result);
        fclose(inputFile);
        return 1;
    }

    // Write the PPM header
    fprintf(outputFile, "P6\n%d %d\n255\n", width, height);

    // Write pixel data
    fwrite(result, sizeof(uint8_t), width * height * 3, outputFile);

    // Close the files and free memory
    fclose(inputFile);
    fclose(outputFile);
    free(pixels);
    free(result);

    return 0;
}
