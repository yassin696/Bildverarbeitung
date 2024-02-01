#define _POSIX_C_SOURCE 199309L
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

void grayscale_naive(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    //initialise the R G and B values
    float* result = (float*) tmp;
    uint8_t R;
    uint8_t G;
    uint8_t B;
    for (size_t i = 0; i < height * width * 3; i+=3) {
        //loading the colours' variables 
        R= img[i];
        G= img[i+1];
        B= img[i+2];
        //calculate the grayscale D value round it and store it to its corresponding place
        result[i/3] = round(a * R + b * G + c * B);
    }
}

// lookup grayscale
// Function to generate a precomputed table for a given coefficient
void gentable(float coeff, float* table) {
    for (size_t i = 0; i < 256; ++i) {
        table[i] = i * coeff;
    }
}

void grayscale_lookup(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    // Generate tables for coefficients
    float tableA[256];
    float tableB[256];
    float tableC[256];
    gentable(a, tableA);
    gentable(b, tableB);
    gentable(c, tableC);
    float* result = (float*)(tmp);
    //initialise the colour values which are also indexes to the prelookup tables
     uint8_t R;
     uint8_t G;
     uint8_t B;
    for (size_t i = 0; i < height * width*3; i+=3) {
        //load R G and B 
         R = img[i];
         G = img[i+ 1];
         B = img[i+ 2];
        // Lookup precomputed values from tables
        float D = round(tableA[R] + tableB[G] + tableC[B]);       
        result[i/3] = D;
    }
}

// simd grayscale
void grayscale_simd(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    float* result = (float*)(tmp);
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
        __m128 D = _mm_add_ps(_mm_add_ps(Ra, Gb), Bc);
        __m128 rounding =_mm_set1_ps(0.5f);
        __m128 D1 = _mm_add_ps(D, rounding);
        // Round the floating-point values to the nearest integer
        __m128 D_rounded = _mm_floor_ps(D1);
        _mm_storeu_ps(result + i/3, D_rounded);
     }

    // Handle any remaining pixels
    for (; i < height * width * 3; i += 3) {
     
        result[i/3] =round(  a *img[i] + b * img[i+1] + c * img[i+2]);
    }
}


// naive interpolation calculation 
void interpolation_calculation_naive(size_t width, size_t height, size_t factor, uint8_t* tmp, uint8_t* result) {
    float* inputArray = (float*) tmp;
    // Calculate scaled width and height
    size_t scaledW = width * factor;  
    size_t scaledH = height * factor; 

    // Inverse factor for interpolation
    float invFactor = 1.0f / factor;

    // Initialize variables outside the loops
    float gy, oneMinusGy, gx, oneMinusGx;
    size_t yIdx, yNextIdx, xIdx, xNextIdx;
    float Q00, Qs0, Q0s, Qss;

    // Iterate over each pixel in the scaled image
    for (size_t yScaled = 0; yScaled < scaledH; yScaled++) {
        //precalculate  intermediate values
        gy = (float)(yScaled % factor) * invFactor;
        oneMinusGy = 1.0f - gy;
        //calculate indexes for the given 
        yIdx = yScaled / factor;
        yNextIdx = (yIdx + 1) % height;

        for (size_t xScaled = 0; xScaled < scaledW; xScaled++) {
            gx = (float)(xScaled % factor) * invFactor;
            oneMinusGx = 1.0f - gx;
            xIdx = xScaled / factor;
            xNextIdx = (xIdx + 1) % width;

            // Retrieve pixel values from the four surrounding corners
            Q00 = inputArray[yIdx * width + xIdx];
            Qs0 = inputArray[yIdx * width + xNextIdx];
            Q0s = inputArray[yNextIdx * width + xIdx];
            Qss = inputArray[yNextIdx * width + xNextIdx];

            // Interpolate and store the result
            result[yScaled * scaledW + xScaled] = (uint8_t)((oneMinusGy * (oneMinusGx * Q00 + gx * Qs0) + gy * (oneMinusGx * Q0s + gx * Qss)) + 0.5f);
        }
    }
}

// algorithmically optimized interpolation calculation
void interpolation_calculation_algorithm_optimized(size_t width, size_t height, size_t factor, uint8_t* tmp, uint8_t* result) {
    float* inputArray = (float*) tmp;
    float Q00, Qs0, Q0s, Qss; // Variables to store pixel values for interpolation
    float val1, val2; // Intermediate values for bilinear interpolation

    // Iterate over height
    for (size_t h = 0; h < height; h++) {
        // Iterate over width
        for (size_t b = 0; b < width; b++) {
            // Get the pixel values for the current and adjacent pixels
            Q00 = inputArray[b + h * width]; // Current pixel 
            Qs0 = inputArray[(b + 1) % width + h * width]; // Right neighbor
            Q0s = inputArray[b + ((h + 1) % height) * width]; // Bottom neighbor
            Qss = inputArray[(b + 1) % width + ((h + 1) % height) * width]; // Diagonal neighbor

            // Perform bilinear interpolation for each subpixel within the current sector
            for (size_t f1 = 0; f1 < factor; f1++) {
                // Calculate intermediate values
                val1 = ((float)(factor - f1) / factor) * Q00 + ((float)f1 / factor) * Qs0;
                val2 = ((float)(factor - f1) / factor) * Q0s + ((float)f1 / factor) * Qss;

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

// simd interpolation calculation
void processEdge( float* inputArray, uint8_t* result, size_t width, size_t height, size_t factor, size_t startWidth, size_t startHeight) {
    float Q00, Qs0, Q0s, Qss; // Variables to store pixel values for interpolation
    float val1, val2; // Intermediate values for bilinear interpolation
    // Iterate over height
    for (size_t h = startHeight; h < height; h++) {
        // Iterate over width
        for (size_t b = startWidth; b < width; b++) {
            // Get the pixel values for the current and adjacent pixels
            Q00 = inputArray[b + h * width]; // Current pixel
            Qs0 = inputArray[(b + 1) % width + h * width]; // Right neighbor
            Q0s = inputArray[b + ((h + 1) % height) * width]; // Bottom neighbor
            Qss = inputArray[(b + 1) % width + ((h + 1) % height) * width]; // Diagonal neighbor

            // Perform bilinear interpolation for each subpixel within the current sector
            for (size_t f1 = 0; f1 < factor; f1++) {
                // Calculate intermediate values
                val1 = ((float)(factor - f1) / factor) * Q00 + ((float)f1 / factor) * Qs0;
                val2 = ((float)(factor - f1) / factor) * Q0s + ((float)f1 / factor) * Qss;

                // Calculate the final interpolated values and assign them to the result
                for (size_t f = 0; f < factor; f++) {
                    result[(h * factor + f) * width * factor + (b * factor + f1)] =
                        (((float)(factor - f) / factor) * val1 +
                        ((float)f / factor) * val2) + 0.5f; // Adding 0.5f for rounding
                }
            }
        }
    }
}


void interpolation_calculation_simd( size_t width, size_t height, size_t factor, uint8_t* tmp, uint8_t* result) {
    float* inputArray = (float*) tmp;
    size_t scaledWidth = width * factor;//calculate new width
    float invFactor = 1.0f / factor;  // Inverse of scaling factor for coefficient calculation

    // Allocate memory for precomputed interpolation coefficients a b c d
    float* coeffTable = (float*)malloc(sizeof(float) * factor * factor * 4); 

    float gx,gy;//initialise intermediate values
    size_t index,x,y,sectorBaseX,coeffIndex,resultIndex;//initlialise indexes
    __m128 Q00,Qs0,Q0s,Qss,a,b,c,d,interpolated,coeff;//initialise intermediate values and final interpolated value
    for (size_t i = 0; i < factor; i++) {
        for (size_t j = 0; j < factor; j++) {
             gx = j * invFactor;
             gy = i * invFactor;
             index = (i * factor + j) * 4;
            // Precompute bilinear interpolation coefficients for each sub-pixel position
            coeffTable[index] = (1 - gx) * (1 - gy);     // Coefficient a
            coeffTable[index + 1] = gx * (1 - gy);       // Coefficient b
            coeffTable[index + 2] = (1 - gx) * gy;       // Coefficient c
            coeffTable[index + 3] = gx * gy;             // Coefficient d
        }
    }

    //  SIMD processing loop over the image
    for (y = 0; y < height - 1; y++) {
        for ( x = 0; x < (width - 1)-((width-1)%4); x += 4) { 
            // Load pixel values for Q00, Q0s, Qs0, Qss using SIMD
             Q00 = _mm_loadu_ps(&inputArray[y * width + x]);
             Q0s = _mm_loadu_ps(&inputArray[(y + 1) * width + x]);
             Qs0 = _mm_loadu_ps(&inputArray[y * width + x + 1]);
             Qss = _mm_loadu_ps(&inputArray[(y + 1) * width + x + 1]);

            // Iterate over each sub-pixel position within the pixel
            for (size_t i = 0; i < factor; i++) {
                for (size_t j = 0; j < factor; j++) {
                     coeffIndex = (i * factor + j) * 4;
                     coeff = _mm_loadu_ps(&coeffTable[coeffIndex]);

                    // Load the precomputed coefficients into SIMD 128bit registers
                     a = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(0, 0, 0, 0));
                     b = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(1, 1, 1, 1));
                     c = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(2, 2, 2, 2));
                     d = _mm_shuffle_ps(coeff, coeff, _MM_SHUFFLE(3, 3, 3, 3));

                    // Perform bilinear interpolation using SIMD
                     interpolated = _mm_add_ps(_mm_add_ps(_mm_mul_ps(a, Q00), _mm_mul_ps(b, Qs0)), _mm_add_ps(_mm_mul_ps(c, Q0s), _mm_mul_ps(d, Qss)));

                    // Store the interpolated results into the result buffer
                    for (size_t k = 0; k < 4; k++) {
                         sectorBaseX = (x + k) * factor;
                         resultIndex = (y * factor + i) * scaledWidth + sectorBaseX + j;
                        result[resultIndex] = (interpolated[k] + 0.5f); // Convert to uint8_t with rounding
                    }
                }
            }
        }
    }

    // Process the edge cases by calling 
    processEdge(inputArray, result, width, height, factor,(width - 1)-((width-1)%4), 0); // Rightmost column
    processEdge(inputArray, result, width, height, factor, 0, height-1); // Bottom row

    free(coeffTable); // Free the allocated memory for coefficients
}


// convert float to uint8_t (by grayscale)
void scale1_handling(const float* tmp, uint8_t* result, size_t width, size_t height) {
    for (size_t i = 0; i < width * height; i++) {
        result[i] = (uint8_t) tmp[i]; // Cast each float value to uint8_t
    }
}

// naive interpolate
void interpolate_naive(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    grayscale_naive(img, width, height, a, b, c, tmp);
    if (scale_factor == 1) {
        scale1_handling((float*) tmp, tmp, width, height);
        return;
    }
    interpolation_calculation_naive(width, height, scale_factor, tmp, result);
}

// algorithmically optimized interpolate
void interpolate_algorithm_optimized(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    grayscale_lookup(img, width, height, a, b, c, tmp);
    if (scale_factor == 1) {
        scale1_handling((float*) tmp, tmp, width, height);
        return;
    }
    interpolation_calculation_algorithm_optimized(width, height, scale_factor, tmp, result);
}

// simd interpolate
void interpolate_simd(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
  //  grayscale_simd(img, width, height, a, b, c, tmp);
    grayscale_simd(img, width, height, a, b, c, tmp);
    if (scale_factor == 1) {
        scale1_handling((float*) tmp, tmp, width, height);
        return;
    }
    interpolation_calculation_simd(width, height, scale_factor, tmp, result);
}

// interpolate (standardversion) (the best interpolation is interpolate_simd)
void interpolate(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    interpolate_simd(img, width, height, a, b, c, scale_factor, tmp, result);
}

// framework
// checking coefficients and precalculate them (a=a/(a+b+c), b=b/(a+b+c), c=c/(a+b+c))
void check_coefficients(float* a, float* b, float* c) {
    float sum = *a + *b + *c;
    if (sum <= 0) {
        *a = 0.299;
        *b = 0.587;
        *c = 0.114;
        sum = 1.0;
        printf("The sum of a, b and c is smaller or equal 0. Therefore, the default value will be used.\n");
    }
    if (sum != 1) {
        *a = *a / sum;
        *b = *b / sum;
        *c = *c / sum;
    }
}

// checking scaling
void check_scaling(int* scaling) {
    if (*scaling < 1) {
        printf("The scale factor is smaller than 1, so the default value 2 will be used.\n");
        *scaling = 2;
    }
}

// checking output filename
char* check_output(char* outputFileName) {
    // Characters like \ / : * ? " < > | . are not allowed in filenames
    // Source: https://verwaltung.uni-koeln.de/stabsstelle01/content/benutzerberatung/it_faq/windows/faqitems163122/index_ger.html
    int length = strlen(outputFileName);
    char* filename = outputFileName;
    if (length == 0 || length > 255) {
        // filename too long or does not exists
        printf("Filename is too long or does not exists. Therefore the default value 'output.pgm' will be used.\n");
        outputFileName = "output.pgm";
        return outputFileName;
    }
    if (filename[0] == '.' || filename[length-1] == '.') {
        // Filename is not allowed to start or end with '.'
        printf("Filename are not allowed to start or end with '.'. Therefore the default value 'output.pgm' will be used.\n");
        outputFileName = "output.pgm";
        return outputFileName;
    }
    int whitespace = 1;
    for (int i = 0; i < length; i++) {
        if (filename[i] != ' ') {
            whitespace = 0;
        }
        if (filename[i] == '/' || filename[i] == '\\' || filename[i] == ':' || filename[i] == '*' || filename[i] == '?' || filename[i] == '"' || filename[i] == '<' || filename[i] == '>' || filename[i] == '|') {
            // Filename contains not allwoed characters
            printf("Filename contains not allwoed characters. Therefore the default value 'output.pgm' will be used.\n");
            outputFileName = "output.pgm";
            return outputFileName;
        }
    }
    if (whitespace == 1) {
        // Filename only contains whitespaces
        printf("Filename only contains whitespaces. Therefore the default value 'output.pgm' will be used.\n");
        outputFileName = "output.pgm";
        return outputFileName;
    }
    if (strcmp(strrchr(outputFileName, '\0') - 4, ".pgm") != 0) {
        // Filename should end with '.pgm' 
        printf("Filename should end with '.pgm'. Therefore '.pgm' will be concatenated on the given filename.\n");
        filename[length] = '.';
        filename[length+1] = 'p';
        filename[length+2] = 'g';
        filename[length+3] = 'm';
        filename[length+4] = '\0';
        outputFileName = filename;
        return outputFileName;
        
    }
    return outputFileName;
}

// checking user inputs
void check_user_input(float* a, float* b, float* c, int* scaling, char* outputFileName) {
    // Checking and calculating of coefficients a, b, c (a=a/(a+b+c), b=b/(a+b+c), c=c/(a+b+c))
    check_coefficients(a, b, c);
    //printf("The coefficients a, b, c: %f %f %f\n", *a, *b, *c); // test coefficients

    // Checking scale_factor
    check_scaling(scaling);
    //printf("Scaling factor: %i\n", *scaling); // test scale_factor 
        
    // Checking output filename
    outputFileName = check_output(outputFileName);
    //printf("Output Filename: %s\n", outputFileName); // test outputfilename
}

// read ppm header - return 1 if the header is false
void skip_comment(FILE *inputFile) {
    // skip possible comments between the values that need to be read for the header
    int ch;
    int itemsRead = fscanf(inputFile, " "); // skip withespace
    if(itemsRead == EOF){
        errno = EINVAL;
        perror("Error by skipping whitespace");
        fclose(inputFile);
        return ;
    }
    while ((ch = fgetc(inputFile)) == '#') { // comment begins
      while ((ch = fgetc(inputFile)) != '\n'); // skip comment 
      itemsRead = fscanf(inputFile, " "); // skip withespace
      if(itemsRead != 1){
        errno = EINVAL;
        perror("Error by skipping whitespace");
        fclose(inputFile);
        return ;
    }
    }
    ungetc(ch, inputFile); // go last readed character back 
}

int read_ppm_header(FILE* inputFile, int* width, int* height, int* maxColorValue) {
    int char1 = fgetc(inputFile);
    int char2 = fgetc(inputFile);
    if (char1 != 'P' || char2 != '6') {
        // Error because wrong picture format 
        errno = EINVAL;
        perror("Invalid Pictureformat. A PPM P6 Picture is expected. For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    skip_comment(inputFile);
    int itemsRead ;
    itemsRead = fscanf(inputFile, "%d", width);
    if(itemsRead != 1){
        errno = EINVAL;
        perror("Invalid width Format. For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    if (*(width) <= 0) {
        // Error because wrong width
        errno = EINVAL;
        perror("Invalue value for width. Width must be greater than 0. For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    skip_comment(inputFile);
    itemsRead = fscanf(inputFile, "%d", height);
    if(itemsRead!=1){
        errno = EINVAL;
        perror("Invalid height Format. For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    if (*(height) <= 0) {
        // Error because wrong height
        errno = EINVAL;
        perror("Invalid value for height. Height must be greater than 0. For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    skip_comment(inputFile);
    itemsRead = fscanf(inputFile, "%d", maxColorValue); 
    if(itemsRead !=1 ){
        errno = EINVAL;
        perror("Invalid maxColorValue Format. For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    if (*(maxColorValue) > 255 || *(maxColorValue) < 128) {
        // Error because wrong maxColorValue
        errno = EINVAL;
        perror("Invalue value for max color. A value between 0 and 255 is expected (the picture should be 24bpp). For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    //printf("Width, height, maxColorValue: %i, %i, %i\n", *(width), *(height), *(maxColorValue));
    return 0;
}

// read ppm
uint8_t* read_ppm(char* inputFileName, int* width, int* height, int* imageSize) {
    // check inputfile 
    if (strcmp(strrchr(inputFileName, '\0') - 4, ".ppm") != 0) {
        // The inputfile does not end with ".ppm" -> false file format
        errno = EIO;
        perror("Error: Wrong file format. A PPM-Picutre is expected. For more information please use the option -h | --help");
        return NULL;
    }
    printf("Inputfile: %s\n", inputFileName);
       
    // read the ppm-header
    FILE* inputFile = fopen(inputFileName, "rb");
    if (inputFile == NULL) {
        // Error while trying open the giving picture
        perror("The giving picutre can not be opened. For more information please use the option -h | --help");
        return NULL;
    }
    int maxColorValue;
    if (read_ppm_header(inputFile, width, height, &maxColorValue) == 1) {
        return NULL;
    }
    //printf("Width, height, maxColorValue: %i, %i, %i\n", *(width), *height, *maxColorValue);

    // read the picture
    *(imageSize) = *(width) * *(height);
    uint8_t* pixels = (uint8_t*)malloc(*(imageSize) * 3 * sizeof(uint8_t));
    if (pixels == NULL) {
        // Error by allocating memory
        perror("Error by allocating memory\n");
        fclose(inputFile);
        return NULL;
    }
    fgetc(inputFile); // new line character after end of the header
    int itemsRead = fread(pixels, sizeof(uint8_t), *(imageSize)*3, inputFile);
    if (itemsRead < *(imageSize)*3) {
        if (feof(inputFile)) {
            // End of file was reached before reading the expected amount of data.
            perror("Error reading input file");
            return NULL;
        } else if (ferror(inputFile)) {
            // An error occurred during reading.
            perror("Error reading input file");
            return NULL;
        }
   }
    fclose(inputFile);
    //uint8_t* check_pixels = pixels;
    //for (int i = 0; i < *(imageSize)*3; i++) {
    //    if (check_pixels[i] > maxColorValue) {
    //        // Error because a pixels is greater than the maxColorValue
    //        errno = EINVAL;
    //        perror("Invalid value for pixels. RGB-Value of a pixel must be smaller than the max color. For more information please use the option -h | --help");
    //        return NULL;
    //    }
    //}
    return pixels;
}

// write ppm
int write_ppm(char* outputFileName, int width, int height, size_t scaling, uint8_t* result) {
    FILE* outputFile = fopen(outputFileName, "wb");
    if (outputFile == NULL) {
        // Error opening the outputfile 
        perror("Error by opening the outputfile. For more information please use the option -h | --help.");
        fclose(outputFile);
        return 1;
    }
    fprintf(outputFile, "P5\n");
    fprintf(outputFile, "%ld %ld\n", (width*scaling), (height*scaling));
    fprintf(outputFile, "255\n");
    fwrite(result, sizeof(uint8_t), (width*scaling)*(height*scaling)*sizeof(uint8_t), outputFile);
    fclose(outputFile);
    return 0;
}

// Description of all options of the program and usage examples are displayed
void print_usage() {
    printf("Usage: interpolate [options] [target] ...\n");
    printf("       default usage: ./interpolate input.pgm\n");
    printf("       using specified implementation: ./interpolate -V 1 input.pgm\n");
    printf("       measuaring runtime: ./interpolate -B 10 input.pgm\n");
    printf("       specify outputfile: ./interpolate -o result.pgm input.pgm\n");
    printf("       using costumized coefficients: ./interpolate --coeffs 0.2,0.5,0.3 input.pgm\n");
    printf("       set a scale factor: ./interpolate -f 3 input.pgm\n");
    printf("       display all options with examples: ./interpolate --help\n");

    printf("Options: \n");
    printf("    -V<Number> :\n");
    printf("                The implementation, that will be used.\n");
    printf("                -V 0 will use the standard implementation (simd implementation), -V 1 the naive implementation, \n");
    printf("                -V 2 the algorithmically optimized implementation and -V 3 the SIMD implementation. \n");
    printf("                If this output is not called, the standard implementation will also be called.\n");
    printf("    -B<Number> :\n");
    printf("                If this option is called, the duration for the specified implementation will be measured. \n")
    printf("                This optional argument indicates the number of the repetitions of the function calling.\n");
    printf("    -o<Filename> :\n");
    printf("                  Output filename If this option is not called, the default value 'output.pgm' will be used.\n");
    printf("    --coeffs<FP Number>,<FP Number>,<FP Number> :\n");
    printf("                                                 The coefficients a, b and c for the grayscale. \n");
    printf("                                                 If this option is not called, the default value 0.299, 0.587 und 0.114 will be used.\n");
    printf("    -f<Number> :\n");
    printf("                Scaling factor. If this option is not called, the default value 2 will be used.\n");
    printf("    -h|--help :\n");
    printf("                A description of all options of the program with usage examples will be displayed.\n");
}

// main - framework programm
int main(int argc, char **argv){
    // declaration of variables
    int option;
    int implementation = 0; // if implementation is 0, then the standard implementation will be used
    int benchmark = 0; // 1, if the programm should measure the duration
    int repetitions = 1; // number of repetitions
    char* inputFileName = NULL; // input file
    char* outputFileName = "output.pgm"; // standard output file, if there is no valide user input
    float a = 0.299; // coefficients a, b, c, if there are no valide user inputs
    float b = 0.587;
    float c = 0.114;
    int scaling = 2;

    // parse command line arguments
    static struct option long_options[] = {
            {"coeffs", required_argument, 0, 'c'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
    };

    while ((option = getopt_long(argc, argv, "V:B:o:f:h", long_options, NULL)) != -1) {
        switch (option) {
            case 0:
                sscanf(optarg, "%f,%f,%f", &a, &b, &c);
                break;
            case 'V':
                implementation = atoi(optarg);
                break;
            case 'B':
                benchmark = 1;
                repetitions = atoi(optarg);
                break;
            case 'o':
                outputFileName = optarg;
                break;
            case 'c':
                sscanf(optarg, "%f,%f,%f", &a, &b, &c);
                break;
            case 'f':
                scaling = atoi(optarg);
                break;
            case 'h':
                // Description of all options of the program and usage examples are displayed (the programm will end after that)
                print_usage();
                exit(0);
            default:
                // Error because of missing entries
                errno = EIO;
                perror("All necessary parameters are missing or wrong use of parameters. For more information please use the option -h | --help.\n");
                exit(1);
        }
    }

    if (optind < argc) {
        inputFileName = argv[optind];
    } else {
        // Error because the inputfile is missing
        errno = EIO;
        perror("Error: The inputfile is missing. For more information please use the option -h | --help");
        return 1;
    }

    // program with checking user inputs and calculation
    // prepare for the grayscale and interpolation

    // check user inputs
    check_user_input(&a, &b, &c, &scaling, outputFileName);
    //printf("%s und %s", outputFileName, inputFileName);

    // read ppm
    int width; int height; int imageSize;
    const uint8_t* pixels = read_ppm(inputFileName, &width, &height, &imageSize);
    if (pixels == NULL) {
        return 1; // error by reading ppm
    }
        
    // allocating memories 
    uint8_t* temp = (uint8_t*)malloc((imageSize * sizeof(float)) + (imageSize * scaling * scaling * sizeof(uint8_t)));
    if (temp == NULL) {
        // Error by allocating memory
        perror("Error by allocating memory");
        return 1;
    }
    uint8_t* result = temp + (imageSize * sizeof(float));
    if (result == NULL) {
        // Error by allocating memory
        perror("Error by allocating memory");
        return 1;
    }
    // result is a part of the allocated memory for temp. The interpolated picture will be saved by result, the grayscale by temp
    
    // Calculation with different implementations
    switch (implementation) {
        case 1:
            // naive implementation
            printf("Naive implementation with coefficients: %f, %f, %f, scaling factor: %i, output filename: %s and inputfile: %s\n", a, b, c, scaling, outputFileName, inputFileName);
            if (benchmark == 1) {
                if (repetitions < 1) { repetitions = 1; }  // repetition can not be smaller than 1
                // measuaring duration
                // using code of gra videos
                struct timespec start;
                clock_gettime(CLOCK_MONOTONIC, &start);
                for (int i = 1; i < repetitions; i++) {
                    // call the functions 
                    interpolate_naive(pixels, width, height, a, b, c, scaling, temp, result);
                }
                struct timespec end;
                clock_gettime(CLOCK_MONOTONIC, &end);
                double time = end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
                double averageTime = time / repetitions;
                printf("Bechmark with %i iterations: %fs. The average time per iteration is: %fs.\n", repetitions, time, averageTime);
                break;
            }
            
            // grayscale and interpolation  
            interpolate_naive(pixels, width, height, a, b, c, scaling, temp, result);
            break;
        case 2:
            // algorithmically optimized
            printf("Algorithmically optimized implementation with coefficients: %f, %f, %f, scaling factor: %i, output filename: %s and inputfile: %s\n", a, b, c, scaling, outputFileName, inputFileName);
            if (benchmark == 1) {
                if (repetitions < 1) { repetitions = 1; }  // repetition can not be smaller than 1
                // measuaring duration
                // using code of gra videos
                struct timespec start;
                clock_gettime(CLOCK_MONOTONIC, &start);
                for (int i = 1; i < repetitions; i++) {
                    // call the functions 
                    interpolate_algorithm_optimized(pixels, width, height, a, b, c, scaling, temp, result);
                }
                struct timespec end;
                clock_gettime(CLOCK_MONOTONIC, &end);
                double time = end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
                double averageTime = time / repetitions;
                printf("Bechmark with %i iterations: %fs. The average time per iteration is: %fs.\n", repetitions, time, averageTime);
                break;
            }
            
            // grayscale and interpolation  
            interpolate_algorithm_optimized(pixels, width, height, a, b, c, scaling, temp, result);
            break;
        case 3:
            // simd implementation
            printf("SIMD implementation with coefficients: %f, %f, %f, scaling factor: %i, output filename: %s and inputfile: %s\n", a, b, c, scaling, outputFileName, inputFileName);
            if (benchmark == 1) {
                if (repetitions < 1) { repetitions = 1; }  // repetition can not be smaller than 1
                // measuaring duration
                // using code of gra videos
                struct timespec start;
                clock_gettime(CLOCK_MONOTONIC, &start);
                for (int i = 1; i < repetitions; i++) {
                    // call the functions 
                    interpolate_simd(pixels, width, height, a, b, c, scaling, temp, result);
                }
                struct timespec end;
                clock_gettime(CLOCK_MONOTONIC, &end);
                double time = end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
                double averageTime = time / repetitions;
                printf("Bechmark with %i iterations: %fs. The average time per iteration is: %fs.\n", repetitions, time, averageTime);
                break;
            }
            
            // grayscale and interpolation  
            interpolate_simd(pixels, width, height, a, b, c, scaling, temp, result);
            break;
        default:
            // standard implementation
            printf("Standard implementation with coefficients: %f, %f, %f, scaling factor: %i, output filename: %s and inputfile: %s\n", a, b, c, scaling, outputFileName, inputFileName);
            if (benchmark == 1) {
                if (repetitions < 1) { repetitions = 1; }  // repetition can not be smaller than 1
                // measuaring duration
                // using code of gra videos
                struct timespec start;
                clock_gettime(CLOCK_MONOTONIC, &start);
                for (int i = 1; i < repetitions; i++) {
                    // call the functions 
                    interpolate(pixels, width, height, a, b, c, scaling, temp, result);
                }
                struct timespec end;
                clock_gettime(CLOCK_MONOTONIC, &end);
                double time = end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
                double averageTime = time / repetitions;
                printf("Bechmark with %i iterations: %fs. The average time per iteration is: %fs.\n", repetitions, time, averageTime);
                break;
            }
            
            // grayscale and interpolation  
            interpolate(pixels, width, height, a, b, c, scaling, temp, result);
            break;
    }


    // Saving the result picture 
    if (scaling == 1) {
        if (write_ppm(outputFileName, width, height, scaling, temp) == 1) {
            free(temp);
            return 1;
        }
    } else {
        if (write_ppm(outputFileName, width, height, scaling, result) == 1) { 
            free(temp);
            return 1;
        }
    }
    free(temp);
    return 0;
}   
