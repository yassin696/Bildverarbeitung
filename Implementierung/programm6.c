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

// naive grayscale-solution
void grayscale_naiv(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    size_t i = 0;
    int j = 0;
    while (i < height * width * 3) {
        uint8_t R = img[i];
        uint8_t G = img[i + 1];
        uint8_t B = img[i + 2];
        uint8_t D = round(R * a + G * b + B * c);
        tmp[j] = D;
        i += 3;
        j += 1;
    }
}

// grayscale optimized with look-up tables
void gentable(float coeff, uint8_t* array) {
    for (int i = 0; i < 256; ++i) {
        array[i] = (uint8_t)(i * coeff);
    }
}

void grayscale_look_up(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    uint8_t* tableA = malloc(256 * sizeof(uint8_t));
    uint8_t* tableB = malloc(256 * sizeof(uint8_t));
    uint8_t* tableC = malloc(256 * sizeof(uint8_t));
    gentable(a, tableA);
    gentable(b, tableB);
    gentable(c, tableC);

    size_t size = height * width;
    for (size_t i = 0; i < size; ++i) {
        uint8_t R = img[i];
        uint8_t G = img[size + i];
        uint8_t B = img[size + size + i];

        // Lookup precomputed values from tables
        uint8_t weightedSum = tableA[R] + tableB[G] + tableC[B];
        tmp[i] = weightedSum;
    }

    free(tableA);
    free(tableB);
    free(tableC);
}


// grayscale optimiyed with SIMD


// naive bilinear interpolation
void bilinear_interpolate_naive(uint8_t* tmp, size_t width, size_t height, size_t scale_factor, uint8_t* result) {
    // Einzelne Sektoren werden bearbeitet
    for(size_t sektorh = 0; sektorh < height;sektorh++) {
        for (size_t sektorb = 0; sektorb < width; sektorb++) {
            //Sektoren a,b,c,d Werte, also (0,0), (0,s), (s,0) und (s,s) werden für alle berechnet
            float a = tmp[sektorh * width + sektorb];
            float b = tmp[sektorh * width + (sektorb + 1) % width];
            float c = tmp[((sektorh + 1) % height) * width + sektorb];
            float d = tmp[((sektorh + 1) % height) * width + (sektorb + 1) % width];

            // Für jeden einzelnen Wert in den Quadraten berechne nach Formel den Wert
            for(size_t y=0;y<scale_factor;y++) {
                for (size_t x = 0; x < scale_factor; x++) {

                    // berechne position im finalen AllozSpeicher
                    int pos = (x + scale_factor * sektorb) + ((sektorh * scale_factor + y) * width * scale_factor);

                    // da wert x=0, y=0 immer den a-Wert ergibt kann dieser auch direkt eingetragen werden
                    // die berechnung füllt ja nur das quadrat für eins kleiner von s aus, dass keine Werte doppelt berechnet werden
                    if(x == 0 && y == 0) {
                        result[pos] = (uint8_t) a;
                        continue;
                    }

                    // berechne Wert
                    float polwert = (a * (scale_factor-y) * (scale_factor-x) ) + ( c * y * ( scale_factor-x) ) + ( b * (scale_factor-y) * x ) + ( d * y * x );
                    // multipliziere mit (1 / s*s)
                    polwert = polwert / (scale_factor * scale_factor);
                    //printf("[%i;%i] %i - %i\n",x,y,pos,polwert);

                    //speichere ab
                    result[pos] = (uint8_t) polwert;
                }
            }
        }
    }
}

// interpolate for the standard implementation
void interpolate(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    // grayscale naive
    if (scale_factor == 1) {
        grayscale_naiv(img, width, height, a, b, c, tmp);
        return;
    } 
    grayscale_naiv(img, width, height, a, b, c, tmp);

    // interpolation naive
    bilinear_interpolate_naive(tmp, width, height, scale_factor, result);

}

// interpolate with optimized grayscale (look-up tables)
void interpolate_optimized1(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    // grayscale optimized with look-up tables
    if (scale_factor == 1) {
        grayscale_look_up(img, width, height, a, b, c, tmp);
        return;
    }
    grayscale_look_up(img, width, height, a, b, c, tmp);

    // interpolation naive
    bilinear_interpolate_naive(tmp, width, height, scale_factor, result);
}

// interpolate with optimiyed grayscale (SIMD)


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
        // Filename are not allowed to start or end with '.'
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
    printf("The coefficients a, b, c: %f %f %f\n", *a, *b, *c); // test coefficients

    // Checking scale_factor
    check_scaling(scaling);
    printf("Scaling factor: %i\n", *scaling); // test scale_factor 
        
    // Checking output filename
    outputFileName = check_output(outputFileName);
    printf("Output Filename: %s\n", outputFileName); // test outputfilename
}

// read ppm header - return 1 if the header is false
void skip_comment(FILE *inputFile) {
    // skip possible comments between the values that need to be read for the header
    int ch;
    fscanf(inputFile, " "); // skip withespace
    while ((ch = fgetc(inputFile)) == '#') { // comment begins
      while ((ch = fgetc(inputFile)) != '\n'); // skip comment 
      fscanf(inputFile, " "); // skip withespace
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
    fscanf(inputFile, "%d", width);
    if (*(width) <= 0) {
        // Error because wrong width
        errno = EINVAL;
        perror("Invalue value for width. Width must be greater than 0. For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    skip_comment(inputFile);
    fscanf(inputFile, "%d", height);
    if (*(height) <= 0) {
        // Error because wrong height
        errno = EINVAL;
        perror("Invalid value for height. Height must be greater than 0. For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    skip_comment(inputFile);
    fscanf(inputFile, "%d", maxColorValue); 
    if (*(maxColorValue) > 255) {
        // Error because wrong maxColorValue
        errno = EINVAL;
        perror("Invalue value for max color. A value smaller or equal 255 is expected. For more information please use the option -h | --help");
        fclose(inputFile);
        return 1;
    }
    if (*(maxColorValue) <= 0) {
        // Error because wrong maxColorValue
        errno = EINVAL;
        perror("Invalid value for max color. Max color must be greater than 0. For more information please use the option -h | --help");
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
        perror("Error by allocating memory");
        fclose(inputFile);
        return NULL;
    }
    fgetc(inputFile); // new line character after end of the header
    fread(pixels, sizeof(uint8_t), *(imageSize)*3, inputFile);
    fclose(inputFile);
    uint8_t* check_pixels = pixels;
    for (int i = 0; i < *(imageSize)*3; i++) {
        if (check_pixels[i] > maxColorValue) {
            // Error because a pixels is greater than the maxColorValue
            errno = EINVAL;
            perror("Invalid value for pixels. RGB-Value of a pixel must be smaller than the max color. For more information please use the option -h | --help");
            fclose(inputFile);
            return NULL;
        }
    }
    fclose(inputFile);
    return pixels;
}

// write ppm
int write_ppm(char* outputFileName, int width, int height, int scaling, uint8_t* result) {
    FILE* outputFile = fopen(outputFileName, "wb");
    if (outputFile == NULL) {
        // Error opening the outputfile 
        perror("Error by opening the outputfile. For more information please use the option -h | --help.");
        fclose(outputFile);
        return 1;
    }
    fprintf(outputFile, "P5\n");
    fprintf(outputFile, "%d %d\n", (width*scaling), (height*scaling));
    fprintf(outputFile, "255\n");
    fwrite(result, sizeof(uint8_t), (width*scaling)*(height*scaling)*sizeof(uint8_t), outputFile);
    fclose(outputFile);
    return 0;
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

                // Verwendungsbeispiele fehlen !!!
                printf("Description of all optiions of the program:\n");
                printf("-V<Number> : The implementation, that will be used. For -V 0, the standard implementation will be used. If this option is not called, the standard implementation also will be used.\n");
                printf("-B<Number> : If this option is called, the duration for the specified implementation will be measured. This optional argument indicates the number of the repetitions of the function calling.\n");
                printf("<Filename> : Positional argument for the input filename. This option must be called, either it will return a error because of missing entries.\n");
                printf("-o<Filename> : Output filename If this option is not called, the default value 'output.pgm' will be used.\n");
                printf("--coeffs<FP Number>,<FP Number>,<FP Number> : The coefficients a, b and c for the grayscale. If this option is not called, the default value 0.299, 0.587 und 0.114 will be used.\n");
                printf("-f<Number> : Scaling factor. If this option is not called, the default value 2 will be used.\n");
                printf("-h|--help : A description of all options of the program with usage examples will be displayed.\n");
                exit(0);
            default:
                // Error because of missing entries
                fprintf(stderr, "All necessary parameters are missing or wrong use of parameters. For more information please use the option -h | --help.\n");
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

    // read ppm
    int width; int height; int imageSize;
    const uint8_t* pixels = read_ppm(inputFileName, &width, &height, &imageSize);
    if (pixels == NULL) {
        return 1; // error by reading ppm
    }
        
    // allocating memories 
    uint8_t* temp = (uint8_t*)malloc((imageSize * sizeof(uint8_t)) + (imageSize * scaling * scaling * sizeof(uint8_t)));
    if (temp == NULL) {
        // Error by allocating memory
        perror("Error by allocating memory");
        return 1;
    }
    uint8_t* result = temp + (imageSize * sizeof(uint8_t));
    if (result == NULL) {
        // Error by allocating memory
        perror("Error by allocating memory");
        return 1;
    }
    // result is a part of the allocated memory for temp. The interpolated picture will be saved by result, the grayscale by temp
    
    // Calculation with different implementations
    switch (implementation) {
        case 0:
            // standard implementation
            printf("Standard implementation\n");
            if (benchmark == 1) {
                // measuaring duration
                // Noch nicht implementiert
                if (repetitions < 1) { repetitions = 1; }
                for (int i = 1; i < repetitions; i++) {
                    // call the functions 
                }
                printf("Times of repetitions: %d \n", repetitions);
                break;
            }
            
            // grayscale and interpolation  
            interpolate(pixels, width, height, a, b, c, scaling, temp, result);
            break;

        default:
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