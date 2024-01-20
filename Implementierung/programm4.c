#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <errno.h>

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
void gentable(float coeff, float* array) {
   // generate look-up tables
    for (int i = 0; i < 256; ++i) {
        array[i] = (i * coeff);
    }
}

void grayscale_look_up(const uint8_t* img, size_t width, size_t height, float a, float b, float c, uint8_t* tmp) {
    // the tables are smaller than 1 MB, which a modern computer has as the stack for a process
    float* tableA = malloc(256 * sizeof(float));
    float* tableB = malloc(256 * sizeof(float));
    float* tableC = malloc(256 * sizeof(float));
    gentable(a, tableA);
    gentable(b, tableB);
    gentable(c, tableC);
    size_t size = height * width;
    for (size_t i = 0; i < size; ++i) {
        uint8_t R = img[i];
        uint8_t G = img[size + i];
        uint8_t B = img[size + size + i];

        // Lookup precomputed values from tables
        float weightedSum = tableA[R] + tableB[G] + tableC[B];
        tmp[i] =(uint8_t)weightedSum;
    }
    free(tableA);
    free(tableB);
    free(tableC);
}

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
    grayscale_naiv(img, width, height, a, b, c, tmp);

    // interpolation naive
    bilinear_interpolate_naive(tmp, width, height, scale_factor, result);
}


// framework
// checking coefficients and precalculate them (a=a/(a+b+c), b=b/(a+b+c), c=c/(a+b+c))
float* check_coefficients(float a, float b, float c) {
    float sum = a + b + c;
    if (sum <= 0) {
        a = 0.299;
        b = 0.587;
        c = 0.114;
        sum = 1.0;
        printf("The sum of a, b and c is smaller or equal 0. Therefore, the default value will be used.\n");
    }
    if (sum != 1) {
         a = a / sum;
         b = b / sum;
          c = c / sum;
    }
    static float  f[3];
    f[0]=a; f[1]=b; f[2]=c;
    float* ptr = f;
    return ptr;
}

// checking scaling
int check_scaling(int scaling) {
    if (scaling < 1) {
        printf("The scale factor is smaller than 1, so the default value 2 will be used.\n");
        scaling = 2;
    }
    return scaling;
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



// read ppm header - return 1 if the header is false
// skip possible comments between the values that need to be read for the header
void skip_comment(FILE *inputFile) {
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
        // Error because wrong pictureformat 
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
                fprintf(stderr, "All necessary parameters are missing. For more information please use the option -h | --help.\n");
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
    if (implementation == 0) {
        // standard implementation
        printf("Standard implementation\n");

        // Checking and calculating of coefficients a, b, c (a=a/(a+b+c), b=b/(a+b+c), c=c/(a+b+c))
        float* f = check_coefficients(a, b, c);
        a=f[0], b=f[1], c=f[2];
        printf("The coefficients a, b, c: %f %f %f\n", a, b, c); // test coefficients

        // Checking scale_factor
        scaling = check_scaling(scaling);
        printf("Scaling factor: %i\n", scaling);
        
        // Checking output filename
        outputFileName = check_output(outputFileName);
        printf("Output Filename: %s\n", outputFileName);

        // check inputfile 
        if (strcmp(strrchr(inputFileName, '\0') - 4, ".ppm") != 0) {
            // The inputfile does not end with ".ppm" -> false file format
            errno = EIO;
            perror("Error: Wrong file format. A PPM-Picutre is expected. For more information please use the option -h | --help");
            return 1;
        }
        printf("Inputfile: %s\n", inputFileName);
       
        // read the ppm-header
        FILE* inputFile = fopen(inputFileName, "rb");
        if (inputFile == NULL) {
            // Error while trying open the giving picture
            perror("The giving picutre can not be opened. For more information please use the option -h | --help");
            return 1;
        }
        int width, height, maxColorValue;
        if (read_ppm_header(inputFile, &width, &height, &maxColorValue) == 1) {
            return 1;
        }
        printf("Width, height, maxColorValue: %i, %i, %i\n", width, height, maxColorValue);

        int imageSize = width * height;
        uint8_t* pixels = (uint8_t*)malloc(imageSize * 3 * sizeof(uint8_t));
        if (pixels == NULL) {
            // Error by allocating memory
            perror("Error by allocating memory");
            fclose(inputFile);
            return 1;
        }
        
        //fgetc(inputFile); // new line character
        fread(pixels, sizeof(uint8_t), imageSize*3, inputFile);
        fclose(inputFile);
        //printf("Berechnung startet\n");

        // Calculation
        if (benchmark == 1) {
            // measuaring duration
            // Noch nicht implementiert
            if (repetitions < 1) { repetitions = 1; }
            for (int i = 1; i < repetitions; i++) {
                // call the functions 
            }
            printf("Times of repetitions: %d \n", repetitions);
            return 0;
        }
        uint8_t* temp = (uint8_t*)malloc(imageSize * sizeof(uint8_t));
        if (temp == NULL) {
            // Error by allocating memory
            perror("Error by allocating memory");
            return 1;
        }
        uint8_t* result = (uint8_t*)malloc(imageSize * scaling * scaling * sizeof(uint8_t));
        if (result == NULL) {
            // Error by allocating memory
            perror("Error by allocating memory");
            return 1;
        }
        
        // Calculation 
        interpolate(pixels, width, height, a, b, c, scaling, temp, result);

        // Saving the result picture 
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
        free(result);
        free(temp);
        free(pixels);

    } else if (implementation == 1) {
        printf("Grayscale optimiert durch Lookup\n");

        // Berechnung und Prüfung der Koeffizienten a, b, c
        float sum = a + b + c;
        if (sum <= 0) {
            a = 0.299;
            b = 0.587;
            c = 0.114;
            sum = 1.0;
            printf("Die Koeffzienten a, b und c sind in der Summe kleiner gleich 0. Deshalb werden Standardwerte verwendet.\n");
        }
        if (sum != 1) {
            a = a / sum;
            b = b / sum;
            c = c / sum;
        }
        printf("Die Koeffizienten a, b, c: %f %f %f\n", a, b, c); // Testen der Koeffizientenberechnung

        // Prüfung der Skalierungsfaktor
        if (scaling < 1) {
            printf("Skalierungsfaktor kleiner als 1, deshalb wird Skalierungsfaktor auf 1 gesetzt.");
            scaling = 1;
        }
        printf("Skalierungsfaktor: %i\n", scaling);
        // oder doch mit Fehlermeldung abbrechen?

        // Prüfung der Ausgabedateiname
        // Die Zeichen \ / : * ? " < > | . sind in Dateinamen nicht erlaubt. 
        // Quelle: https://verwaltung.uni-koeln.de/stabsstelle01/content/benutzerberatung/it_faq/windows/faqitems163122/index_ger.html
        char* filename = outputFileName;
        int length = strlen(outputFileName);
        if (length == 0 || length > 255) {
            // Dateiname zu lang oder nicht existierend
            printf("Dateiname ist zu lang oder nicht existierend, deshalb wird output.pgm als Ausgabedateiname verwendet.");
            outputFileName = "output.pgm";
        }
        if (filename[0] == '.' || filename[length-1] == '.') {
            // Dateiname darf nicht mit . beginnen oder enden
            printf("Dateiname darf nicht mit . beginnen oder enden, deshalb wird output.pgm als Ausgabedateiname verwendet.");
            outputFileName = "output.pgm";
        }
        int whitespace = 1;
        for (int i = 0; i < length; i++) {
            if (filename[i] != ' ') {
                whitespace = 0;
            }
            if (filename[i] == '/' || filename[i] == '\\' || filename[i] == ':' || filename[i] == '*' || filename[i] == '?' || filename[i] == '"' || filename[i] == '<' || filename[i] == '>' || filename[i] == '|') {
                // Dateiname enthält nicht erlaubte Zeichen
                printf("Dateiname enthaelt nicht erlaubte Zeichen, deshalb wird output.pgm als Ausgabedateiname verwendet.");
                outputFileName = "output.pgm";
            }
        }
        if (whitespace == 1) {
            // Dateiname enthält nur Leerzeichen
            printf("Dateiname enthaelt nur Leerzeichen, deshalb wird output.pgm als Ausgabedateiname verwendet.");
            outputFileName = "output.pgm";
        }
        printf("Ausgabedatei: %s\n", outputFileName);
        // oder doch mit Fehlermeldung abbrechen?

        // Prüfen des Bildes 
        if (strcmp(strrchr(inputFileName, '\0') - 4, ".ppm") != 0) {
            // Die Eingabedatei endet nicht mit  ".ppm" -> falsches Eingabeformat
            printf("Test failed here \n");
            fprintf(stderr, "Ungueltiges Dateiformat. Es wird ein PPM-Bild erwartet.\n");
            return 1;
        }
        printf("Eingabedatei: %s\n", inputFileName);

        FILE* inputFile = fopen(inputFileName, "rb");
        if (inputFile == NULL) {
            // Fehler beim Öffnen des Bildes
            fprintf(stderr, "Das angegebene Bild kann nicht geoeffnet werden.");
        }
        //printf("Bild kann geoeffnet werden!\n");
        char magicNumber[3];
        fscanf(inputFile, "%s", magicNumber);
        if (strcmp(magicNumber, "P6") != 0) {
            // Fehlermeldung wegen falsches Bildformat
            printf("here\n");
            fprintf(stderr, "Ungueltiges Dateiformat. Es wird ein P6 PPM-Bild erwartet.\n");
            fclose(inputFile);
            return 1;
        }
        //printf("Header erste Zeile stimmt\n");
        int ch;
        while(1) {
            while((ch=fgetc(inputFile)) == ' '); // whitespace überspringen
            if (ch == '#') {
                while((ch=fgetc(inputFile)) != '\n'); // Kommentar überspringen
            } else {
                ungetc(ch, inputFile); //letzte gelesene Zeichen zurückgehen
                break;
            }
        }
        int width;
        fscanf(inputFile, "%d", &width);
        while(1) {
            while((ch=fgetc(inputFile)) == ' '); // whitespace überspringen
            if (ch == '#') {
                while((ch=fgetc(inputFile)) != '\n'); // Kommentar überspringen
            } else {
                ungetc(ch, inputFile); //letzte gelesene Zeichen zurückgehen
                break;
            }
        }
        int height;
        fscanf(inputFile, "%d", &height);
        while(1) {
            while((ch=fgetc(inputFile)) == ' '); // whitespace überspringen
            if (ch == '#') {
                while((ch=fgetc(inputFile)) != '\n'); // Kommentar überspringen
            } else {
                ungetc(ch, inputFile); //letzte gelesene Zeichen zurückgehen
                break;
            }
        }
        int maxColorValue;
        fscanf(inputFile, "%d", &maxColorValue);
        printf("Width, height, maxColorValue: %i, %i, %i\n", width, height, maxColorValue);
        if (maxColorValue > 255) {
            // Fehlermeldung wegen falsches Bildformat
            fprintf(stderr, "Ungueltiger Maximalwert für die Farben. Es wird ein Wert von 255 erwartet.\n");
            fclose(inputFile);
            return 1;
        }
        if (maxColorValue == 0) {
            // Fehlermeldung wegen falsches Bildformat
            fprintf(stderr, "Ungueltiger Maximalwert für die Farben. Null als Wer ist nicht erlaubt.\n");
            fclose(inputFile);
            return 1;
        }
        //printf("Header stimmt vollständig\n");
        int imageSize = width * height;
        typedef struct {
            uint8_t r;
            uint8_t g;
            uint8_t b;
        } Pixel;
        Pixel* pixels = (Pixel*)malloc(imageSize * sizeof(Pixel));
        fgetc(inputFile); // new line character
        fread(pixels, sizeof(Pixel), imageSize, inputFile);
        fclose(inputFile);
        //printf("Berechnung startet\n");

        // Berechnung
        if (benchmark == 1) {
            // Laufzeit wird gemessen und ausgegeben
            // Noch nicht implementiert
            if (repetitions < 1) { repetitions = 1; }
            // oder doch mit Fehlermeldung abbrechen?
            for (int i = 1; i < repetitions; i++) {
                // Programm ausführen
            }
            printf("Anzahl der Wiederholungen: %d \n", repetitions);
            return 0;
        }
        uint8_t* temp = (uint8_t*)malloc(imageSize * sizeof(uint8_t));
        if (temp == NULL) {
            // Fehler beim Allizieren des Speichers
            printf("Fehler beim Allozieren des Speichers.\n");
            return 1;
        }
        uint8_t* result = (uint8_t*)malloc(imageSize * scaling * scaling * sizeof(uint8_t));
        if (result == NULL) {
            // Fehler beim Allizieren des Speichers
            printf("Fehler beim Allozieren des Speichers.\n");
            return 1;
        }
        uint8_t* img = (uint8_t*)malloc(imageSize * 3 * sizeof(uint8_t));
        if (img == NULL) {
            // Fehler beim Allizieren des Speichers
            printf("Fehler beim Allozieren des Speichers.\n");
            return 1;
        }
        // img speichert zuerst alle r-value, dann alle g-value und zuletzt alle b-value
        for (int i = 0; i < imageSize; i++) {
            img[i] = (uint8_t) pixels[i].r;
            img[i+imageSize] = (uint8_t) pixels[i].g;
            img[i+imageSize+imageSize] = (uint8_t) pixels[i].b;
        }
        
        // Berechnung
        interpolate(img, width, height, a, b, c, scaling, temp, result);

        // Abspeichern
        FILE* outputFile = fopen(outputFileName, "wb");
        if (outputFile == NULL) {
            fprintf(stderr, "Ungueltiges Ausgabedateiformat. Es wird ein P5 PGM-Bild erwartet.\n");
            fclose(outputFile);
            return 1;
        }
        fprintf(outputFile, "P5\n");
        fprintf(outputFile, "%d %d\n", (width*scaling), (height*scaling));
        fprintf(outputFile, "255\n");
        fwrite(result, sizeof(uint8_t), (width*scaling)*(height*scaling)*sizeof(uint8_t), outputFile);
        fclose(outputFile);
        free(result);
        free(temp);
        free(pixels);
    } else {
        //another implementation ??? 
    }

    return 0;
}