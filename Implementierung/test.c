#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>


void grayscale(uint8_t* img, uint8_t* temp, int width, int height, float a, float b, float c) {
    // Graustufenkonvertierung
    printf("grayscale\n");
    
    int i = 0;
    while (i < height * width) {
        uint8_t R = img[i];
        uint8_t G = img[i + width*height];
        uint8_t B = img[i + width*height*2];
        uint8_t D = round(R * a + G * b + B * c);
        temp[i] = D;
        i += 1;
    }
}



int main(int argc, char **argv){
    //Rahmenprogramm

    // Variable Deklaration
    int option;
    int implementation = 0;  // Hauptimplementierung
    int benchmark = 0;// 1, falls Laufzeit gemessen werden soll
    int repetitions = 1; // Anzahl der Wiederholungen
    char* inputFileName = NULL;
    char* outputFileName = "output.pgm";
    float a = 0.299;
    float b = 0.587;
    float c = 0.114;
    int scaling = 1; // Soll die Skalierungsfaktor einen Integer sein?  
        
    if (optind < argc) {
        inputFileName = argv[optind];
    } else {
        // Fehlermeldung wegen Fehlen des Bildes
        fprintf(stderr, "Das Bild, das interpoliert werden soll, ist nicht eingegeben.");
        return 1;
    }

     //Eigentliches Programm
    if (implementation == 0) {
        // Hauptimplementation
        printf("Hauptimplementation\n");
        
        printf("Eingabedatei: %s\n", inputFileName);
    
        FILE* inputFile = fopen(inputFileName, "rb");
        if (inputFile == NULL) {
            // Fehler beim Öffnen des Bildes
            fprintf(stderr, "Das angegebene Bild kann nicht geöffnet werden.");
        }
        printf("Bild kann geöffnet werden!\n");
        char magicNumber[3];
        fscanf(inputFile, "%s", magicNumber);
        if (strcmp(magicNumber, "P6") != 0) {
            // Fehlermeldung wegen falsches Bildformat
            fprintf(stderr, "Ungültiges Dateiformat. Es wird ein P6 PPM-Bild erwartet.\n");
            fclose(inputFile);
            return 1;
        }
        int c;
        while(1) {
            while((c=fgetc(inputFile)) == ' '); // whitespace überspringen
            if (c == '#') {
                while((c=fgetc(inputFile)) != '\n'); // Kommentar überspringen
            } else {
                ungetc(c, inputFile); //letzte gelesene Zeichen zurückgehen
                break;
            }
        }
        int width = 0;
        fscanf(inputFile, "%d", &width);
        while(1) {
            while((c=fgetc(inputFile)) == ' '); // whitespace überspringen
            if (c == '#') {
                while((c=fgetc(inputFile)) != '\n'); // Kommentar überspringen
            } else {
                ungetc(c, inputFile); //letzte gelesene Zeichen zurückgehen
                break;
            }
        }
        int height = 0;
        fscanf(inputFile, "%d", &height);
        while(1) {
            while((c=fgetc(inputFile)) == ' '); // whitespace überspringen
            if (c == '#') {
                while((c=fgetc(inputFile)) != '\n'); // Kommentar überspringen
            } else {
                ungetc(c, inputFile); //letzte gelesene Zeichen zurückgehen
                break;
            }
        }
        int maxColorValue;
        fscanf(inputFile, "%d", &maxColorValue);
        printf("Width, height, maxColorValue: %i, %i, %i\n", width, height, maxColorValue);
        

        int imageSize = width * height;
        typedef struct {
            uint8_t r;
            uint8_t g;
            uint8_t b;
        } Pixel;
        Pixel* pixels = (Pixel*)malloc(imageSize * 3 * sizeof(Pixel));
        //fgetc(inputFile); // new line character
        fread(pixels, sizeof(Pixel), imageSize*3, inputFile);
        fclose(inputFile);
        printf("Berechnung startet\n");
        uint8_t* img = (uint8_t*)malloc(imageSize * 3 * sizeof(uint8_t));
        if (img == NULL) {
            // Fehler beim Allizieren des Speichers
            printf("Fehler beim Allozieren des Speichers.\n");
            return 1;
        }
        // img speichert zuerst alle r-value, dann alle g-value und zuletzt alle b-value
        for (int i = 0; i < imageSize; i++) {
            img[i] = (uint8_t) pixels[i].r;
            img[i+imageSize] = pixels[i].g;
            img[i+imageSize+imageSize] = pixels[i].b;
        }

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
        uint8_t* result = (uint8_t*)malloc(imageSize * scaling * scaling);
        if (result == NULL) {
            // Fehler beim Allizieren des Speichers
            printf("Fehler beim Allozieren des Speichers.\n");
            return 1;
        }
        grayscale(img, temp, width, height, a, b, c);
        printf("grayscale erfolgreich\n");
        if (scaling == 1) {
            FILE* outputFile = fopen(outputFileName, "wb");
            if (outputFile == NULL) {
                fprintf(stderr, "Ungültiges Ausgabedateiformat. Es wird ein P5 PGM-Bild erwartet.\n");
                fclose(outputFile);
                return 1;
            }
            fprintf(outputFile, "P5\n");
            fprintf(outputFile, "%d %d\n", (width*scaling), (height*scaling));
            fprintf(outputFile, "255\n");
            fwrite(temp, sizeof(uint8_t), (width*scaling)*(height*scaling)*sizeof(uint8_t), outputFile);
            fclose(outputFile);
        } else {
            // Abspeichern
            FILE* outputFile = fopen(outputFileName, "wb");
            if (outputFile == NULL) {
                fprintf(stderr, "Ungültiges Ausgabedateiformat. Es wird ein P5 PGM-Bild erwartet.\n");
                fclose(outputFile);
                return 1;
            }
            fprintf(outputFile, "P5\n");
            fprintf(outputFile, "%d %d\n", (width*scaling), (height*scaling));
            fprintf(outputFile, "255\n");
            fwrite(result, sizeof(uint8_t), (width*scaling)*(height*scaling)*sizeof(uint8_t), outputFile);
            fclose(outputFile);
        }
        free(result);
        free(temp);
        free(pixels);

    } else {
        // Andere Version ???
    }

    return 0;
}