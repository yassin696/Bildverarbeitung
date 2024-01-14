#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>

int main(int argc, char **argv){
    //Rahmenprogramm

    // Variable Deklaration
    int implementation = 0;  // Hauptimplementierung
    char* inputFileName = NULL;
    char* outputFileName = "output.ppm";

    // Kommandozeilen-Argumente parsen
    
        
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
            printf("here\n");
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
        int height;
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
        
        printf("Header stimmt vollständig\n");

        int imageSize = width * height;
        uint8_t* pixels = (uint8_t*)malloc(imageSize * 3 * sizeof(uint8_t));
        if (pixels == NULL) { 
            // Fehler beim Allizieren des Speichers
            printf("Fehler beim Allozieren des Speichers.\n");
            fclose(inputFile);
            return 1;
        }
        fgetc(inputFile); // new line character
        fread(pixels, sizeof(uint8_t), imageSize*3, inputFile);
        fclose(inputFile);
        printf("Berechnung startet\n");     
        
            // Abspeichern
            FILE* outputFile = fopen(outputFileName, "wb");
            if (outputFile == NULL) {
                fprintf(stderr, "Ungültiges Ausgabedateiformat. Es wird ein P5 PGM-Bild erwartet.\n");
                fclose(outputFile);
                return 1;
            }
            fprintf(outputFile, "P6\n");
            fprintf(outputFile, "%d %d\n", (width), (height));
            fprintf(outputFile, "255\n");
            fwrite(pixels, sizeof(uint8_t), (width)*(height)*sizeof(uint8_t)*3, outputFile);
            fclose(outputFile);
            
        free(pixels);

    } else {
        // Andere Version ???
    }

    return 0;
}