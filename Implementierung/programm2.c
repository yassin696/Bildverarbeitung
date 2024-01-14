#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>

void interpolate(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    // Graustufenkonvertierung
    size_t i = 0;
    int j = 0;
    if (scale_factor == 1) {
        while (i < height * width * 3) {
            uint8_t R = img[i];
            uint8_t G = img[i + 1];
            uint8_t B = img[i + 2];
            uint8_t D = round(R * a + G * b + B * c);
            result[j] = D;
            i += 3;
            j += 1;
        }
        return;
    } 
    while (i < height * width * 3) {
        uint8_t R = img[i];
        uint8_t G = img[i + 1];
        uint8_t B = img[i + 2];
        uint8_t D = round(R * a + G * b + B * c);
        tmp[j] = D;
        i += 3;
        j += 1;
    }

    // Interpolation
     // Einzelne Sektoren werden bearbeitet
    for(size_t sektorh = 0; sektorh < height;sektorh++) {
        for (size_t sektorb = 0; sektorb < width; sektorb++) {
            //Sektoren a,b,c,d Werte, also (0,0), (0,s), (s,0) und (s,s) werden für alle berechnet
            int a = tmp[sektorh * width + sektorb];
            int b = tmp[sektorh * width + (sektorb + 1) % width];
            int c = tmp[((sektorh + 1) % height) * width + sektorb];
            int d = tmp[((sektorh + 1) % height) * width + (sektorb + 1) % width];

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
                    int polwert = (a * (scale_factor-y) * (scale_factor-x) ) + ( c * y * ( scale_factor-x) ) + ( b * (scale_factor-y) * x ) + ( d * y * x );
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


void gentable(float coeff, int16_t* table) {
    for (int i = 0; i < 256; ++i) {
        table[i] = (int16_t)(i * coeff);
    }
}

void interpolate1(const uint8_t* img, size_t width, size_t height, float a, float b, float c, size_t scale_factor, uint8_t* tmp, uint8_t* result) {
    // Graustufenkonvertierung
    int16_t tableA[256];
    int16_t tableB[256];
    int16_t tableC[256];
    gentable(a, tableA);
    gentable(b, tableB);
    gentable(c, tableC);
    if (scale_factor == 1) {
        for (size_t i = 0; i < height * width; ++i) {
            uint8_t R = img[i * 3];
            uint8_t G = img[i * 3 + 1];
            uint8_t B = img[i * 3 + 2];
            
            // Lookup precomputed values from tables
            int16_t weightedSum = tableA[R] + tableB[G] + tableC[B];
            result[i] = (uint8_t)(weightedSum);
        }
        return;
    }
    for (size_t i = 0; i < height * width; ++i) {
        uint8_t R = img[i * 3];
        uint8_t G = img[i * 3 + 1];
        uint8_t B = img[i * 3 + 2];
            
        // Lookup precomputed values from tables
        int16_t weightedSum = tableA[R] + tableB[G] + tableC[B];
        tmp[i] = (uint8_t)(weightedSum);
    }

    // Interpolation

    // Einzelne Sektoren werden bearbeitet
    for(size_t sektorh = 0; sektorh < height;sektorh++) {
        for (size_t sektorb = 0; sektorb < width; sektorb++) {
            //Sektoren a,b,c,d Werte, also (0,0), (0,s), (s,0) und (s,s) werden für alle berechnet
            int a = tmp[sektorh * width + sektorb];
            int b = tmp[sektorh * width + (sektorb + 1) % width];
            int c = tmp[((sektorh + 1) % height) * width + sektorb];
            int d = tmp[((sektorh + 1) % height) * width + (sektorb + 1) % width];

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
                    int polwert = (a * (scale_factor-y) * (scale_factor-x) ) + ( c * y * ( scale_factor-x) ) + ( b * (scale_factor-y) * x ) + ( d * y * x );
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

    // Kommandozeilen-Argumente parsen
    static struct option long_options[] = {
            {"coeffs", required_argument, 0, 'c'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
    };
    while ((option = getopt_long(argc, argv, "VBof:h", long_options, NULL)) != -1) {
        switch (option) {
            case 0:
                sscanf(optarg, "%f,%f,%f", &a, &b, &c);
                break;
            case 'V':
                implementation = atoi(optarg);
                printf("Implementierung: %i\n", implementation);
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
                //printf("a, b, c: %f %f %f\n", a, b, c);
                break;
            case 'f':
                // Skalierungsfaktor
                scaling = atoi(optarg);
                break;
            case 'h':
                // Beschreibung aller Optionen des Programms und Verwendungsbeispiele werden ausgegeben und das Programm danach beendet
                // Verwendungsbeispiele fehlen !!!
                printf("Beschreibung der Optionen:\n");
                printf("-V<Zahl> : Die Implementierung, die verwendet werden soll. Hierbei wird mit -V 0 die Hauptimplementierung verwendet. Wenn diese Option nicht gesetzt wird, wird ebenfalls die Hauptimplementierung ausgefuehrt werden.\n");
                printf("-B<Zahl> : Falls gesetzt, wird die Laufzeit der angegebenen Implementierung gemessen und ausgegeben. Das optionale Argument dieser Option gibt die Anzahl an Wiederholungen des Funktionsaufrufs an.\n");
                printf("<Dateiname> : Positionales Argument für die Eingabedatei.\n");
                printf("-o<Dateiname> : Ausgabedatei.\n");
                printf("--coeffs<FP Zahl>,<FP Zahl>,<FP Zahl> : Die Koeffizienten der Graustufenkonvertierung a, b und c. Falls diese Option nicht gesetzt wird, werden die Standardwerte 0.299, 0.587 und 0.114 verwendet.\n");
                printf("-f<Zahl> : Skalierungsfaktor.\n");
                printf("-h|--help : Eine Beschreibung aller Optionen des Programms und Verwendungsbeispiele.\n");
                exit(0);
            default:
                // Fehlermehldung wegen fehlender Eingaben
                fprintf(stderr, "Es fehlen alle noetigen Parameter. Fuer weitere Informationen verwenden Sie bitte die Option -h|--help.");
                exit(1);
        }
    }


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
        
        // Berechnung
        interpolate(pixels, width, height, a, b, c, scaling, temp, result);

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
        Pixel* pixels = (Pixel*)malloc(imageSize * 3 * sizeof(Pixel));
        fgetc(inputFile); // new line character
        fread(pixels, sizeof(Pixel), imageSize*3, inputFile);
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
            img[i+imageSize] = pixels[i].g;
            img[i+imageSize+imageSize] = pixels[i].b;
        }
        
        // Berechnung
        interpolate1(img, width, height, a, b, c, scaling, temp, result);

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
        // Andere Implementierung ???
    }

    return 0;
}