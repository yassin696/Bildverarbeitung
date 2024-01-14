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
    int j = 0;
    while (i < height * width) {
        uint8_t R = img[i];
        uint8_t G = img[i + width*height];
        uint8_t B = img[i + width*height+width*height];
        uint8_t D = round(R * a + G * b + B * c);
        temp[j] = D;
        i += 1;
        j += 1;
    }
}
void interpolation(uint8_t* intArray, uint8_t* allozspeicher, int breite, int hoehe, int factor) {
    // Interpolation
    
    // Einzelne Sektoren werden bearbeitet
    for(int sektorh = 0; sektorh < hoehe;sektorh++) {
        for (int sektorb = 0; sektorb < breite; sektorb++) {
            //Sektoren a,b,c,d Werte, also (0,0), (0,s), (s,0) und (s,s) werden für alle berechnet
            int a = intArray[sektorh * breite + sektorb];
            int b = intArray[sektorh * breite + (sektorb + 1) % breite];
            int c = intArray[((sektorh + 1) % hoehe) * breite + sektorb];
            int d = intArray[((sektorh + 1) % hoehe) * breite + (sektorb + 1) % breite];

            // Für jeden einzelnen Wert in den Quadraten berechne nach Formel den Wert
            for(int y=0;y<factor;y++) {
                for (int x = 0; x < factor; x++) {


                    // berechne position im finalen AllozSpeicher
                    int pos = (x + factor * sektorb) + ((sektorh * factor + y) * breite * factor);

                    // da wert x=0, y=0 immer den a-Wert ergibt kann dieser auch direkt eingetragen werden
                    // die berechnung füllt ja nur das quadrat für eins kleiner von s aus, dass keine Werte doppelt berechnet werden
                    if(x == 0 && y == 0) {
                        allozspeicher[pos] = (uint8_t) a;
                        continue;
                    }

                    // berechne Wert
                    int polwert = (a * (factor-y) * (factor-x) ) + ( c * y * ( factor-x) ) + ( b * (factor-y) * x ) + ( d * y * x );
                    // multipliziere mit (1 / s*s)
                    polwert = polwert / (factor * factor);
                    //printf("[%i;%i] %i - %i\n",x,y,pos,polwert);

                    //speichere ab
                    allozspeicher[pos] = (uint8_t) polwert;
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
    int* option_index = 0;
    while ((option = getopt_long(argc, argv, "VBof:h", long_options, option_index)) != -1) {
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
            printf("-V<Zahl> : Die Implementierung, die verwendet werden soll. Hierbei soll mit -V 0 Ihre Hauptimplementierung verwendet werden. Wenn diese Option nicht gesetzt wird, soll ebenfalls die Hauptimplementierung ausgeführt werden.\n");
            printf("-B<Zahl> : Falls gesetzt, wird die Laufzeit der angegebenen Implementierung gemessen und ausgegeben. Das optionale Argument dieser Option gibt die Anzahl an Wiederholungen des Funktionsaufrufs an.\n");
            printf("<Dateiname> : Positionales Argument für die Eingabedatei.\n");
            printf("-o<Dateiname> : Ausgabedatei.\n");
            printf("--coeffs<FP Zahl>,<FP Zahl>,<FP Zahl> : Die Koeffizienten der Graustufenkonvertierung a, b und c. Falls diese Option nicht gesetzt wird, werden die Standardwerte 0.299, 0.587 und 0.114 verwendet.\n");
            printf("-f<Zahl> : Skalierungsfaktor.\n");
            printf("-h|--help : Eine Beschreibung aller Optionen des Programms und Verwendungsbeispiele.\n");
            exit(0);            
        default:
            // Fehlermehldung wegen fehlender Eingaben
            fprintf(stderr, "Es fehlen alle nötigen Parameter. Für weitere Informationen verwenden Sie bitte die Option -h|--help.");
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
        if (sum != 1) {
            a = a / sum;
            b = b / sum;
            c = c / sum;
        }
        printf("Die Koeffizienten a, b, c: %f %f %f\n", a, b, c); // Testen der Koeffizientenberechnung

        // Prüfung der Skalierungsfaktor
        if (scaling < 1) { scaling = 1; }
        printf("Skalierungsfaktor: %i\n", scaling);
        // oder doch mit Fehlermeldung abbrechen?

        // Prüfung der Ausgabedateiname
        // Die Zeichen \ / : * ? " < > | . sind in Dateinamen nicht erlaubt. 
        // Quelle: https://verwaltung.uni-koeln.de/stabsstelle01/content/benutzerberatung/it_faq/windows/faqitems163122/index_ger.html
        char* filename = outputFileName;
        int length = strlen(outputFileName);
        if (length == 0 || length > 255) {
            // Dateiname zu lang oder nicht existierend
            outputFileName = "output.pgm";
        }
        if (filename[0] == '.' || filename[length-1] == '.') {
            // Dateiname darf nicht mit . beginnen oder enden
            outputFileName = "output.pgm";
        }
        int whitespace = 1;
        for (int i = 0; i < length; i++) {
            if (filename[i] != ' ') {
                whitespace = 0;
            }
            if (filename[i] == '/' || filename[i] == '\\' || filename[i] == ':' || filename[i] == '*' || filename[i] == '?' || filename[i] == '"' || filename[i] == '<' || filename[i] == '>' || filename[i] == '|') {
                // Dateiname enthält nicht erlaubte Zeichen
                outputFileName = "output.pgm";
            }
        }
        if (whitespace == 1) {
            // Dateiname enthält nur Leerzeichen
            outputFileName = "output.pgm";
        }
        printf("Ausgabedatei: %s\n", outputFileName);
        // oder doch mit Fehlermeldung abbrechen?

        // Prüfen des Bildes 
        if (strcmp(strrchr(inputFileName, '\0') - 4, ".ppm") != 0) {
            // Die Eingabedatei endet nicht mit  ".ppm" -> falsches Eingabeformat
            printf("Test failed here \n");
            fprintf(stderr, "Ungültiges Dateiformat. Es wird ein PPM-Bild erwartet.\n");
            return 1;
        }
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
        printf("Header erste Zeile stimmt\n");
        int width = 0;
        int height = 0;
        int maxColorValue;
        int ch; 
        while ((fgetc(inputFile)) == '#') {
            while ((ch = fgetc(inputFile)) != '\n' && ch != EOF); // Alle Kommentare ignorieren
        }
        fseek(inputFile, -1, SEEK_CUR);
        fscanf(inputFile, "%d %d", &width, &height);
        fscanf(inputFile, "%d", &maxColorValue);
        if (maxColorValue > 255) {
            // Fehlermeldung wegen falsches Bildformat
            fprintf(stderr, "Ungültiger Maximalwert für die Farben. Es wird ein Wert von 255 erwartet.\n");
            fclose(inputFile);
            return 1;
        }
        printf("Header stimmt vollständig\n");
        int imageSize = width * height;
        typedef struct {
            uint8_t r;
            uint8_t g;
            uint8_t b;
        } Pixel;
        Pixel* pixels = (Pixel*)malloc(imageSize * 3 * sizeof(Pixel));
        if (pixels == NULL) {
            // Fehler beim Allizieren des Speichers
            printf("Fehler beim Allozieren des Speichers.\n");
            fclose(inputFile);
            return 1;
        }
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
        uint8_t* img = (uint8_t*)malloc(imageSize * 3 * sizeof(uint8_t));
        if (img == NULL) {
            // Fehler beim Allizieren des Speichers
            printf("Fehler beim Allozieren des Speichers.\n");
            return 1;
        }
        // img speichert zuerst alle r-value, dann alle g-value und zuletzt alle b-value
        for (int i = 0; i < imageSize; i++) {
            img[i] = pixels[i].r;
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