#include <stdbool.h>
#include "stdio.h"
#include <stdint.h>
#include <stdlib.h>

float* grayscale(const uint8_t* img, size_t width, size_t height, float a, float b, float c) {
    float* result = malloc((width ) * (height ) * sizeof(int));

    if (result == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < height * width*3; i += 3) {
        uint8_t R = img[i];
        uint8_t G = img[i + 1];
        uint8_t B = img[i + 2];

       
        float gray = (R * a + G * b + B * c)/(a+b+c);
        result[i / 3] = gray;
    }
return result;
}

bool intpolcalculation(float* intArray, float* allozspeicher, int hoehe, int breite, int factor) {
    // Einzelne Sektoren werden bearbeitet
    for(int sektorh = 0; sektorh < hoehe;sektorh++) {
        for (int sektorb = 0; sektorb < breite; sektorb++) {
            //Sektoren a,b,c,d Werte, also (0,0), (0,s), (s,0) und (s,s) werden für alle berechnet
            float a = intArray[sektorh * breite + sektorb];
            float b = intArray[sektorh * breite + (sektorb + 1) % breite];
            float c = intArray[((sektorh + 1) % hoehe) * breite + sektorb];
            float d = intArray[((sektorh + 1) % hoehe) * breite + (sektorb + 1) % breite];

            // Für jeden einzelnen Wert in den Quadraten berechne nach Formel den Wert
            for(int y=0;y<factor;y++) {
                for (int x = 0; x < factor; x++) {


                    // berechne position im finalen AllozSpeicher
                    int pos = (x + factor * sektorb) + ((sektorh * factor + y) * breite * factor);

                    // da wert x=0, y=0 immer den a-Wert ergibt kann dieser auch direkt eingetragen werden
                    // die berechnung füllt ja nur das quadrat für eins kleiner von s aus, dass keine Werte doppelt berechnet werden
                    if(x == 0 && y == 0) {
                        allozspeicher[pos] = a;
                        continue;
                    }

                    // berechne Wert
                    float polwert = (a * (factor-y) * (factor-x) ) + ( c * y * ( factor-x) ) + ( b * (factor-y) * x ) + ( d * y * x );
                    // multipliziere mit (1 / s*s)
                    polwert = polwert / (factor * factor);
                    //printf("[%i;%i] %i - %i\n",x,y,pos,polwert);

                    //speichere ab
                    allozspeicher[pos] = polwert;
                }
            }
        }
    }
    return true;
}

int main(){
    size_t width = 1;
    size_t height = 3;
uint8_t inputImage[] = {
    255, 0, 0,0, 255, 0,0, 0, 255 
};
    // Test grayscale conversion
    float* result = grayscale(inputImage, width, height, 0.299, 0.587, 0.114);

    if (result == NULL) {
        // Handle memory allocation failure
        return 1;
    }

    // Print the resulting grayscale values
    printf("Grayscale values:\n");
    for (size_t i = 0; i < width*height; ++i) {
        printf("pixel number %d %f \n",i, result[i]);
    }
    printf("\n");


    // Don't forget to free the allocated memory
 
int factor = 3;
float allozspeicher[height * width * factor * factor];
intpolcalculation(result,allozspeicher,height,width,factor);
    for (int i = 0; i < height * width * factor * factor; ++i) {
        if(i % (width * factor) == 0) {
            printf("\n");
        }
        printf("%d ", allozspeicher[i]);
    }
   free(result);
    return 0;
}
/*void main() {
    int hoehe = 2;
    int breite = 4;
    int factor = 3;
    int intArray[8] = {34, 67, 123, 89, 45, 210, 156, 78};
    int allozspeicher[height * width * factor * factor];
    //for (int i = 0; i < hoehe * breite * factor * factor; ++i) {
        // -1 als Fehlerwert, dass mans schöner sehen kann
        //allozspeicher[i] = -1;
    //}
    intpolcalculation(intArray,allozspeicher,hoehe,breite,factor);
    for (int i = 0; i < hoehe * breite * factor * factor; ++i) {
        if(i % (breite * factor) == 0) {
            printf("\n");
        }
        printf("%d ", allozspeicher[i]);
    }
}*/