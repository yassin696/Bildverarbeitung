//
// Created by Constantin Carste on 20.01.24.
//

//
// Created by Constantin Carste on 27.12.23.
//

#include <stdbool.h>
#include "stdio.h"
#include "inttypes.h"

void intpolcalculation(uint8_t* intArray, uint8_t* allozspeicher, int breite, int hoehe, int factor) {
    // Berechne die einzelnen horizontalen Linien vom neuen Bild

    //Berechne die vertikalen Werte
    for(int h=0;h<hoehe;h++) {
        for (int b = 0; b < breite; b++) {
            float v1_1 = (float) intArray[b + (h * breite)];
            float v1_2 = (float) intArray[((b + 1) % breite) + (h * breite)];
            float v2_1 = (float) intArray[b + ((h+1) % hoehe * breite)];
            float v2_2 = (float) intArray[((b + 1) % breite) + ((h+1) % hoehe * breite)];
            for (int f1 = 0; f1 < factor; f1++) {
                float val1 = ((float) (factor - f1) / (float) factor) * v1_1 + ((float) f1 / (float) factor) * v1_2;
                float val2 = ((float) (factor - f1) / (float) factor) * v2_1 + ((float) f1 / (float) factor) * v2_2;
                for (int f = 0; f < factor; f++) {
                    allozspeicher[(h * factor + f) * breite * factor + (b * factor + f1) ] = (uint8_t) (
                            ((float) (factor - f) / (float) factor) * val1 +
                            ((float) f / (float) factor) * val2);
                }
            }
        }
    }
}

void main() {
    int hoehe = 2;
    int breite = 4;
    int factor = 6;
    u_int8_t intArray[8] = {1, 75, 200, 75, 45, 210, 156, 78};
    u_int8_t allozspeicher[hoehe * breite * factor * factor];
    for (int i = 0; i < hoehe * breite * factor * factor; ++i) {
    // -1 als Fehlerwert, dass mans schöner sehen kann
        allozspeicher[i] = 0;
    }

    /*
     *    1  75  200 75
     *    45 210 156 78
     */


    intpolcalculation(intArray,allozspeicher,hoehe,breite,factor);
    for (int i = 0; i < hoehe * breite * factor * factor; ++i) {
        if(i % (breite * factor) == 0) {
            printf("\n");
        }
        printf("%d ", allozspeicher[i]);
    }
}