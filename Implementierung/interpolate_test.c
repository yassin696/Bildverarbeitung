//
// Created by Constantin Carste on 27.12.23.
//

#include <stdbool.h>
#include "stdio.h"
#include "inttypes.h"

void intpolcalculation(uint8_t* intArray, uint8_t* allozspeicher, int breite, int hoehe, int factor) {
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

void main() {
    int hoehe = 2;
    int breite = 4;
    int factor = 3;
    u_int8_t intArray[8] = {1, 75, 200, 75, 45, 210, 156, 78};
    u_int8_t allozspeicher[hoehe * breite * factor * factor];
    //for (int i = 0; i < hoehe * breite * factor * factor; ++i) {
        // -1 als Fehlerwert, dass mans schöner sehen kann
        //allozspeicher[i] = -1;
    //}

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