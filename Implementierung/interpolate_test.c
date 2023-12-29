//
// Created by Constantin Carste on 27.12.23.
//

#include <stdbool.h>
#include "stdio.h"

bool intpolcalculation(int* intArray, int* allozspeicher, int hoehe, int breite, int factor) {
    for(int sektorh = 0; sektorh < hoehe;sektorh++) {
        for (int sektorb = 0; sektorb < breite; sektorb++) {
            // sektor
            //  a (0,0) -- c (s,0)
            //  |              |
            //  b (s,0) -- d (s,s)

            int a = intArray[sektorh * breite + sektorb];
            int b = intArray[sektorh * breite + (sektorb + 1) % breite];
            int c = intArray[((sektorh + 1) % hoehe) * breite + sektorb];
            int d = intArray[((sektorh + 1) % hoehe) * breite + (sektorb + 1) % breite];

            for(int y=0;y<factor;y++) {
                for (int x = 0; x < factor; x++) {
                    //         -            h o e h e         -
                    //printf("(%i * %i + %i ) * %i\n",sektorh,factor,y,breite);
                    int pos = (x + factor * sektorb) + ((sektorh * factor + y) * breite * factor);
                    int polwert = (a * (factor-y) * (factor-x) ) + ( c * y * ( factor-x) ) + ( b * (factor-y) * x ) + ( d * y * x );
                    polwert = polwert / (factor * factor);
                    //printf("[%i;%i] %i - %i\n",x,y,pos,polwert);
                    allozspeicher[pos] = polwert;
                }
            }
            //return true;
        }
    }
    return true;
}

void main() {
    int hoehe = 2;
    int breite = 4;
    int factor = 3;
    int intArray[8] = {34, 67, 123, 89, 45, 210, 156, 78};
    int allozspeicher[hoehe * breite * factor * factor];
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
}