CC=gcc
CFLAGS=-Wall -Wextra -std=c11 -msse4.1 -O3
LDFLAGS=-lm
TARGET=interpolate
OBJ=implementierung.o rahmenprogramm.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET) $(LDFLAGS)

implementierung.o: implementierung.c implementierung.h
	$(CC) $(CFLAGS) -c implementierung.c -o implementierung.o

rahmenprogramm.o: rahmenprogramm.c implementierung.h
	$(CC) $(CFLAGS) -c rahmenprogramm.c -o rahmenprogramm.o

clean:
	rm -f $(TARGET) $(OBJ)
