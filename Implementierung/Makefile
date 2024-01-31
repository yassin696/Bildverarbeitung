CC=gcc
CFLAGS=-Wall -Wextra -std=c11 -msse4.1-O3
LDFLAGS=-lm -pthread

TARGET=interpolate
OBJ=finalfile.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(OBJ)
