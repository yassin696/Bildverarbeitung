
.PHONY: all clean run

CC = gcc
CFLAGS = -Wall -Wextra
LDFLAGS = -lm -lpthread
-std=c11

all: interpolate

interpolate: programm6.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -f interpolate

run: interpolate
	./interpolate
