
.PHONY: all clean run

CC = gcc
CFLAGS = -Wall -Wextra
LDFLAGS = -lm -lpthread


all: interpolate

interpolate: finalfile.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -f interpolate


