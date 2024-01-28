
.PHONY: all clean run

CC = gcc
CFLAGS = -Wall -Wextra
LDFLAGS = -lm -lpthread
-std=c11

all: interpolate

interpolate: finalfile.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -f interpolate


