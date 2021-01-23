CC = gcc
CFLAGS = -lgomp
SRC=$(wildcard *.c)


main: $(SRC)
	gcc -o $@ $^ $(CFLAGS)

clean:
	rm main