CC = gcc

rho: main.c rho.o
	${CC} -std=c99 -O0 -g -o rho main.c rho.o -lm

rho-o: rho.c rho.h
	${CC} -c -Wall -Werror -fpic rho.c rho.h
	${CC} -shared -o librho.so rho.o

clean:
	@rm -fr *.o *.so rho
	@echo "==================="
	@echo "clean complete !"
	@echo "==================="
