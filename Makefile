all: 
	gcc test.c -o test -lcblas -latlas 
clean:
	rm test

