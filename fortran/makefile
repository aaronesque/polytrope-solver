CFLAGS = -fdefault-real-8

solve:
	gfortran ${CFLAGS} solver.f90 main.f90 -o object/main -J object
	./object/main

clean:
	rm ./object/*
