CC=g++
CFLAGS=`pkg-config --cflags eigen3 opencv4`
LFLAGS=`pkg-config --libs eigen3 opencv4` -pthread
DEPS=

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hellomake: main.o
	$(CC) -o mandelbrot main.o $(LFLAGS)

