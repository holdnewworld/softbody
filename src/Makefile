OBJS = Clock.o SoftBody.o Main.o
CC = clang++
CFLAGS = -Wall -c $(DEBUG)
LDFLAGS = -Wall $(DEBUG) -framework GLUT -framework OpenGL -framework Cocoa

softbody : $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o softbody

Clock.o : Clock.hpp Clock.cpp
	$(CC) $(CFLAGS) Clock.cpp

SoftBody.o : SoftBody.hpp SoftBody.cpp FloatUtil.hpp
	$(CC) $(CFLAGS) SoftBody.cpp

Main.o : Main.cpp
	$(CC) $(CFLAGS) Main.cpp

clean:
	rm *.o softbody
	rm -f *~
