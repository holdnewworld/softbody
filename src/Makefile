OBJS = Clock.o SoftBody.o Main.o
CC = clang++
CFLAGS = -Wall -c $(DEBUG)
LDFLAGS = -Wall $(DEBUG) -framework GLUT -framework OpenGL -framework Cocoa

softbody : $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o softbody

Clock.o : Clock.h Clock.cpp
	$(CC) $(CFLAGS) Clock.cpp

SoftBody.o : SoftBody.h SoftBody.cpp FloatUtil.h
	$(CC) $(CFLAGS) SoftBody.cpp

Main.o : Main.cpp
	$(CC) $(CFLAGS) Main.cpp

clean:
	\rm *.o *~ softbody
