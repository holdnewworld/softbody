OBJS = Clock.o SoftBody.o Main.o
CC = clang++
CFLAGS = -Wall -c $(DEBUG)
LDFLAGS = -Wall $(DEBUG) -framework GLUT -framework OpenGL -framework Cocoa

softbody : $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o softbody

Clock.o : Clock.h++ Clock.c++
	$(CC) $(CFLAGS) Clock.c++

SoftBody.o : SoftBody.h++ SoftBody.c++ FloatUtil.h++
	$(CC) $(CFLAGS) SoftBody.c++

Main.o : Main.c++
	$(CC) $(CFLAGS) Main.c++

clean:
	\rm *.o *~ softbody
