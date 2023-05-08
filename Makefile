CXX = g++
LIBS = -L/home/roko15/CBZ/CBG/src
INCPATH = -I/home/roko15/CBZ/CBG/src
RM = rm -f
TARGET=a.out
OBJS=main.o
OUTPUTS=
BACKUPS = *~
a.out:${OBJS} ${OBJS_PLOSE}
	${CXX} ${LIBS} ${OBJS} -lCBG
main.o: main.cxx
	${CXX} ${CXXFLAGS} ${INCPATH} -c main.cxx
clean:
	${RM} ${OBJS}
distclean:
	${RM} ${OBJS} ${TARGET} ${OUTPUTS} ${BACKUPS}
