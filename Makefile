
CFLAGS = -Wall -Wextra -Werror -O3 -g
LDLIBS = -ljpeg -lfftw3 -lconfig -lm
TARGETS=align
SCRIPTS=align_dir make_mp4

all: ${TARGETS} ${SCRIPTS}

align:	align.c

tar:
	rm -f align.tgz; tar cvfz align.tgz Makefile *.[ch] ${SCRIPTS}

tags:
	ctags -R *.[ch]

check:
	cppcheck -q *.[ch]

install:	${TARGETS} ${SCRIPTS}
	cp ${TARGETS} ${SCRIPTS} ~/bin

clean:
	rm -f *.o ${TARGETS} align.tgz tags
