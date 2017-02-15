
CFLAGS = -W -Wall -O3
LDLIBS = -lm

.PHONY : clean

test_window_functions : test_window_functions.o window_functions.o

test_window_functions.o : test_window_functions.c window_functions.h

window_functions.o : window_functions.c window_functions.h

clean :
	$(RM) test_window_functions test_window_functions.o window_functions.o
