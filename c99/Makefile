
CFLAGS = -W -Wall -O3
LDLIBS = -lm

.PHONY : clean default

default : test_window_functions test_fft time_window_functions

test_window_functions : test_window_functions.o window_functions.o

test_window_functions.o : test_window_functions.c window_functions.h

test_fft : test_fft.o window_functions.o

test_fft.o : test_fft.c window_functions.h

time_window_functions : time_window_functions.o window_functions.o

time_window_functions.o : time_window_functions.c window_functions.h

window_functions.o : window_functions.c window_functions.h

clean :
	$(RM) test_window_functions test_window_functions.o
	$(RM) test_fft test_fft.o
	$(RM) time_window_functions time_window_functions.o
	$(RM) window_functions.o
