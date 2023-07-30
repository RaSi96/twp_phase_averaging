#directory variables setup
SRCS	= $(wildcard *.c)
HEADERS = $(wildcard *.h)

#compiler setup
CC	= gcc
CFLAGS	= -Wall -Wextra -Wformat=2 -Wnull-dereference -Wpedantic \
	  -Wuninitialized -Warray-bounds=2 -Wunused-macros -Wcast-align \
	  -Wpadded

LDLIBS	= -lm

#make commands
all:	$(SRCS) $(HEADERS)
		$(CC) $(CCFLAGS) $^ $(LDLIBS) -o twp_avg_single

clean:
		rm -rf all
