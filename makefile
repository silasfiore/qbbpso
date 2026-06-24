CC = gcc
CFLAGS = -std=c99 -I. -Wall -Werror -I.
LDFLAGS = -lm
TARGET = qbbpso
SRCS = qbbpso.c
OBJS = $(SRCS:.c=.o)

all:$(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)