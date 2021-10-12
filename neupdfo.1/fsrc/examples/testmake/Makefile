TEST = test_gcc test_clang
test_gcc: CC = gcc -Wall -Wextra
test_clang: CC = clang -Wall -Wextra
.PHONY: all

all: $(TEST)

test_%: hello_%
	./$<

hello_%: hello.c
	$(CC) $< -o $@
