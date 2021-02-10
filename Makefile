CXXFLAGS=-std=c++14 -D_GLIBCXX_DEBUG -g -W -Wall -Wundef

run_main: main
	./main

main:

clean:
	rm main
