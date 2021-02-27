CXXFLAGS=-std=c++14 -D_GLIBCXX_DEBUG -g -W -Wall -Wundef -MD -MP

all:
	$(MAKE) run_unit_tests
	$(MAKE) run_main

run_unit_tests: optional_test.pass

optional_test: optional_test.o
	$(CXX) -o $@ $^

optional_test.pass: optional_test
	./optional_test
	touch $@

run_main: main
	./main

main: main.o
	$(CXX) -o $@ $^

clean:
	rm -f main *_test *.o *.pass *.d

-include *.d
