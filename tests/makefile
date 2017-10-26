all:
	c++ -Wall -std=c++14 -c factorial.c++ 
	c++ -Wall -std=c++14 -I/Users/natj/libs/ tests-main.c++ -c
	c++ -Wall -std=c++14 -I/Users/natj/libs/ tests-main.o factorial.o tests-factorial.c++ -o tests 
	./tests 
