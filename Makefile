all: rcfsearch.exe krdecomp.exe

rcfsearch.exe: rcfsearch.o rcftable.o finfield.o
	g++ -Wall -Wno-sign-compare -Wno-comment -std=c++11 -O2 -o rcfsearch.exe rcftable.o rcfsearch.o finfield.o

rcftable.o: rcftable.cpp
	g++ -Wall -Wno-sign-compare -Wno-comment -std=c++11 -O2 -c rcftable.cpp

rcfsearch.o: rcfsearch.cpp
	g++ -Wall -Wno-sign-compare -Wno-comment -std=c++11 -O2 -c rcfsearch.cpp 

finfield.o: finfield.cpp
	g++ -Wall -Wno-sign-compare -Wno-comment -std=c++11 -O2 -c finfield.cpp

krdecomp.exe: krdecomp.o glfqchar.o finfield.o
	g++ -Wall -Wno-sign-compare -Wno-comment -std=c++11 -O2  -o krdecomp.exe krdecomp.o glfqchar.o finfield.o

krdecomp.o: krdecomp.cpp
	g++ -Wall -Wno-sign-compare -Wno-comment -std=c++11 -O2 -I "C:/Users/sheil/eigen/eigen-git-mirror" -c krdecomp.cpp

glfqchar.o: glfqchar.cpp
	g++ -Wall -Wno-sign-compare -Wno-comment -std=c++11 -O2 -I "C:/Users/sheil/eigen/eigen-git-mirror" -c glfqchar.cpp


