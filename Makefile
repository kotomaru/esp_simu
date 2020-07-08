CXX = g++ -g

all: src/resl.cpp
	$(CXX) -o simu_resl src/resl.cpp -lm
