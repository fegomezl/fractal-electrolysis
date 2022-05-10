#Compiling parameters
PROCCESORS = 4
CXX = g++ #mpic++
FLAGS = -std=c++17 -O3 -Wall #-g -fsanitize=address -fsanitize=leak -fsanitize=undefined
RUN = ./#mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all main clean oclean

all: main

main: main.x
	@echo -e 'Running program ... \n'
	@$(RUN)$< 
	@echo -e '\nDone!\n'

main.x: $(DEPENDENCIES)
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ -o $@
	@echo -e 'Done!\n'

.objects/%.o: code/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< -o $@
	@echo -e 'Done!\n'

plot:
	@python3 settings/plot.py

clean:
	@rm -rf *.x results/*.txt results/*.pdf results/*.dat

oclean:
	@rm -rf .objects/*.o
