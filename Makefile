#Compiling parameters
PROCCESORS = 4
CXX = g++ #mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
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

clean:
	@rm -rf *.x

oclean:
	@rm -rf .objects/*.o
