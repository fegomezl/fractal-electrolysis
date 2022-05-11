#Compiling parameters
PROCCESORS = 4
CXX = g++ #mpic++
FLAGS = -std=c++17 -O3 -Wall #-g -fsanitize=address -fsanitize=leak -fsanitize=undefined
RUN = ./#mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

#Arguments for 'plot'
ifeq (plot,$(firstword $(MAKECMDGOALS)))
  PLOT_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(PLOT_ARGS):;@:)
endif

.PHONY: all main plot clean oclean

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
	@python3 settings/plot.py $(PLOT_ARGS)

clean:
	@rm -rf *.x results/*.pdf results/data/*.dat

oclean:
	@rm -rf .objects/*.o
