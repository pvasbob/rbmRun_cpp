ALL =  MultiDimArraySetToValue.h rbm_variables.h MultiDimArrayAllocate.o  ReadComplex.o rbm_methods.o  rbm_main.o
Target = main
CXX = g++


$(Target): $(ALL)
	$(CXX) -o $(Target) $(ALL)


clean:
	rm -f *.o 



