ALL = 	ReadInputToCol.o MultiDimArraySetToValue.h  MultiDimArrayPrint.h DotProduct.o \
	rbm_variables.h MultiDimArrayAllocate.o  ReadComplex.o rbm_methods.o  rbm_main.o 
		
Target = RBM
CXX = g++


$(Target): $(ALL)
	$(CXX) -o $(Target) $(ALL)


clean:
	rm -f *.o 



