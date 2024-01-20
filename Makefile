ALL = 	ReadInputToCol.o MultiDimArraySetToValue.h  MultiDimArrayPrint.h DotProduct.o \
	MsgToScreen.o DiagMat.o \
	rbm_variables.h MultiDimArrayAllocate.h  ReadComplex.o rbm_methods.o  rbm_main.o 
		
Target = RBM
CXX = g++
Flag = -llapack -lblas


$(Target): $(ALL)
	$(CXX) -o $(Target) $(ALL) $(Flag)


clean:
	rm -f *.o 



