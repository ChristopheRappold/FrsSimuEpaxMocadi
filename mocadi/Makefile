
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBSDIR  := $(shell root-config --libdir)
ROOTLIBS     := -L$(ROOTLIBSDIR) -lCore -lRIO -lHist -lTree -lMatrix -lPhysics -lMathCore -lm -ldl -rdynamic
#$(shell root-config --libs)
LIBS         += $(ROOTLIBS)

.PHONY : clean

FuncLoadTree : FuncLoadTree.cc
	g++ -fPIC $(ROOTCFLAGS) -I$(FULLRECO)/Ana_Event -c FuncLoadTree.cc $(LIBS) -L$(FULLRECO)/Ana_Event -lMCAnaEvent 
	g++ -shared -O FuncLoadTree.o -Wl,-R$(FULLRECO)/Ana_Event,-R$(ROOTLIBSDIR),-E,-soname,libFuncLoadTree.so -o libFuncLoadTree.so $(LIBS) -L$(FULLRECO)/Ana_Event -lMCAnaEvent 

FuncLoadTreeMC : FuncLoadTreeMC.cc
	g++ -fPIC $(ROOTCFLAGS) -I$(FULLRECO)/Ana_Event -c FuncLoadTreeMC.cc $(LIBS) -L$(FULLRECO)/Ana_Event -lMCAnaEventG4Sol 
	g++ -shared -O FuncLoadTreeMC.o -o libFuncLoadTreeMC.so -Wl,-R$(FULLRECO)/Ana_Event,-R$(ROOTLIBSDIR),-E,-soname,libFuncLoadTreeMC.so $(LIBS) -L$(FULLRECO)/Ana_Event -lMCAnaEventG4Sol 

FuncLoadTreeMC_nnL : FuncLoadTreeMC_nnL.cc
	g++ -fPIC $(ROOTCFLAGS) -I$(FULLRECO)/Ana_Event -c FuncLoadTreeMC_nnL.cc $(LIBS) -L$(FULLRECO)/Ana_Event -lMCAnaEventG4Sol 
	g++ -shared -O FuncLoadTreeMC_nnL.o -o libFuncLoadTreeMC_nnL.so -Wl,-R$(FULLRECO)/Ana_Event,-R$(ROOTLIBSDIR),-E,-soname,libFuncLoadTreeMC_nnL.so $(LIBS) -L$(FULLRECO)/Ana_Event -lMCAnaEventG4Sol 


all : TestFunc.c FuncLoadTree FuncLoadTreeMC
	g++ -o TestFunc TestFunc.c -L./ -lFuncLoadTree
#$(LIBS) -L/home/christophe/Workspace/FullReco_Mainz/src/Ana_Event -lMCAnaEvent 

clean :
	rm -f libFuncLoadTree.so libFuncLoadTreeMC.so FuncLoadTree.o TestFunc
