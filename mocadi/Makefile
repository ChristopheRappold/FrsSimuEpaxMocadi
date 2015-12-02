
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBSDIR  := $(shell root-config --libdir)
ROOTLIBS     := -L$(ROOTLIBSDIR) -lCore -lRIO -lHist -lTree -lMatrix -lPhysics -lMathCore -lm -ldl -rdynamic
#$(shell root-config --libs)
LIBS         += $(ROOTLIBS)

.PHONY : clean

FuncLoadTree : FuncLoadTree.cc
	g++ -fPIC $(ROOTCFLAGS) -I$(FULLRECO)/Ana_Event -c FuncLoadTree.cc $(LIBS) -L$(FULLRECO)/Ana_Event -lMCAnaEvent 
	g++ -shared -O FuncLoadTree.o -Wl,-R$(FULLRECO)/Ana_Event,-R$(ROOTLIBSDIR),-E,-soname,libFuncLoadTree.so -o libFuncLoadTree.so $(LIBS) -L$(FULLRECO)/Ana_Event -lMCAnaEvent 

all : TestFunc.c FuncLoadTree
	g++ -o TestFunc TestFunc.c -L./ -lFuncLoadTree
#$(LIBS) -L/home/christophe/Workspace/FullReco_Mainz/src/Ana_Event -lMCAnaEvent 

clean :
	rm -f libFuncLoadTree.so FuncLoadTree.o TestFunc
