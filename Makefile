CXX = g++

BIG_INLINE_LIMIT = 20000
OTHERFLAGS = -DNumInsertStates=1 -DVERSION="2.0" 
CFLAGS1 =  -DUSE_CTENSOR -DUSE_MEA_NUSSINOV -DUSE_TMPL_ON_TRANSITION -DUSE_FLOAT_SCORE
VIENNA   = ./vienna
CXXFLAGS = -O3 -DNDEBUG $(OTHERFLAGS) -funroll-loops $(CFLAGS1) -finline-limit=$(BIG_INLINE_LIMIT) -I$(VIENNA) -L$(VIENNA) -I./ -L./ 

TARGETS = picxaa-r

OBJS = Main.o AlifoldMEA.o  McCaskill.o vienna/energy_param.o 

.PHONY : all

all : 
	make $(TARGETS)

picxaa-r : $(OBJS)
	$(CXX) $(CXXFLAGS) -lm -o $@ $(OBJS)

Main.o :  MultiSequence.h ProbabilisticModel.h ScoreType.h Sequence.h FileBuffer.h SparseMatrix.h Defaults.h SafeVector.h AlignGraph.h Sequence.h BPPMatrix.hpp AlifoldMEA.h  
McCaskill.o: McCaskill.hpp $(VIENNA)/energy_param.hpp Util.hpp Beta.hpp ScoreType.h
energy_param.o: $(VIENNA)/energy_param.hpp 
AlifoldMEA.o: nrutil.h Util.hpp Beta.hpp BPPMatrix.hpp MultiSequence.h Sequence.h SafeVector.h


.PHONY : clean
clean:
	rm -f $(TARGETS) *.o *.o~ *.h~ *.hpp~ *.cpp~ *.cc~ $(VIENNA)/*.o $(VIENNA)/*.o~

