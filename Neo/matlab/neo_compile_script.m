%% compile the neo_mex code
mex -c CFLAGS=-fPIC -DNUMBITS=32 -largeArrayDims ../multilevelLib/mlkkm.cpp -I../multilevelLib -I../metisLib -cxx
mex -c CFLAGS=-fPIC -DNUMBITS=32 -largeArrayDims ../multilevelLib/wkkm.cpp -I../multilevelLib -I../metisLib -cxx
mex -c CFLAGS=-fPIC -DNUMBITS=32 -largeArrayDims ../programs/neo_mex.cpp -I../multilevelLib -I../metisLib -cxx
mex -c CFLAGS=-fPIC -DNUMBITS=32 -largeArrayDims ../programs/io.cpp -I../multilevelLib -I../metisLib -cxx
mex -largeArrayDims CFLAGS=-fPIC -DNUMBITS=32 -L/p/lib -L. -L../ neo_mex.o mlkkm.o wkkm.o io.o -lmetis -lm -cxx