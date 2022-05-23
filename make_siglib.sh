swig -Isiglib/include -lua -c++ siglib.i
gcc -fmax-errors=1 -Isiglib/include -O2 -fPIC -march=native -mavx2 -shared -o siglib.so siglib_wrap.cxx libsiglib.a -lstdc++ -lm -lluajit
