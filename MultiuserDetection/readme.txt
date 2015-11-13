The file "example.m" contains a simple simulation of the iFDM for the
multiuser blind channel estimation problem. For description of the
variables and structs, see the file "descr.m".


IMPORTANT NOTE:
Please make sure that you compile using mex the file "pgas_C.cpp" (or
"pgas_C_parallel.cpp if you intend to use several cores using
OpenMP) under the folder ./code/pgas/mex. Move the compiled file
to ./code/pgas/compiled. Make sure that you use the -lgsl option to
compile, since the code makes use of the GSL library.
