icc -lgsl -qopenmp $1.c
install_name_tool -change @rpath/libiomp5.dylib /opt/intel/lib/libiomp5.dylib a.out
