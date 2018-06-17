CXX =mpicxx

export CFLAGS="-fno-stack-protector"
#DEFINES = -DMPICH_IGNORE_CXX_SEEK -DDIM_TWO -DPARALLEL   -lz 
DEFINES = -DQED_BLOCK  -lz
#INC = -I/home/clapa003/HDF-1.6.5/hdf5-1.6.5/hdf5/include
#LIB = -L/home/clapa003/HDF-1.6.5/hdf5-1.6.5/hdf5/lib -lhdf5
#SOURCES = main.cpp oerror.cpp readfile.cpp parameter.cpp random.cpp  fieldinfo.cpp boundaryEB.cpp laser.cpp network.cpp  domain.cpp  density.cpp  diagnostic.cpp collision.cpp
SOURCE = main.cpp readfile.cpp #parameter.cpp

#OBJECTS = main.o   oerror.o readfile.o parameter.o random.o  fieldinfo.o boundaryEB.o laser.o network.o  domain.o  density.o  diagnostic.o collision.o
OBJECTS = main.o readfile.o #parameter.o

CFLAGS = -O3 $(DEFINES)  

all:	$(OBJECTS) 
	$(CXX) -o ../bin/gztrace $(CFLAGS) $(OBJECTS) $(LIB) 

$(OBJECTS):%.o:%.cpp
	$(CXX) $(CFLAGS) $(INC) -c $<  -o $@
	@ echo "$@ created sucessfully"

clean:
	-rm -f *.o ../bin/gztrace ../t_txt/*.txt ../Data/*.txt  ../Data/photon/*.txt

