HW3_CXX_SRC=$(wildcard ./src/*.cpp)
HW3_HEADER=$(wildcard ./src/*.h)
HW3_OBJ=$(notdir $(patsubst %.cpp,%.o,$(HW3_CXX_SRC)))

IMAGE_LIB_SRC=$(wildcard ./lib/imageIO/*.cpp)
IMAGE_LIB_HEADER=$(wildcard ./lib/imageIO/*.h)
IMAGE_LIB_OBJ=$(notdir $(patsubst %.cpp,%.o,$(IMAGE_LIB_SRC)))

HEADER=$(HW3_HEADER) $(IMAGE_LIB_HEADER)
CXX_OBJ=$(HW3_OBJ) $(IMAGE_LIB_OBJ)

CXX=g++
TARGET=raytrace
CXXFLAGS=-DGLM_FORCE_RADIANS
OPT=-O3

PLATFORM=Linux
INCLUDE=-I./lib/imageIO -fopenmp
LIB=-ljpeg -fopenmp
LDFLAGS=

all: $(TARGET)

$(TARGET): $(CXX_OBJ)
	$(CXX) $(LDFLAGS) $^ $(OPT) $(LIB) -o $@

$(HW3_OBJ):%.o: ./src/%.cpp $(HW3_HEADER)
	$(CXX) -c $(CXXFLAGS) $(OPT) $(INCLUDE) $< -o $@

$(IMAGE_LIB_OBJ):%.o: ./lib/imageIO/%.cpp $(IMAGE_LIB_HEADER)
	$(CXX) -c $(CXXFLAGS) $(OPT) $(INCLUDE) $< -o $@

clean:
	rm -rf *.o $(TARGET)
