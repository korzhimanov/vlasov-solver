PYTHONVERSION = 2.7
CXX = g++
CFLAGS = -Wall -g -O2 -fopenmp -I./include -I/usr/include/python$(PYTHONVERSION)
LIBS = -lz -lgsl -lgslcblas -lpython$(PYTHONVERSION)
INCLUDE_DIR = include
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
OUT = vlasov

default: all

$(OBJ_DIR)/mymath.o: $(INCLUDE_DIR)/mymath.h $(SRC_DIR)/mymath.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/mymath.cpp -o $(OBJ_DIR)/mymath.o

$(OBJ_DIR)/file_saving.o: $(INCLUDE_DIR)/file_saving.h $(SRC_DIR)/file_saving.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/file_saving.cpp -o $(OBJ_DIR)/file_saving.o

$(OBJ_DIR)/mesh.o: $(INCLUDE_DIR)/mesh.h $(SRC_DIR)/mesh.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/mesh.cpp -o $(OBJ_DIR)/mesh.o

$(OBJ_DIR)/fdtd.o: $(INCLUDE_DIR)/mesh.h $(INCLUDE_DIR)/fdtd.h $(SRC_DIR)/fdtd.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/fdtd.cpp -o $(OBJ_DIR)/fdtd.o

$(OBJ_DIR)/pfc.o: $(INCLUDE_DIR)/pyinput.h $(INCLUDE_DIR)/mesh.h $(INCLUDE_DIR)/fdtd.h $(INCLUDE_DIR)/pfc.h $(SRC_DIR)/pfc.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/pfc.cpp -o $(OBJ_DIR)/pfc.o

$(OBJ_DIR)/particle.o: $(INCLUDE_DIR)/particle.h $(SRC_DIR)/particle.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/particle.cpp -o $(OBJ_DIR)/particle.o

$(OBJ_DIR)/pyinput.o: $(INCLUDE_DIR)/pfunc.h $(INCLUDE_DIR)/pyinput.h $(SRC_DIR)/pyinput.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/pyinput.cpp -o $(OBJ_DIR)/pyinput.o

$(OBJ_DIR)/pfunc.o: $(INCLUDE_DIR)/pfunc.h $(SRC_DIR)/pfunc.cpp
	$(CXX) $(CFLAGS) $(LFLAGS) -c $(SRC_DIR)/pfunc.cpp -o $(OBJ_DIR)/pfunc.o

$(BIN_DIR)/$(OUT): $(OBJ_DIR)/pfunc.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/mymath.o $(OBJ_DIR)/file_saving.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/fdtd.o $(OBJ_DIR)/pfc.o $(OBJ_DIR)/particle.o $(INCLUDE_DIR)/pfc.h $(INCLUDE_DIR)/mymath.h $(INCLUDE_DIR)/file_saving.h $(INCLUDE_DIR)/fdtd.h $(INCLUDE_DIR)/particle.h $(INCLUDE_DIR)/solver.h $(SRC_DIR)/main.cpp
	$(CXX) $(SRC_DIR)/main.cpp -o $(BIN_DIR)/$(OUT) $(OBJ_DIR)/pfunc.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/mymath.o $(OBJ_DIR)/file_saving.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/fdtd.o $(OBJ_DIR)/pfc.o $(OBJ_DIR)/particle.o $(CFLAGS) $(LFLAGS) $(LIBS)

all: $(OBJ_DIR)/pfunc.o $(OBJ_DIR)/pyinput.o $(OBJ_DIR)/mymath.o $(OBJ_DIR)/file_saving.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/fdtd.o $(OBJ_DIR)/pfc.o $(OBJ_DIR)/particle.o $(BIN_DIR)/$(OUT)

clean:
	rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/$(OUT)
