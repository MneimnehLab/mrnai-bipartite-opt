SOURCE_DIR=src
BUILD_DIR=build

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	CC=g++
	C++FLAGS= -std=c++11 -pedantic -O3 -Wall #: temporarily remove this to silence unused-var warnings
endif
ifeq ($(UNAME_S),Darwin)
	CC=clang++
	C++FLAGS= -std=c++11 -stdlib=libc++ -pedantic -O3 -Wall #: temporarily remove this to silence unused-var warnings
endif

ALL_OBJ = $(BUILD_DIR)/utils/parser.o 
ALL_OBJ+= $(BUILD_DIR)/components/Weights.o 
ALL_OBJ+= $(BUILD_DIR)/components/PRBDPCore.o 
ALL_OBJ+= $(BUILD_DIR)/components/Window.o 
ALL_OBJ+= $(BUILD_DIR)/components/Config.o
ALL_OBJ+= $(BUILD_DIR)/components/RNAProperties.o 

PROG_NAME=complete
ALL_PERMS_PROG=allperms
CONVG_PROG=converge

$(BUILD_DIR)/%.o : $(SOURCE_DIR)/%.cc
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/utils
	@mkdir -p $(BUILD_DIR)/components
	$(CC) $(C++FLAGS) -c  $< -o $@ 

all:
	make $(PROG_NAME)
	# make $(ALL_PERMS_PROG)
	make $(CONVG_PROG)
	cd rnaup_weights/ && make
	mkdir -p output

$(PROG_NAME): $(ALL_OBJ) $(BUILD_DIR)/complete.o
	$(CC) $(C++FLAGS) -o  $(BUILD_DIR)/$@ $(ALL_OBJ) $(BUILD_DIR)/complete.o

$(ALL_PERMS_PROG): $(ALL_OBJ) $(BUILD_DIR)/allperms.o
	$(CC) $(C++FLAGS) -o  $(BUILD_DIR)/$@ $(ALL_OBJ) $(BUILD_DIR)/allperms.o

$(CONVG_PROG): $(ALL_OBJ) $(BUILD_DIR)/convg.o
	$(CC) $(C++FLAGS) -o  $(BUILD_DIR)/$@ $(ALL_OBJ) $(BUILD_DIR)/convg.o

clean:
	(rm -f $(BUILD_DIR)/*.o;)
	(rm -f $(BUILD_DIR)/components/*.o;)
	(rm -f $(BUILD_DIR)/utils/*.o;)

