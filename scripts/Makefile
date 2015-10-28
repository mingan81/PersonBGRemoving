# Home made makefile


HOST := $(shell hostname)

PAT1=/
PAT2=/ /
LOCDIRNAME=$(subst $(PAT1),,$(lastword $(subst $(PAT1),$(PAT2),$(CURDIR))))
# $(warning $(LOCDIRNAME))

export OPENCV_INC=/usr/include
#export MYSQL_INC=/usr/include
export BOOST_INC=/usr/include
#export CPPUNIT_INC=/usr/include
#export PUGIXML_INC=/home/qgan/Software/pugixml/src/



export CC= g++
export LD= g++ 

ifeq ($(PROFILE), Debug)
	#export CXXFLAGS = -Wall -D_TRACES_ -fpermissive
	CXXFLAGS = -D_TRACES_ -fpermissive
	export CFLAGS = -g -O3 
	EXTENSION = _deb
else
	#export CXXFLAGS = -Wall -fpermissive
	export CXXFLAGS = -fpermissive
	export CFLAGS = -O3
	EXTENSION = _rel
endif



#export IFLAGS= -I/usr/local/include/ -I$(OPENCV_INC)/opencv -I$(MYSQL_INC)/mysql -I$(USR_INC)/mysql++ -I$(CPPUNIT_INC)/cppunit -I$(BOOST_INC)/boost -I$(PUGIXML_INC)
export IFLAGS= -I/usr/local/include/ -I$(OPENCV_INC)/opencv -I$(BOOST_INC)/boost 
export LDPATHS= -L/usr/local/lib/ -L/usr/lib/
export LDLIBS= -lopencv_core -lopencv_contrib -lopencv_highgui -lopencv_imgproc -lopencv_legacy -lopencv_objdetect -lopencv_video -lboost_system -lboost_filesystem -lboost_thread

LOCAL_SRC=../src/
LOCAL_INC=../include/
LOCAL_OBJ=../obj/


LFILENAME=$(basename $(notdir $(wildcard $(LOCAL_SRC)*.cpp)))
LTARGETNAME=$(addsuffix $(EXTENSION), $(LFILENAME))
LSRC=$(addprefix $(LOCAL_SRC), $(addsuffix .cpp, $(LTARGETNAME)))
LOBJ=$(addprefix $(LOCAL_OBJ), $(addsuffix .o, $(LTARGETNAME)))



TARGET=Scripts$(EXTENSION)

EXE=Scripts$(EXTENSION)\

all: $(EXE)
	@echo fin

$(TARGET): $(LOBJ) $(TARGET).o
	@echo "\nlinking $@..."
	@$(LD) $(CFLAGS) $^ $(LDLIBS) -o $@ 


%.o: $(subst obj,src,$(subst $(EXTENSION).o,.cpp,$@)) 
	@echo "\nCompiling $(subst obj,src,$(subst $(EXTENSION).o,.cpp,$@)) ..."
	@$(CC) $(IFLAGS) $(CFLAGS) $(CXXFLAGS) -c $(subst obj,src,$(subst $(EXTENSION).o,.cpp,$@)) -o $@ 


$(TARGET).o: Scripts.cpp Makefile 
	@echo "\ncompiling $*.cpp ..."
	@echo $(DIRECTORY)
	$(CC) $(IFLAGS) $(CFLAGS) $(CXXFLAGS) -c $< -o $@ 
#	@echo "end of $*.cpp compilation...\n"

test: $(TARGET)
	./$(TARGET)

clean:
	rm -rf *$(EXTENSION).o $(LOBJ) $(EXE)
