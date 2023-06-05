CXX = g++
TARGET = offline_analyzer
CXXFLAGS = -Wall -m64 -Wuninitialized -O2 -pthread -march=native -std=c++11 -Iinclude #-fsanitize=address
#CXXFLAGS += -DHAVE_ROOT_XML -DHAVE_ROOT_HTTP -DHAVE_THTTP_SERVER
SRCDIR   = src
HEADDIR  = include
OBJDIR   = obj
BUILD	 = build

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(HEADDIR)/*.hpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS =  $(shell root-config --libs)
ROOTLIBS += -lSpectrum -lutil -lrt -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic -lThread -lXMLParser -lXMLIO -lRHTTP -lm -lz -lpthread



# $(BUILD)/$(TARGET): $(OBJECTS)
$(TARGET): $(OBJECTS)
	@$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $(OBJECTS) $(ROOTLIBS)
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@ $(ROOTLIBS)
	@echo "Compiled "$<" successfully!"

clean:
	rm -rf $(OBJDIR)/*.o $(TARGET)

debug: CXXFLAGS += -g

debug: $(TARGET)
