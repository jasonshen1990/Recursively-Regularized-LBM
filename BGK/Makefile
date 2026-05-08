NVCC = nvcc
# NVCCFLAGS = -O3 -arch=sm_86 -std=c++14 -diag-suppress 20044
# NVCCFLAGS = -O3 -std=c++14 -arch=sm_89 -fmad=false --prec-sqrt=true --prec-div=true -rdc=true
# NVCCFLAGS = -O3 -std=c++14 -arch=sm_89 -rdc=true
# NVCCFLAGS = -O3 -std=c++14 -arch=sm_80 -fmad=false --prec-sqrt=true --prec-div=true -rdc=true
# NVCCFLAGS = -O3 -std=c++14 -arch=sm_80 -rdc=true
# NVCCFLAGS = -O3 -std=c++14 -arch=sm_80 -fmad=false -rdc=true
# NVCCFLAGS = -O3 -std=c++14 -arch=sm_80 --prec-sqrt=true -rdc=true
NVCCFLAGS = -O3 -std=c++14 -arch=sm_89 -rdc=true --prec-div=true

# NVCCFLAGS = -O3 -std=c++14 -arch=sm_80 -rdc=true
# SRCS = main.cu globals.cu initial.cu collision.cu streaming_old.cu macrovar.cu output.cu 
SRCS = main.cu globals.cu initial.cu collision.cu streaming.cu macrovar.cu output.cu 
OBJS = $(SRCS:.cu=.o)
TARGET = chann

all: $(TARGET)

$(TARGET): $(OBJS)
	$(NVCC) $(NVCCFLAGS) $(OBJS) -o $(TARGET)
	rm -f $(OBJS)

%.o: %.cu rb3d.h parameters.h
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

run: $(TARGET)
	./$(TARGET) 
