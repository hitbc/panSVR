################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/PanSVgenerateVCF/SignalAssembly.cpp \
../src/PanSVgenerateVCF/deBGA_index.cpp \
../src/PanSVgenerateVCF/getSignalRead.cpp \
../src/PanSVgenerateVCF/read_realignment.cpp 

CPP_DEPS += \
./src/PanSVgenerateVCF/SignalAssembly.d \
./src/PanSVgenerateVCF/deBGA_index.d \
./src/PanSVgenerateVCF/getSignalRead.d \
./src/PanSVgenerateVCF/read_realignment.d 

OBJS += \
./src/PanSVgenerateVCF/SignalAssembly.o \
./src/PanSVgenerateVCF/deBGA_index.o \
./src/PanSVgenerateVCF/getSignalRead.o \
./src/PanSVgenerateVCF/read_realignment.o 


# Each subdirectory must supply rules for building sources it contributes
src/PanSVgenerateVCF/%.o: ../src/PanSVgenerateVCF/%.cpp src/PanSVgenerateVCF/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++17 -I../src/htslib -Ipthread -Im -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-PanSVgenerateVCF

clean-src-2f-PanSVgenerateVCF:
	-$(RM) ./src/PanSVgenerateVCF/SignalAssembly.d ./src/PanSVgenerateVCF/SignalAssembly.o ./src/PanSVgenerateVCF/deBGA_index.d ./src/PanSVgenerateVCF/deBGA_index.o ./src/PanSVgenerateVCF/getSignalRead.d ./src/PanSVgenerateVCF/getSignalRead.o ./src/PanSVgenerateVCF/read_realignment.d ./src/PanSVgenerateVCF/read_realignment.o

.PHONY: clean-src-2f-PanSVgenerateVCF

