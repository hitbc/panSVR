################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/analysis.cpp \
../src/deBGA_index.cpp \
../src/getSignalRead.cpp \
../src/jlra_aln.cpp \
../src/main.cpp 

OBJS += \
./src/analysis.o \
./src/deBGA_index.o \
./src/getSignalRead.o \
./src/jlra_aln.o \
./src/main.o 

CPP_DEPS += \
./src/analysis.d \
./src/deBGA_index.d \
./src/getSignalRead.d \
./src/jlra_aln.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"../src/htslib" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


