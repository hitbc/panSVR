################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/cpp_lib/Assembler/assembler.cpp \
../src/cpp_lib/Assembler/mantaAssembler.cpp 

CPP_DEPS += \
./src/cpp_lib/Assembler/assembler.d \
./src/cpp_lib/Assembler/mantaAssembler.d 

OBJS += \
./src/cpp_lib/Assembler/assembler.o \
./src/cpp_lib/Assembler/mantaAssembler.o 


# Each subdirectory must supply rules for building sources it contributes
src/cpp_lib/Assembler/%.o: ../src/cpp_lib/Assembler/%.cpp src/cpp_lib/Assembler/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++17 -I../src/htslib -Ipthread -Im -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-cpp_lib-2f-Assembler

clean-src-2f-cpp_lib-2f-Assembler:
	-$(RM) ./src/cpp_lib/Assembler/assembler.d ./src/cpp_lib/Assembler/assembler.o ./src/cpp_lib/Assembler/mantaAssembler.d ./src/cpp_lib/Assembler/mantaAssembler.o

.PHONY: clean-src-2f-cpp_lib-2f-Assembler

