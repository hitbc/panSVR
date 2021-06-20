################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/cpp_lib/Assembler/mantaAssembler.cpp 

OBJS += \
./src/cpp_lib/Assembler/mantaAssembler.o 

CPP_DEPS += \
./src/cpp_lib/Assembler/mantaAssembler.d 


# Each subdirectory must supply rules for building sources it contributes
src/cpp_lib/Assembler/%.o: ../src/cpp_lib/Assembler/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"/home/fenghe/eclipse-workspace/panSV/src/htslib" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


