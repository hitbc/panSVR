################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/cpp_lib/graph.cpp 

OBJS += \
./src/cpp_lib/graph.o 

CPP_DEPS += \
./src/cpp_lib/graph.d 


# Each subdirectory must supply rules for building sources it contributes
src/cpp_lib/%.o: ../src/cpp_lib/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"../src/htslib" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


