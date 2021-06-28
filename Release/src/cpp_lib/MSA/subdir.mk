################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/cpp_lib/MSA/abPOA_handler.cpp 

OBJS += \
./src/cpp_lib/MSA/abPOA_handler.o 

CPP_DEPS += \
./src/cpp_lib/MSA/abPOA_handler.d 


# Each subdirectory must supply rules for building sources it contributes
src/cpp_lib/MSA/%.o: ../src/cpp_lib/MSA/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"../src/htslib" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


