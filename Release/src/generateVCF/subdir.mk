################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/generateVCF/SignalAssembly.cpp 

OBJS += \
./src/generateVCF/SignalAssembly.o 

CPP_DEPS += \
./src/generateVCF/SignalAssembly.d 


# Each subdirectory must supply rules for building sources it contributes
src/generateVCF/%.o: ../src/generateVCF/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"/home/fenghe/eclipse-workspace/panSV/src/htslib" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


