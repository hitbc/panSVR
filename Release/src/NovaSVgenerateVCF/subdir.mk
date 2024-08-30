################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NovaSVgenerateVCF/ReadHandler.cpp \
../src/NovaSVgenerateVCF/SveHandler.cpp \
../src/NovaSVgenerateVCF/sv_main.cpp \
../src/NovaSVgenerateVCF/sve.cpp 

CPP_DEPS += \
./src/NovaSVgenerateVCF/ReadHandler.d \
./src/NovaSVgenerateVCF/SveHandler.d \
./src/NovaSVgenerateVCF/sv_main.d \
./src/NovaSVgenerateVCF/sve.d 

OBJS += \
./src/NovaSVgenerateVCF/ReadHandler.o \
./src/NovaSVgenerateVCF/SveHandler.o \
./src/NovaSVgenerateVCF/sv_main.o \
./src/NovaSVgenerateVCF/sve.o 


# Each subdirectory must supply rules for building sources it contributes
src/NovaSVgenerateVCF/%.o: ../src/NovaSVgenerateVCF/%.cpp src/NovaSVgenerateVCF/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++17 -I../src/htslib -Ipthread -Im -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-NovaSVgenerateVCF

clean-src-2f-NovaSVgenerateVCF:
	-$(RM) ./src/NovaSVgenerateVCF/ReadHandler.d ./src/NovaSVgenerateVCF/ReadHandler.o ./src/NovaSVgenerateVCF/SveHandler.d ./src/NovaSVgenerateVCF/SveHandler.o ./src/NovaSVgenerateVCF/sv_main.d ./src/NovaSVgenerateVCF/sv_main.o ./src/NovaSVgenerateVCF/sve.d ./src/NovaSVgenerateVCF/sve.o

.PHONY: clean-src-2f-NovaSVgenerateVCF

