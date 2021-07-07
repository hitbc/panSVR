################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/clib/bam_file.c \
../src/clib/binarys_qsort.c \
../src/clib/kthread.c \
../src/clib/utils.c \
../src/clib/vcf_file.c 

OBJS += \
./src/clib/bam_file.o \
./src/clib/binarys_qsort.o \
./src/clib/kthread.o \
./src/clib/utils.o \
./src/clib/vcf_file.o 

C_DEPS += \
./src/clib/bam_file.d \
./src/clib/binarys_qsort.d \
./src/clib/kthread.d \
./src/clib/utils.d \
./src/clib/vcf_file.d 


# Each subdirectory must supply rules for building sources it contributes
src/clib/%.o: ../src/clib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c++17 -I../src/htslib -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


