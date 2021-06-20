################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/kswlib/kalloc.c \
../src/kswlib/ksw2_dispatch.c \
../src/kswlib/ksw2_extd2_sse.c \
../src/kswlib/ksw2_exts2_sse.c \
../src/kswlib/ksw2_extz2_sse.c \
../src/kswlib/ksw2_ll_sse.c 

OBJS += \
./src/kswlib/kalloc.o \
./src/kswlib/ksw2_dispatch.o \
./src/kswlib/ksw2_extd2_sse.o \
./src/kswlib/ksw2_exts2_sse.o \
./src/kswlib/ksw2_extz2_sse.o \
./src/kswlib/ksw2_ll_sse.o 

C_DEPS += \
./src/kswlib/kalloc.d \
./src/kswlib/ksw2_dispatch.d \
./src/kswlib/ksw2_extd2_sse.d \
./src/kswlib/ksw2_exts2_sse.d \
./src/kswlib/ksw2_extz2_sse.d \
./src/kswlib/ksw2_ll_sse.d 


# Each subdirectory must supply rules for building sources it contributes
src/kswlib/%.o: ../src/kswlib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I"/home/fenghe/eclipse-workspace/panSV/src/htslib" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


