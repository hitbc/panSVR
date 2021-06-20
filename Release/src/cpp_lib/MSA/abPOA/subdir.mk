################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/cpp_lib/MSA/abPOA/abPOA_utils.c \
../src/cpp_lib/MSA/abPOA/abpoa.c \
../src/cpp_lib/MSA/abPOA/abpoa_align.c \
../src/cpp_lib/MSA/abPOA/abpoa_graph.c \
../src/cpp_lib/MSA/abPOA/abpoa_plot.c \
../src/cpp_lib/MSA/abPOA/abpoa_seq.c \
../src/cpp_lib/MSA/abPOA/simd_abpoa_align.c \
../src/cpp_lib/MSA/abPOA/simd_check.c \
../src/cpp_lib/MSA/abPOA/sub_example.c 

OBJS += \
./src/cpp_lib/MSA/abPOA/abPOA_utils.o \
./src/cpp_lib/MSA/abPOA/abpoa.o \
./src/cpp_lib/MSA/abPOA/abpoa_align.o \
./src/cpp_lib/MSA/abPOA/abpoa_graph.o \
./src/cpp_lib/MSA/abPOA/abpoa_plot.o \
./src/cpp_lib/MSA/abPOA/abpoa_seq.o \
./src/cpp_lib/MSA/abPOA/simd_abpoa_align.o \
./src/cpp_lib/MSA/abPOA/simd_check.o \
./src/cpp_lib/MSA/abPOA/sub_example.o 

C_DEPS += \
./src/cpp_lib/MSA/abPOA/abPOA_utils.d \
./src/cpp_lib/MSA/abPOA/abpoa.d \
./src/cpp_lib/MSA/abPOA/abpoa_align.d \
./src/cpp_lib/MSA/abPOA/abpoa_graph.d \
./src/cpp_lib/MSA/abPOA/abpoa_plot.d \
./src/cpp_lib/MSA/abPOA/abpoa_seq.d \
./src/cpp_lib/MSA/abPOA/simd_abpoa_align.d \
./src/cpp_lib/MSA/abPOA/simd_check.d \
./src/cpp_lib/MSA/abPOA/sub_example.d 


# Each subdirectory must supply rules for building sources it contributes
src/cpp_lib/MSA/abPOA/%.o: ../src/cpp_lib/MSA/abPOA/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I"/home/fenghe/eclipse-workspace/panSV/src/htslib" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


