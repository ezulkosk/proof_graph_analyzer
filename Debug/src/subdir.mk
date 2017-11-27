################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/merge_resolution_dimacs_check.cpp \
../src/proof_graph_analyzer.cpp 

OBJS += \
./src/merge_resolution_dimacs_check.o \
./src/proof_graph_analyzer.o 

CPP_DEPS += \
./src/merge_resolution_dimacs_check.d \
./src/proof_graph_analyzer.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


