#!/bin/bash

# ANSI escape codes for colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}Compiling the OpenMP project into the debug binary for the current macOS machine...${NC}"

mkdir -p build

clang++ \
  -std=c++23 \
  -Xclang \
  -fopenmp \
  -L/opt/homebrew/opt/libomp/lib \
  -I/opt/homebrew/opt/libomp/include \
  -lomp \
  -Wall -Wextra \
  ./main.cpp -o ./build/main-debug

if [[ $? -eq 0 ]]; then
  echo -e "${GREEN}Compilation was successful.${NC} ${YELLOW}Executing...${NC}"

  chmod +x ./build/main-debug 2> /dev/null
  ./build/main-debug
else
  echo -e "${RED}Compilation failed.${NC}"
  exit 1
fi
