#!/bin/bash

# ANSI escape codes for colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

project_name="$1"
if [[ -z "$project_name" ]]; then
  echo -e "${RED}Error: Project name is required as the first argument.${NC}"
  exit 1
fi

script_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

src_dir="${script_path}/../src/$project_name"
if [[ ! -d "$src_dir" ]]; then
  echo -e "${RED}Error: Source directory not found: $src_dir${NC}"
  exit 1
fi

echo -e "${GREEN}Compiling the OpenMP project '$project_name' into the debug binary for the current macOS machine...${NC}"

build_dir="build/$project_name"
mkdir -p "$build_dir"

clang++ \
  -std=c++23 \
  -O0 \
  -Xclang \
  -fopenmp \
  -L/opt/homebrew/opt/libomp/lib \
  -I/opt/homebrew/opt/libomp/include \
  -lomp \
  -Wall -Wextra \
  "$src_dir"/main.cpp -o "$build_dir"/main-debug

if [[ $? -eq 0 ]]; then
  echo -e "${GREEN}Compilation was successful.${NC}"
  chmod +x "$build_dir"/main-debug 2> /dev/null
else
  echo -e "${RED}Compilation failed.${NC}"
  exit 1
fi
