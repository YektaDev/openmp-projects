echo "Compiling the OpenMP project into the release binary for the current macOS machine, with aggressive optimizations..."

clang++ \
  -std=c++26 \
  -O3 -march=native \
  -Xclang \
  -fopenmp \
  -L/opt/homebrew/opt/libomp/lib \
  -I/opt/homebrew/opt/libomp/include \
  -lomp \
  ./main.cpp -o ./main

echo "Compilation was successful."
