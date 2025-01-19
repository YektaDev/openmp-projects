echo "Compiling the release binary for the current ARM machine in Release mode, with aggressive optimizations..."

clang++ -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp ./main.cpp -o ./main

