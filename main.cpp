#include <omp.h>
#include <iostream>

int main() {
#pragma omp parallel num_threads(4)
    {
        std::cout << "Hello and welcome\n";
    }
    return 0;
}
