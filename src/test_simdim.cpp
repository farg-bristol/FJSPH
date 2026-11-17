#include "Var.h"
#include "VarDefs.h"
#include <iostream>

struct part
{
    StateVecD xi;
};

int main(int argc, char* argv[])
{
    std::cout << std::boolalpha << std::is_standard_layout_v<SPHPart> << '\n';

    std::cout << "sizeof(SPHPart) = " << sizeof(SPHPart) << '\n';

    std::cout << "xi     " << offsetof(SPHPart, xi) << '\n';
    std::cout << "v      " << offsetof(SPHPart, v) << '\n';
    std::cout << "acc    " << offsetof(SPHPart, acc) << '\n';
    std::cout << "rho    " << offsetof(SPHPart, rho) << '\n';
    std::cout << "p      " << offsetof(SPHPart, p) << '\n';
}