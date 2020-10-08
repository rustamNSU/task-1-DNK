#include "DNK.h"

#include <iostream>

int main()
{
    RNK rnk(100, Nucleotide::G);
    std::cout << rnk.at(19) << rnk[20] << rnk[18] << rnk.capacity() << std::endl;
    for (int i = 0; i < rnk.size(); ++i)
    {
        rnk[i] = Nucleotide::C;
    }
    rnk[20] = Nucleotide::T;
    std::cout << rnk.at(19) << rnk[20] << rnk[18]<< std::endl;
    return 0;
}
