#include "../RNK.h"

#include <iostream>
#include "gtest/gtest.h"
#include <bitset>

TEST(DNK, BaseMethods)
{
    DNK dnk(10);
    dnk[8] = Nucleotide::G;
    dnk[9] = Nucleotide::T;

    std::cout << dnk;

    ASSERT_EQ(dnk.GetFirstChain().StringRepresent(),  "AAAAAAAAGT");
    ASSERT_EQ(dnk.GetSecondChain().StringRepresent(), "TTTTTTTTCA");
}

