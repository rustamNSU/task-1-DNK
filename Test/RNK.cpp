#include "../RNK.h"

#include <iostream>
#include "gtest/gtest.h"
#include <bitset>

TEST(RNK, Constructors)
{
    RNK rnk(10);
    ASSERT_EQ(rnk.size(), 10);
}

TEST(RNK, OperatorsTest)
{
    RNK rnk(5, Nucleotide::G);
    rnk.resize(10, Nucleotide::C);
    rnk[9] = Nucleotide::T;
    ASSERT_EQ(rnk.StringRepresent(), std::string{"GGGGGCCCCT"});
}

TEST(RNK, OperatorSum)
{
    RNK rnk1(10, Nucleotide::G);
    RNK rnk2(10, Nucleotide::C);

    rnk1.resize(3);
    rnk2.resize(3);

    rnk1 = rnk1 + rnk2;
    rnk1.resize(10, Nucleotide::T);

    ASSERT_EQ(rnk1.StringRepresent(), std::string{"GGGCCCTTTT"});
}

TEST(RNK, Cardinality)
{
    RNK rnk(10);
    rnk[3] = Nucleotide::G;
    rnk[5] = Nucleotide::C;

    rnk = rnk + rnk + rnk;
    auto card = rnk.cardinality();

    ASSERT_EQ(card[Nucleotide::A], 24);
    ASSERT_EQ(rnk.cardinality(Nucleotide::G), 3);
    ASSERT_EQ(card[Nucleotide::C], 3);
    ASSERT_EQ(card[Nucleotide::T], 0);
}

TEST(RNK, HugeSize)
{
    RNK rnk1(1000, Nucleotide::A);
    RNK rnk2(1000, Nucleotide::G);
    RNK rnk3(1000, Nucleotide::C);
    RNK rnk4(1000, Nucleotide::T);

    RNK rnk = rnk1 + rnk2 + rnk3 + rnk4;  // size = 4000
    for (int i = 0; i < 10; ++i){         // should be 4000*1024
        rnk = rnk + rnk;
    }
    auto card = rnk.cardinality();

    ASSERT_EQ(rnk.size(), 4000*1024);

    ASSERT_EQ(card[Nucleotide::A], 1000*1024);
    ASSERT_EQ(card[Nucleotide::G], 1000*1024);
    ASSERT_EQ(card[Nucleotide::C], 1000*1024);
    ASSERT_EQ(card[Nucleotide::T], 1000*1024);
}

TEST(RNK, AdditionalOperators)
{
    RNK rnk1(10, Nucleotide::G);
    RNK rnk2(10, Nucleotide::C);

    ASSERT_TRUE(rnk1 != rnk2);
    ASSERT_TRUE(~rnk1 == rnk2);
    ASSERT_TRUE(rnk1.isComplementary(rnk2));
}

TEST(RNK, Split)
{
    RNK rnk1(10);
    rnk1[8] = Nucleotide::G;
    rnk1[9] = Nucleotide::G;

    auto rnk2 = rnk1.split(5);
    ASSERT_EQ(rnk1.StringRepresent(), "AAAAA");
    ASSERT_EQ(rnk2.StringRepresent(), "AAAGG");
}

TEST(RNK, TemplateParametr)
{
    RNK<uint64_t> rnk(10, Nucleotide::G);
    rnk.resize(40, Nucleotide::C);

    ASSERT_EQ(rnk.StringRepresent(), "GGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    ASSERT_EQ((~rnk).StringRepresent(), "CCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
}

