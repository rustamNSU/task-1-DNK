#pragma once

#include <iostream>
#include <string>
#include <sstream>

enum Nucleobase {
    A = 0,
    G = 1,
    C = 2,
    U = 3
};

enum Nucleobase GetNucleobase(size_t* data, size_t size, uint32_t N)
{
    if (N > size * 31)
    {
        std::cout << " return " << std::endl;
        return A;
    }

    size_t mask = 0x0;
    switch (N % 4)
    {
        case 0:
            mask = 0b11000000;
            break;
        case 1:
            mask = 0b00110000;
            break;
        case 2:
            mask = 0b00001100;
            break;
        case 3:
            mask = 0b00000011;
            break;
    }


//    uint8_t *ndata  = (uint8_t*)(data + N/8) ;
//    uint8_t bt = &ndata;

    uint8_t shift = N % 8;
    size_t* elem = data + N / 8;


    size_t mElem = *elem;

    uint8_t ret = ((mElem >> shift * 8) & mask);
    enum Nucl rett = (enum Nucl) ret;


    stringstream ss;
    ss << ret;
    printf("\r\n elem = %x shift = %i mask = %x ret = %i \r\n", (mElem >> shift * 8), shift, mask, ret);
    return rett;




//    uint8_t ret = &ndata & mask;
//    enum Nucl ret = &ndata + &(ndata+1);

}


class RNA {
    class reference {
        size_t nuclNum;

        reference(size_t num) : nuclNum(num)
        {
        }
    };

public:
    RNA();

    RNA(RNA &);

    RNA(RNA &&);

    reference operatorp[](size_t);
};


int main()
{

    size_t nucl = {0x000000C0};
    GetNucl(&nucl, 1, 0);

