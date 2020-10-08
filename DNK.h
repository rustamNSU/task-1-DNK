#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <exception>  // std::runtime_error
#include <cstddef>    // size_t (but can do without it)
using std::size_t;

enum class Nucleotide {
    A = 0,
    G = 1,
    C = 2,
    U = 3,
    T = 3
};

std::ostream &operator<<(std::ostream &os, Nucleotide nucleotide)
{
    switch (nucleotide)
    {
        case Nucleotide::A:
            return os << "A";
        case Nucleotide::G:
            return os << "G";
        case Nucleotide::C:
            return os << "C";
        case Nucleotide::T:
            return os << "T";
    }
}

/**********************************************************************//**
 * @class RNK
 * @param DataType is type of element that is stored in data array
 * @param size_ is number of nucleotides in RNK
 * @param nElement is number of nucleotides in one DataType element
 * @param data_ is storage for DataType elements
 *************************************************************************/
template< class DataType = uint8_t>
class RNK{
private:

    std::vector<DataType> data_;             ///< storage
    size_t nElement_ = sizeof(DataType) * 4; ///< nucleotides number in one element
    size_t size_;                            ///< number of nucleotides

public:

    /* Constructors */
    RNK() = default;
    RNK(size_t size, Nucleotide nucleotide);

    /* Methods */
    DataType GetMask(size_t i) const;
    void InsertNucleotide(size_t i, Nucleotide nucleotide, DataType &element);

    /**********************************************************************//**
     * const access method, so return Nucleotide
     * @param index
     * @return Nucleotide
     *************************************************************************/
    Nucleotide at(size_t index) const;

    /**********************************************************************//**
     * const access method, so return Nucleotide
     * @param index
     * @return Nucleotide
     *************************************************************************/
    Nucleotide operator[](size_t index) const;

    /**********************************************************************//**
     * Proxy-reference for single nucleotide in data
     *************************************************************************/
    class NucleotideReference{
        friend class RNK<DataType>;
    private:
        size_t index_;            ///< Nucleotide index (in element), index < nElement_
        DataType& element_;
        Nucleotide nucleotide_;

    public:
        NucleotideReference() = delete;
        NucleotideReference(size_t index, DataType& element);

        operator Nucleotide();
        NucleotideReference& operator=(Nucleotide nucleotide);

        DataType GetMask(size_t index);
        void InsertNucleotide(Nucleotide nucleotide);
    };

    NucleotideReference operator[](size_t index);


    /* Getters */
    size_t size() const;
    size_t capacity() const;
};


template<class DataType>
RNK<DataType>::RNK(size_t size, Nucleotide nucleotide)
{
    size_ = size;
    size_t data_size = static_cast<size_t>(size / nElement_) + 1;

    DataType fillingElement;
    for (int i = 0; i < nElement_; ++i){
        InsertNucleotide(i, nucleotide, fillingElement);
    }
    data_ = std::vector<DataType>(data_size, fillingElement);
}


template<class DataType>
size_t RNK<DataType>::size() const
{
    return size_;
}


template<class DataType>
size_t RNK<DataType>::capacity() const
{
    return data_.capacity() * nElement_;
}


template<class DataType>
DataType RNK<DataType>::GetMask(size_t i) const{
    size_t index = i % nElement_;
    DataType result = 03;                // 00...0011
    result <<= (2*index);
    return result;
}


template<class DataType>
void RNK<DataType>::InsertNucleotide(size_t i, Nucleotide nucleotide, DataType &element)
{
    auto mask = GetMask(i);
    mask = ~mask;
    auto shiftingNucleotide = static_cast<DataType>(nucleotide);
    size_t index = i % nElement_;
    shiftingNucleotide = shiftingNucleotide << 2 * index;
    element = (element & mask) | shiftingNucleotide;
}


template<class DataType>
Nucleotide RNK<DataType>::at(size_t index) const
{
    if (index > size_){
        std::cerr << "Wrong index = " << index << " - out of range (size = " << size_ << ")" << std::endl;
        std::exit(100);
    }
    size_t data_index = index / nElement_;
    auto mask = GetMask(index);
    DataType nucleotide = (data_[data_index] & mask) >> 2 * (index - data_index * nElement_);
    return static_cast<Nucleotide>(nucleotide);
}


template<class DataType>
Nucleotide RNK<DataType>::operator[](size_t index) const
{
    size_t data_index = index / nElement_;
    auto mask = GetMask(index);
    DataType nucleotide = (data_[data_index] & mask) >> 2 * (index - data_index * nElement_);
    return static_cast<Nucleotide>(nucleotide);
}

template<class DataType>
typename RNK<DataType>::NucleotideReference RNK<DataType>::operator[](size_t index)
{
    size_t data_index = index / nElement_;
    return NucleotideReference(index - data_index*nElement_, data_[data_index]);
}


template<class DataType>
RNK<DataType>::NucleotideReference::NucleotideReference(size_t index, DataType &element):
        index_(index), element_(element)
{
    auto mask = GetMask(index);
    DataType nucleotide = (element_ & mask) >> 2 * index;
    nucleotide_ = static_cast<Nucleotide>(nucleotide);
}


template<class DataType>
RNK<DataType>::NucleotideReference::operator Nucleotide()
{
    return nucleotide_;
}

template<class DataType>
typename RNK<DataType>::NucleotideReference& RNK<DataType>::NucleotideReference::operator=(Nucleotide nucleotide)
{
    InsertNucleotide(nucleotide);
    return *this;
}

template<class DataType>
DataType RNK<DataType>::NucleotideReference::GetMask(size_t index)
{
    DataType result = 03;                // 00...0011
    return result << 2 * index;
}

template<class DataType>
void RNK<DataType>::NucleotideReference::InsertNucleotide(Nucleotide nucleotide)
{
    nucleotide_ = nucleotide;
    auto mask = GetMask(index_);
    mask = ~mask;
    auto shiftingNucleotide = static_cast<DataType>(nucleotide);
    shiftingNucleotide = shiftingNucleotide << 2 * index_;
    element_ = (element_ & mask) | shiftingNucleotide;
}

