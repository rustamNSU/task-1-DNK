#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <functional> // std::hash<int>

#include <exception>  // std::runtime_error
#include <cstddef>    // size_t (but can do without it)
#include <cassert>    // assert
using std::size_t;

enum class Nucleotide {
    A = 0,
    G = 1,
    C = 2,
    U = 3,
    T = 3
};

inline std::ostream& operator<<(std::ostream &os, Nucleotide nucleotide)
{
    switch (nucleotide)
    {
        case Nucleotide::A : return os << "A";
        case Nucleotide::G : return os << "G";
        case Nucleotide::C : return os << "C";
        case Nucleotide::T : return os << "T";
    }
    return os;
}

inline char NucleotideLetter(Nucleotide nucleotide)
{
    switch (nucleotide)
    {
        case Nucleotide::A : return 'A';
        case Nucleotide::G : return 'G';
        case Nucleotide::C : return 'C';
        case Nucleotide::T : return 'T';
    }
    return 'A';
}

/* Forward declare for friend operators */
template<class DataType = uint8_t> class RNK;

template<class DataType>
RNK<DataType> operator+(const RNK<DataType> &rnk1, const RNK<DataType> &rnk2);

template<class DataType>
std::ostream& operator<<(std::ostream &os, const RNK<DataType> &rnk);

template<class DataType>
bool operator==(const RNK<DataType> &rnk1, const RNK<DataType> &rnk2);

template<class DataType>
bool operator!=(const RNK<DataType> &rnk1, const RNK<DataType> &rnk2);








/**********************************************************************//**
 * @class RNK
 * @param DataType is type of element that is stored in data array
 * @param size_ is number of nucleotides in RNK
 * @param nElement is number of nucleotides in one DataType element
 * @param data_ is storage for DataType elements
 *************************************************************************/
template< class DataType>
class RNK{
private:

    std::vector<DataType> data_;             ///< storage
    size_t nElement_ = sizeof(DataType) * 4; ///< nucleotides number in one element
    size_t size_;                            ///< number of nucleotides

public:
    /* Constructors */
    RNK() = default;
    RNK(size_t size, Nucleotide nucleotide = Nucleotide::A);

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

    friend RNK<DataType> operator+ <> (const RNK &rnk1, const RNK &rnk2);
    friend std::ostream& operator<< <> (std::ostream &os, const RNK &rnk);
    friend bool          operator== <> (const RNK &rnk1, const RNK &rnk2);
    friend bool          operator!= <> (const RNK &rnk1, const RNK &rnk2);
    RNK<DataType>        operator~() const;

    std::string StringRepresent() const;
    bool isComplementary (const RNK<DataType> &rnk) const;
    RNK<DataType> split(size_t index);

    void resize(size_t new_size, Nucleotide nucleotide = Nucleotide::A);
    size_t cardinality(Nucleotide nucleotide) const;
    std::unordered_map< Nucleotide, int> cardinality() const;
    void trim(size_t last_index);

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


template<class DataType>
void RNK<DataType>::resize(size_t new_size, Nucleotide nucleotide)
{
    if (new_size < size_){
        size_ = new_size;
        return;
    }
    size_t new_data_size = new_size / nElement_ + 1;
    size_t size = size_;

    size_t tail_length = (size_ % nElement_) ? nElement_ - size_ % nElement_ : 0;
    while (tail_length !=0 && size < new_size)
    {
        (*this)[size] = nucleotide;
        ++size;
        --tail_length;
    }
    if (size == new_size)
    {
        size_ = new_size;
        return;
    }

    DataType filling_element;
    for (int i = 0; i < nElement_; ++i){
        InsertNucleotide(i, nucleotide, filling_element);
    }

    if (data_.size() < new_data_size)
    {
        data_.resize(new_data_size);
    }
    for (int i = size / nElement_; i < data_.size(); ++i){
        data_[i] = filling_element;
    }
    size_ = new_size;
}

template<class DataType>
size_t RNK<DataType>::cardinality(Nucleotide nucleotide) const
{
    size_t result = 0;
    for (size_t i = 0; i < size_; ++i)
    {
        if ((*this)[i] == nucleotide) ++result;
    }
    return result;
}

template<class DataType>
std::unordered_map<Nucleotide, int> RNK<DataType>::cardinality() const
{
    std::unordered_map<Nucleotide, int> result =
            {
                    {Nucleotide::A, 0},
                    {Nucleotide::G, 0},
                    {Nucleotide::C, 0},
                    {Nucleotide::T, 0}
            };
    for (size_t i = 0; i < size_; ++i)
    {
        result[(*this)[i]]++;
    }
    return result;
}

template<class DataType>
void RNK<DataType>::trim(size_t last_index)
{
    if (last_index < size_){ size_ = last_index; }
}

template<class DataType>
RNK<DataType> operator+(const RNK<DataType> &rnk1, const RNK<DataType> &rnk2)
{
    auto result = rnk1;
    result.resize(rnk1.size() + rnk2.size());
    for (int i = rnk1.size(), j = 0; i < result.size(); ++i, ++j)
    {
        result[i] = rnk2[j];
    }
    return result;
}

template<class DataType>
std::ostream& operator<<(std::ostream &os, const RNK<DataType> &rnk)
{
    for (int i = 0; i < rnk.size(); ++i){
        os << rnk[i];
    }
    os << '\n';
    return os;
}

template<class DataType>
std::string RNK<DataType>::StringRepresent() const
{
    std::string result;
    for (int i = 0; i < size_; ++i){
        result.push_back(NucleotideLetter((*this)[i]));
    }
    return result;
}

template<class DataType>
bool operator==(const RNK<DataType> &rnk1, const RNK<DataType> &rnk2)
{
    if (rnk1.size() != rnk2.size()) { return false; }
    size_t full_elements_number = rnk1.size() / rnk1.nElement_;
    for (size_t i = 0; i < full_elements_number; ++i)
    {
        if (rnk1.data_[i] != rnk2.data_[i]) { return false; }
    }

    size_t tail_start = full_elements_number * rnk1.nElement_;
    for (size_t i = tail_start; i < rnk1.size(); ++i)
    {
        if (rnk1[i] != rnk2[i]) { return false;}
    }

    return true;
}

template<class DataType>
RNK<DataType> RNK<DataType>::operator~() const
{
    auto result = *this;
    for (size_t i = 0; i < (size_ / nElement_) + 1; ++i){
        result.data_[i] = ~result.data_[i];
    }
    return result;
}

template<class DataType>
bool operator!=(const RNK<DataType> &rnk1, const RNK<DataType> &rnk2)
{
    return !(rnk1 == rnk2);
}

template<class DataType>
bool RNK<DataType>::isComplementary(const RNK<DataType> &rnk) const
{
    return *this == ~rnk;
}

template<class DataType>
RNK<DataType> RNK<DataType>::split(size_t index)
{
    if (index >= size_) { return RNK<DataType>(); }
    size_t tail_size = size_ - index;
    RNK<DataType> tail(tail_size);
    for (size_t i = tail_size, j = 0; i < size_; ++i, ++j){
        tail[j] = static_cast<Nucleotide>((*this)[i]);
    }
    this->trim(index);
    return tail;
}



/* forward declare for friend function */
template< class DataType> class DNK;

template< class DataType>
std::ostream& operator<<(std::ostream &os, const DNK<DataType> &dnk);

template< class DataType>
DNK<DataType> operator+(const DNK<DataType> &dnk1, const DNK<DataType> &dnk2);

template< class DataType>
DNK<DataType> operator==(const DNK<DataType> &dnk1, const DNK<DataType> &dnk2);

template< class DataType>
DNK<DataType> operator!=(const DNK<DataType> &dnk1, const DNK<DataType> &dnk2);



template< class DataType = uint8_t >
class DNK
{
private:
    RNK<DataType> rnk_;  // first chain in DNK. Second chain is (~rnk)

public:
    DNK();
    DNK(size_t size, Nucleotide nucleotide = Nucleotide::A);

    /* Getters */
    RNK<DataType> GetFirstChain() const;
    RNK<DataType> GetSecondChain() const;  // ~rnk_
    size_t size() const{ return rnk_.size(); }

    friend std::ostream& operator<< <> (std::ostream &os, const DNK &dnk);
    friend DNK<DataType> operator+  <> (const DNK &dnk1, const DNK &dnk2);
    friend DNK<DataType> operator== <> (const DNK &dnk1, const DNK &dnk2);
    friend DNK<DataType> operator!= <> (const DNK &dnk1, const DNK &dnk2);

    typename RNK<DataType>::NucleotideReference operator[](size_t index);
    Nucleotide operator[](size_t index) const;


};


template< class DataType>
DNK<DataType>::DNK() = default;


template< class DataType>
DNK<DataType>::DNK(size_t size, Nucleotide nucleotide) :
        rnk_(RNK<DataType>(size, nucleotide))
{}


template< class DataType>
RNK<DataType> DNK<DataType>::GetFirstChain() const
{
    return rnk_;
}


template< class DataType>
RNK<DataType> DNK<DataType>::GetSecondChain() const
{
    return ~rnk_;
}


template< class DataType>
std::ostream& operator<<(std::ostream &os, const DNK<DataType> &dnk)
{
    os << "First chain: "  << dnk.GetFirstChain()  << "\n"
       << "Second chain: " << dnk.GetSecondChain() << std::endl;
    return os;
}


template< class DataType>
DNK<DataType> operator+(const DNK<DataType> &dnk1, const DNK<DataType> &dnk2)
{
    auto result = dnk1;
    result.rnk_ = dnk1.rnk_ + dnk2.rnk_;
    return result;
}


template< class DataType>
DNK<DataType> operator==(const DNK<DataType> &dnk1, const DNK<DataType> &dnk2)
{
    return dnk1.rnk_ == dnk2.rnk_;
}


template< class DataType>
DNK<DataType> operator!=(const DNK<DataType> &dnk1, const DNK<DataType> &dnk2)
{
    return !(dnk1 == dnk2);
}

template<class DataType>
typename RNK<DataType>::NucleotideReference DNK<DataType>::operator[](size_t index)
{
    return rnk_[index];
}

template<class DataType>
Nucleotide DNK<DataType>::operator[](size_t index) const
{
    return rnk_[index];
}


