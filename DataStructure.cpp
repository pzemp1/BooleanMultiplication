#include <deque>
#include <memory>
#include <string>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <chrono>
#include <list>
#include <vector>
#include <unordered_map>
#include <bitset>
#include <experimental/iterator>
#include <climits>

//https://www.usna.edu/Users/cs/roche/courses/s12si335/u04/index.html

//All numbers are going to be treated as positive numbers
//Currently this class is strictly used for mutliplication
//optimisations.

namespace BigNum {
    class BoolArray {
        public:
            mutable bool sign_;
            inline static bool map_[60][4];
            inline static bool output_[2][2][2];
            inline static bool carry_[2][2][2];

            BoolArray() = default;
            ~BoolArray() = default;

            //Creates
            explicit BoolArray(std::string const &str) {
                sign_ = false;
                dimensions_ = 4*str.size();
                auto pow10 = BoolArray('1',dimensions_);
                
                arr_ = std::make_unique<bool[]>(dimensions_);
                auto start = str.rbegin();
                int power = 0;

                while (start != str.rend()) {
                    auto num = BoolArray(*start,dimensions_);
                    //auto num2 = BoolArray(*start,dimensions_);
                    MultiplyBy10(pow10,power);
                    auto copy10 = pow10;

                    if (power > 0) {
                        num *= pow10;
                    }
                    *this += num;
                    //std::cout << *this << std::endl;
                    ++power;
                    ++start;
                }
                //std::cout << *this << std::endl;
                int i = 0;
                if (arr_[i] != 1) {
                    ++i;
                    while (arr_[i] == 0) {
                        ++i;
                    }


                    dimensions_ = dimensions_-i+1;
                    auto new_arr_ = std::make_unique<bool[]>(dimensions_);
                    std::memcpy(new_arr_.get(),arr_.get()+(i-1),dimensions_*sizeof(bool));
                    std::swap(new_arr_, arr_);
                }
            }

            explicit BoolArray(char const &a,unsigned int const &size) {
                sign_ = false;
                dimensions_ = size;
                arr_ = std::make_unique<bool[]>(dimensions_);
                std::memcpy(arr_.get()+(dimensions_-4),&map_[a][0],sizeof(bool)*4);
            }

            explicit BoolArray(unsigned int const &size) : dimensions_{size} {
                sign_ = false;
                arr_ = std::make_unique<bool[]>(dimensions_);
            }

            explicit BoolArray(BoolArray const &ref, int start) : dimensions_{ref.dimensions_-1} {
                //std::cout << ref << std::endl;
                //std::cout << "Start" << std::endl;
                //std::cout << ref.dimensions_-1 << std::endl;
                //std::cout << dimensions_ << std::endl;
                //std::cout << start << std::endl;
                sign_ = false;
                //std::cout << "Make_unique" << std::endl;
                arr_ = std::make_unique<bool[]>(dimensions_);
                //std::cout << "Make_unique Accomplished" << std::endl;
                std::memcpy(arr_.get(),ref.arr_.get()+start,dimensions_*sizeof(bool));
                //std::cout << *this << std::endl;
            }

            explicit BoolArray(BoolArray const &ref, unsigned int const &start, unsigned int const &n) : dimensions_{n} {
                sign_ = false;
                arr_ = std::make_unique<bool[]>(dimensions_);
                std::memcpy(arr_.get(),ref.arr_.get()+start,n*sizeof(bool));
            }

            BoolArray(BoolArray const &ref) : dimensions_{ref.dimensions_}{
                //Only doing memory allocation here! No setting anything to 0
                sign_ = false;
                arr_ = std::make_unique<bool[]>(dimensions_);
                std::memcpy(arr_.get(),ref.arr_.get(),dimensions_*sizeof(bool));
            }

            //Copy Assignment
            BoolArray& operator=(BoolArray const& ev)
            {
                sign_ = false;
                if (this == &ev) {
                    return *this;
                } else {
                    dimensions_ = ev.dimensions_;
                    arr_ = std::make_unique<bool[]>(dimensions_);
                    std::memcpy(arr_.get(),ev.arr_.get(),dimensions_*sizeof(bool));
                }
                return *this;
            }

            BoolArray(BoolArray &&ref) {
                sign_ = false;
                std::swap(sign_,ref.sign_);
                dimensions_ = 0;
                std::swap(dimensions_,ref.dimensions_);
                arr_ = std::make_unique<bool[]>(0);
                std::swap(arr_, ref.arr_);
            }

            BoolArray operator=(BoolArray &&Orig) noexcept
            {
                sign_ = false;
                if (this != &Orig) {
                    dimensions_ = 0;
                    std::swap(dimensions_,Orig.dimensions_);
                arr_ = std::make_unique<bool[]>(0);
                    std::swap(arr_, Orig.arr_);
                }
                return *this;
            }

            BoolArray operator>>=(unsigned int const &x) {
                auto const start = dimensions_ - x;
                auto new_arr_ = std::make_unique<bool[]>(dimensions_);
                std::memcpy(new_arr_.get()+x,arr_.get(),start*sizeof(bool));
                arr_.swap(new_arr_);                
                return *this;
            }

            BoolArray operator<<=(unsigned int const &x) {
                auto const start = dimensions_ - x;
                auto new_arr_ = std::make_unique<bool[]>(dimensions_);
                std::memcpy(new_arr_.get(),arr_.get()+x,start*sizeof(bool));
                arr_.swap(new_arr_);
                return *this;
            }

            bool operator*() {
                for (unsigned int i = 0; i < dimensions_; ++i) {
                    if (arr_[i]) {
                        return true;
                    }
                }
                return false;
            }
            template<typename T>
            BoolArray operator+=(T &&x) const {
                auto const size = x.dimensions_;
                bool remainder = 0;
                for (int i = size-1; i >= 0; --i) {
                    bool const tmp = arr_[i];
                    arr_[i] = output_[remainder][arr_[i]][x.arr_[i]];
                    remainder = carry_[remainder][tmp][x.arr_[i]];
                }
                return *this;
            }

            BoolArray operator-=(BoolArray &x) const{
                if (x == *this) {
                    for (int i = 0; i < dimensions_; ++i) {
                        arr_[i] = 0;
                    }
                    return *this;
                }

                for (int i = 0; i < x.dimensions_; ++i) {
                    x.arr_[i] = !(x.arr_[i]);
                }

                bool remainder = 0;
                int i = x.dimensions_-1;
                int j = dimensions_-1;
                for (; i >= 0; --i) {
                    bool const tmp = arr_[j];
                    arr_[j] = output_[remainder][arr_[j]][x.arr_[i]];
                    remainder = carry_[remainder][tmp][x.arr_[i]];
                    --j;
                }

                while (j>=0) {
                    bool const tmp = arr_[j];
                    arr_[j] = output_[remainder][arr_[j]][1];
                    remainder = carry_[remainder][tmp][1];
                    --j;
                }

                j = dimensions_-1;
                if (remainder) {
                    while (remainder && j >= 0) {
                        bool const tmp = arr_[j];
                        arr_[j] = output_[remainder][arr_[j]][0];
                        remainder = carry_[remainder][tmp][0];
                        --j;
                    }
                } else {
                    sign_ = true;

                    for (int i = 0; i < dimensions_; ++i) {
                        arr_[i] = !(arr_[i]);
                    }
                }
                
                return *this;
            }


            BoolArray operator*=(BoolArray x) {
                //Need to create 2 copies
                auto tmp = BoolArray(*this);
                auto result = BoolArray(dimensions_);

                while (*tmp) {
                    if (tmp.arr_[tmp.dimensions_-1]) {
                        result += x;
                    }
                    x <<= 1;
                    tmp >>= 1;
                }
                std::swap(arr_,result.arr_);
                return result;
            }

            bool operator[](unsigned int const &x) const {
                return *(arr_.get()+x);
            }


            bool& operator[](unsigned int const &x) {
                return *(arr_.get() + x);
            }

            unsigned int getDimensions() const {
                return dimensions_;
            }

            unsigned int size() const {
                return dimensions_;
            }
                // << n
            void SpecialisedCopy(BoolArray const &a, int const &n) const {
                std::memcpy(arr_.get()+n,a.arr_.get(), a.dimensions_ * sizeof(bool));
            }

            void ChangeSize(unsigned int &x) const {
                dimensions_ = x;
            }


            void SpecialisedSubtract(BoolArray const &x, int const &n) const {
                //Follow Normal subtraction steps
                //1. Flip the bits of x
                for (int i = 0; i < x.dimensions_; ++i) {
                    x.arr_[i] = !(x.arr_[i]);
                }
                //Figure out the order of traversal
                bool remainder = 0;
                int i = (dimensions_ - 1) - n;
                int count = x.dimensions_ - 1;
               // std::cout << count << std::endl;
               // std::cout << i << std::endl;
                for (; count >= 0; --i) {
                    bool const tmp = arr_[i];
                    arr_[i] = output_[remainder][arr_[i]][x.arr_[count]];
                    remainder = carry_[remainder][tmp][x.arr_[count]];
                    --count;
                }
                
                i = dimensions_ - 1 - n;
                count = x.dimensions_ - 1;
                if (remainder) {
                    while (remainder) {
                        bool const tmp = arr_[i];
                        arr_[i] = output_[remainder][arr_[i]][x.arr_[count]];
                        remainder = carry_[remainder][tmp][x.arr_[count]];
                        --count;
                        --i;
                    }
                } else {   
                    for (; count >= 0; --i) {
                        arr_[i] = !(arr_[i]);
                        --count;
                    }
                }
                      
                //If there is a remainder then we add

            }

            void SpecialisedAdd(BoolArray const &x, int const &n) const {
                /* n = 3
                   [0,1,2,3,4,5,6,7,8,9,10,11] = this
                   [0,1,2,3,4,5] = x
                   Lets say we are going to shift 3 times
                   add 3 0's to x
                   [0,1,2,3,4,5,0,0,0]
                   We start the addition
                   [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
                            [0, 1, 2, 3, 4, 5, 0,   0, 0]
                    Our starting "add" is a t index = 12 - 1 - 3
                    i = 8, count = 5, we start the addition from here.
                */
                bool remainder = 0;
                int i = (dimensions_ - 1)- n;
                int count = x.dimensions_-1;
                for (; count >= 0; --i) {
                    bool const tmp = arr_[i];
                    arr_[i] = output_[remainder][arr_[i]][x.arr_[count]];
                    remainder = carry_[remainder][tmp][x.arr_[count]];
                    --count;
                }
                while (remainder && i >= 0) {
                    bool const tmp = arr_[i];
                    arr_[i] = output_[remainder][arr_[i]][0];
                    remainder = carry_[remainder][tmp][0];
                    --i;
                }
            }
            /*
            void AddWithShift(KaratsubaArray const &a, int const &n) {

            }
            */


            /*
            KaratsubaArray& operator|=(KaratsubaArray &x) {

            }

            KaratsubaArray& operator/=(KaratsubaArray &x) {

            }

            KaratsubaArray& operator+=(KaratsubaArray &x) {

            }

            KaratsubaArray& operator-=(KaratsubaArray &x) {

            }
            */




            //auto AddWithShift(BinaryArray &ref, int &start) ->void;

            //auto OrWithShift(BinaryArray &ref, int &start) ->void;
            //auto AddWithNoShift(BinaryArray &red)-> void;

            //auto SubtractWithNoShift(BinaryArray &ref)->void;

            inline auto static ConstructStaticMap() -> void {
                std::string const x = "0123456789";
                int count = 0;
                for (int i = 48; i < 58; ++i) {
                    char a = i;
                    auto b1 = std::bitset<4>(x[count]);

                    for (int j = 0; j < 4; ++j) {
                        map_[i][3-j] = int(b1[j]);
                    }
                    ++count;
                }
                bool const output[2][2][2] = {   {  {0,1}, {1,0}  },
                                                 {  {1,0}, {0,1}   }  };

                bool const carry[2][2][2] = {    {  {0,0}, {0,1}  },
                                                 {  {0,1}, {1,1}  } };

                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        for (int k = 0; k < 2; ++k) {
                            output_[i][j][k] = output[i][j][k];
                            carry_[i][j][k] = carry[i][j][k];
                        }
                    }
                }
            }

            friend std::ostream& operator<<(std::ostream& os, BoolArray const& a) {
                //os << "[";
                std::copy(a.arr_.get(), a.arr_.get()+a.dimensions_,
                        std::experimental::make_ostream_joiner(os, ""));
                //os << "]";
                return os;
            }

            friend bool operator==(BoolArray const &x, BoolArray const &y) {
                return std::equal(x.arr_.get(), x.arr_.get()+x.dimensions_
                                  ,y.arr_.get(), y.arr_.get()+y.dimensions_);
            }

            friend bool operator!=(BoolArray const &x, BoolArray const &y) {
                return !(x==y);
            }

            friend BoolArray operator+(BoolArray &x, BoolArray  &y) {
                if (x.size() <= y.size()) {
                    int max = y.size();
                    int mSize = max+1;
                    max--;
                    int min = x.size();
                    min--;
                    auto ret = BoolArray(mSize);
                    mSize--;
                    bool remainder = 0;

                    for (; min >= 0; --min) {
                        ret.arr_[mSize] = output_[remainder][ y.arr_[max] ][x.arr_[min]];
                        remainder = carry_[remainder][ y.arr_[max] ][ x.arr_[min] ];
                        --max;
                        --mSize;
                    }
                    
                    for (; max >= 0; --max) {
                        ret.arr_[mSize] = output_[remainder][ y.arr_[max] ][ 0 ];
                        remainder = carry_[remainder][ y.arr_[max] ][ 0 ];
                        --mSize;
                    }

                    if (remainder) {
                        ret.arr_[mSize] = 1;
                    } else {
                        auto newRet = BoolArray(ret,1); //[0]
                        return newRet;
                    }
                    return ret;

                } else {
                    int max = x.size();
                    int mSize = max+1;
                    max--;
                    int min = y.size();
                    min--;
                    auto ret = BoolArray(mSize);
                    mSize--;
                    bool remainder = 0;

                    for (; min >= 0; --min) {
                        ret[mSize] = output_[remainder][ x[max] ][y[min] ];
                        remainder = carry_[remainder][ x[max] ][ y[min] ];
                        --max;
                        --mSize;
                    }
                    
                    for (; max >= 0; --max) {
                        ret[mSize] = output_[remainder][ x[max] ][ 0 ];
                        remainder = carry_[remainder][ x[max] ][ 0 ];
                        --mSize;
                    }

                    if (remainder) {
                        ret[mSize] = 1;
                    } else {
                        auto newRet = BoolArray(ret,1);
                        return newRet;
                    }
                    return ret;
                }
            }

            void ReSize(unsigned int &x) const {
                if (x > dimensions_) {
                    auto new_arr_ = std::make_unique<bool[]>(x);
                    //std::cout << "A1" << std::endl;
                    std::memcpy(new_arr_.get()+(x-dimensions_),arr_.get(),dimensions_*sizeof(bool));
                    dimensions_ = x;
                    //std::cout << "A2" << std::endl;
                    //arr_.swap(new_arr_);
                    new_arr_.swap(arr_);
                    new_arr_.reset();
                    //std::cout << "A3" << std::endl;
                }
            }
            template <typename T>
            static auto KaratsubaMultiplication(BoolArray const &a, BoolArray const &b, T &&x) -> BoolArray {
                //need to add some sort of padding. 
                if (a.getDimensions() == 1 && b.getDimensions() == 1) {
                    auto ret = BoolArray(2);
                    if (a[0] & b[0]) {
                        ret[1] = 1;
                    }
                    return ret;

                } else {
                    auto size = std::max(a.size(), b.size());
                    auto half = (size+2-1)/2;
                    if (size != a.size()) {
                        a.ReSize(size);
                    }
                    else if (size != b.size()) {
                        b.ReSize(size);
                    }

                    BoolArray A_left;
                    BoolArray B_left;
                    BoolArray A_right;
                    BoolArray B_right;
                    if (a.size() == 1) {
                        A_left = BoolArray(1);
                        A_right = BoolArray(a);
                    } else {
                        A_left = BoolArray(a,0,a.size()/2); 
                        A_right = BoolArray(a,a.size()/2,(a.size() + 2 - 1)/2);
                    }
                    if (b.size() == 1) {
                        B_left = BoolArray(1);
                        B_right = BoolArray(b);
                    } else {
                        B_left = BoolArray(b,0,b.size()/2); 
                        B_right = BoolArray(b,b.size()/2,(b.size() + 2 - 1)/2);
                    }

                    auto d1 = KaratsubaMultiplication(A_left,B_left,x+1);
                    auto d0 = KaratsubaMultiplication(A_right,B_right,x+1);
                    auto d_01 = KaratsubaMultiplication(A_right + A_left,B_right + B_left,x+1);

                    auto C = BoolArray(2*size);
                    C.SpecialisedAdd(d1,2*half); //
                    C.SpecialisedAdd(d0,0);

                    auto tmp = d1 + d0;
                    //if (x == 0) std::cout << tmp << std::endl;
                    d_01 -= tmp; //tmp can be bigger than d_01;
                    if (d_01.sign_) {
                        C.SpecialisedSubtract(d_01,(size+2-1)/2);
                    } else {
                        C.SpecialisedAdd(d_01,(size+2-1)/2);
                    }
                    return C;
                }
            }
            
        private:
            auto MultiplyBy10(BoolArray &pow10, int &power) -> void {
                if (power > 0) {
                    auto prev = BoolArray(pow10);
                    pow10 <<= 3;
                    prev <<= 1;
                    pow10 += prev;
                }
            }
            mutable unsigned int dimensions_;
            mutable std::unique_ptr<bool[]> arr_;
    };
    
    
    
    //00000000000000000100111001011100
    //00000000000000000010111100010001

}
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

std::chrono::time_point<std::chrono::high_resolution_clock> t;

template<typename T>
void timer(T &&i,T &&j) {
    duration<double, std::milli> ms_double = j - i;
    std::cout << ms_double.count() << "ms" << std::endl;
}
/*
    Testing Rationale
*/
//Testing Arithmetic operators and modifiers etc
/*
void TestAM() {
    auto x = BigNum::BoolArray("1532");
    //Check Left and right
    auto left = BigNum::BoolArray(x,0,8);
    auto right = BigNum::BoolArray(x,8,8);
    //Now lets attach them together
    auto empty = BigNum::BoolArray(16);
    empty.SpecialisedCopy(left,0);
    empty.SpecialisedCopy(right,0);
    if (!(empty == x)) {
        std::cout << "Check Line 469 and 470, test failed" << std::endl;
        return;
    }

}
    */



auto main() -> int {
    BigNum::BoolArray::ConstructStaticMap();
    /*
    //Size Constructor works
    //Char constructor works
    //String Constructor works
    auto x = BigNum::BoolArray("9999");
    auto x3 = BigNum::BoolArray("9999");
    std::cout << x << std::endl;
    //Lets try and split up the BoolArray
    //[00000101] [11111100]
    auto left = BigNum::BoolArray(x,0,8);
    auto right = BigNum::BoolArray(x,8,8);
    std::cout << left << right << std::endl;
    //Now lets put them back together
    auto empty = BigNum::BoolArray(16);
    //[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    //[0,1,2,3,4,5,6,7] = Left
    //[8,9,10,11,12,13,14,15] = Right
    //copies from left to right
    empty.SpecialisedCopy(left,0);
    empty.SpecialisedCopy(right,8);
    std::cout << empty << std::endl;
    //0000010111111100
    //Now we want to move it to the right 4 times
    empty >>= 4;
    //Expected 0000000001011111
    std::cout << empty << std::endl;
    empty <<= 3;
    //Expected 0000001011111000
    std::cout << empty << std::endl;
    auto x1 = BigNum::BoolArray("9973");
    //add x and x2
    std::cout << "Before addition : "<< x << std::endl;
    std::cout << "Before addition : "<< x1 << std::endl;
    x += x1;
    std::cout << "After addition  : "<< x << std::endl;
    //Expected = 0010110010001101
    //Works well
    //Now we will try addition with negative numbers
    std::cout << x3 << " - " << x1 << " = ";
    x3 -= x1;
    std::cout << x3 << std::endl; //this works.
    //010111111100
    //010111111100
    // Now we need to test the functionality of
    //specialAdd, specifies the starting and ending position
    //of where the addition occurs.
    // 0000100111000100
    auto simple = BigNum::BoolArray("2500");
    auto simple2500 = BigNum::BoolArray("2500");
    //01011111
    //01011111000
    auto simple95 = BigNum::BoolArray("95");

    auto num = BigNum::BoolArray(16);
    auto stest = BigNum::BoolArray("95");
    auto compare = BigNum::BoolArray("760");
    std::cout << "A" << std::endl;
    std::cout << stest << std::endl;
    std::cout << "..." << std::endl;
    num.SpecialisedCopy(stest,8);
    num <<= 3;
    std::cout << num << std::endl;


    // Answer is = 0000110010111100
    std::cout << simple << " Add but shift 3 times " << simple95 << " = ";
    simple.SpecialisedAdd(simple95,3);
    std::cout << simple << std::endl;
    //Works perfectly
    //Now need to check SpecialisedSubtract

    std::cout << "subtraction" << std::endl;
    simple2500 -= num;
    std::cout << simple2500 << std::endl;
    */
    
    auto m1 = BigNum::BoolArray("94567899876345673");
    auto m2 = BigNum::BoolArray("98945678998763456");
    unsigned long long x = 945;
    unsigned long long y = 9223372036854775807;
    //auto m3 = BigNum::BoolArray("21434567");
    //std::cout << m3 << std::endl;

    std::cout << m1.size() << std::endl;
    std::cout << m2.size() << std::endl;
    //std::cout << "m3 = " << m3 << std::endl;
    //std::cout << m1 + m2 << std::endl; //[10]

    auto t1 = high_resolution_clock::now();
    for(int i = 0; i < 100; ++i) {
        auto a = x*y;
    }
    //auto a = BigNum::BoolArray::KaratsubaMultiplication(m1,m2,0);
    auto t2 = high_resolution_clock::now();
    timer(t1,t2);
    //std::cout << BigNum::KaratsubaMultiplication(m1,m2,0) << std::endl; 

    return 0;
}
