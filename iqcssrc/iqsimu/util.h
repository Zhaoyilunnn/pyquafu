#pragma once

#include <iostream>
#include <vector>
#include <time.h>
#include <regex>
#include <type_traits>
#include <bitset>
#include <random>
#include <chrono>
#include <set>
#define TICK(x) auto bench_##x = std::chrono::steady_clock::now();
#define TOCK(x) std::cout << #x ": " << std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - bench_##x).count() << "s" << std::endl;


namespace Utils{

/*----------------stl function------------------*/

std::vector<int> randomArr(size_t length, size_t max){
    srand((unsigned)time(NULL));
    std::vector<int> arr(length);
    for (size_t i=0;i < arr.size();i++){
        arr[i] = rand()%max;
    }
    return arr;
}

std::vector<double> randomDoubleArr(size_t length){
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0., 1.);
    
    std::vector<double> randarr;
    for (auto n = 0; n < length;++n){
        randarr.push_back(distr(eng));
    }
    return randarr;
}

int randomint(int min, int max){
    srand((unsigned)time(NULL));
    return (rand()%(max - min + 1)) + min;
}

static uint32_t randomize(uint32_t i){
    i = (i^61)^(i>>16);
    i *=9;
    i ^= i<<4;
    i *= 0x27d4eb2d;
    i ^= i >> 15;
    return i;
}

template<class T>
void printVector(std::vector<T> const &arr){
    for (auto i : arr){
        std::cout << i << " ";
    }
    std::cout << std::endl;
}


template <class T>
vector<T> vectors_intersection(vector<T> v1,vector<T> v2){
	vector<T> v;
	sort(v1.begin(),v1.end());   
	sort(v2.begin(),v2.end());   
	set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v));
	return v;
}

template <class T>
vector<T> vectors_union(vector<T> v1,vector<T> v2){
	vector<T> v;
	sort(v1.begin(),v1.end());   
	sort(v2.begin(),v2.end());   
	set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v));
	return v;
}

template<class T>
void printSet(std::set<T> const &arr){
    for (auto it = arr.begin();it !=arr.end();++it){
        std::cout << (*it) << " ";
    }
    std::cout << std::endl;
}

template <class T>
std::set<T> intersection(std::set<T> const& s1, std::set<T> const& s2){
    std::set<T> s;
    set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter(s, s.begin()));
    return s;
}

template <class T>
bool vectors_includes(vector<T> v1,vector<T> v2){
	vector<T> v;
	sort(v1.begin(),v1.end());   
	sort(v2.begin(),v2.end());   
	return includes(v1.begin(),v1.end(),v2.begin(),v2.end());;
}


template <class T>
vector<size_t> sort_indices(vector<T> &v)
{
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
    return idx;
}



std::vector<std::string> split_string(const std::string& str, char delim) {
    std::vector<std::string> elems;
    auto lastPos = str.find_first_not_of(delim, 0);
    auto pos = str.find_first_of(delim, lastPos);
    while (pos != std::string::npos || lastPos != std::string::npos) {
        elems.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delim, pos);
        pos = str.find_first_of(delim, lastPos);
    }
    return elems;
}

std::vector<std::string> split_string(const std::string& str, char delim, uint num){
    auto end = str.length();
    std::vector<std::string> elems;
    auto lastPos = str.find_first_not_of(delim, 0);
    auto pos = str.find_first_of(delim, lastPos);
    while ((pos != std::string::npos || lastPos != std::string::npos) && elems.size() < num) {
        elems.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delim, pos);
        pos = str.find_first_of(delim, lastPos);
    }

    if ((pos != std::string::npos || lastPos != std::string::npos)){
        elems.push_back(str.substr(lastPos, end));
    }
    return elems;
}

template<class real_t>
std::vector<real_t> find_numbers(const std::string &str){
    std::smatch matchs;
    std::vector<real_t> res;
    std::regex pattern;
    if (std::is_unsigned<real_t>::value){
        pattern = std::regex("\\d+");
    }else if(std::is_floating_point<real_t>::value){
        pattern = std::regex("-?(([1-9]\\d*\\.\\d*)|(0\\.\\d*[1-9]\\d*))|\\d+");
    }

    auto begin = std::sregex_iterator(str.begin(), str.end(), pattern);
    const std::sregex_iterator end;

    for (std::sregex_iterator i = begin; i != end; ++i) {
        std::string match_str = i->str();
         if (std::is_unsigned<real_t>::value){
            res.push_back(std::stoi(match_str));
        }else if (std::is_floating_point<real_t>::value){
            res.push_back(std::stod(match_str));
        }
    }

    return res;
}


/*----------------bit function------------------*/
const std::complex<double> PHASE_YZ[4] = {1, imag_I, -1, -imag_I};

inline static uint popcount(uint x) {
    x = ((x & 0xaaaaaaaaaaaaaaaaUL) >> 1) + (x & 0x5555555555555555UL);
    x = ((x & 0xccccccccccccccccUL) >> 2) + (x & 0x3333333333333333UL);
    x = ((x & 0xf0f0f0f0f0f0f0f0UL) >> 4) + (x & 0x0f0f0f0f0f0f0f0fUL);
    x = ((x & 0xff00ff00ff00ff00UL) >> 8) + (x & 0x00ff00ff00ff00ffUL);
    x = ((x & 0xffff0000ffff0000UL) >> 16) + (x & 0x0000ffff0000ffffUL);
    x = ((x & 0xffffffff00000000UL) >> 32) + (x & 0x00000000ffffffffUL);
    return (uint)x;
}



/*----------------Eigengate---------------------*/
RowMatrixXcd CNOT(4, 4);
RowMatrixXcd H(2, 2);
RowMatrixXcd X(2, 2);
RowMatrixXcd Y(2, 2);
RowMatrixXcd Z(2, 2);

void init_gate(){
    CNOT << 1, 0, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0,
            0, 1, 0, 0;

    H << 1/sqrt(2), 1/sqrt(2), 
        1/sqrt(2), -1/sqrt(2);

    X << 0, 1,
        1, 0;

    Y << 0, -imag_I,
        imag_I, 0;

    Z << 1, 0,
        0, -1;
}

RowMatrixXcd rx(double theta){
    RowMatrixXcd RX(2, 2);
    RX << cos(0.5 * theta), -imag_I * sin(0.5 * theta),
        -imag_I * sin(0.5 * theta), cos(0.5 * theta);
    return RX; 
}

}//namespace Utils