#pragma once
#include <string>
#include <map>
#include <complex>
#include <cctype> 

struct case_insensitive_compare 
{
    bool operator()(const std::string& a, const std::string& b) const {
        auto ita = a.begin();
        auto itb = b.begin();
        
        while (ita != a.end() && itb != b.end()) {
            char ca = std::tolower(static_cast<unsigned char>(*ita));
            char cb = std::tolower(static_cast<unsigned char>(*itb));
            if (ca < cb) return true;
            if (ca > cb) return false;
            ++ita;
            ++itb;
        }
        return (ita == a.end()) && (itb != b.end());
    }
};

template <class T>
using case_insensitive_map = std::map<std::string, T, case_insensitive_compare>;
template<class T> class material_lib
{
    static inline case_insensitive_map<std::complex<T>> mapping
    {
        {"Air",  std::complex<T>(1)           },
        {"TaN",  std::complex<T>(0.926, 0.044)},
        {"Mo",   std::complex<T>(0.919, 0.007)},
        {"Si",   std::complex<T>(0.999, 0.002)},
        {"Ru",   std::complex<T>(0.883, 0.018)},
        {"LTEM", std::complex<T>(0.977, 0.009)},
    };
public:
    static void set(const std::string& name, std::complex<T> val)
    {
        mapping[name] = val;
    }
    static std::complex<T> get(const std::string& name)
    {
        auto it = mapping.find(name);
        if(it == mapping.end()) throw std::runtime_error(name + " not exists in material lib");
        return it->second;
    } 
};