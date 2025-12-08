#ifndef __TYPES_H__
#define __TYPES_H__

#include<map>
#include<string>

// signed primitive types
typedef char s8;
typedef short s16;
typedef int s32;
typedef long long s64;

// unsigned primitive types
typedef unsigned char u8;
typedef unsigned short u16;
typedef unsigned int u32;
typedef unsigned long long u64;

// for matrix class tolerances
constexpr static double tol = 1e-12;

// for parsing cmdline args
typedef std::map<std::string, std::string> Params_t;

#endif // __TYPES_H__
