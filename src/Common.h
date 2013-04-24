#ifndef __Common_h__
#define __Common_h__

//==============================================================================
// System Inclusions
//==============================================================================
#include <cassert>

//==============================================================================
// Type Definitions
//==============================================================================
typedef unsigned char       uint8_t;
typedef unsigned short      uint16_t;
typedef unsigned int        uint32_t;
typedef unsigned long long  uint64_t;

//==============================================================================
// Macros
//==============================================================================
#define UNREFERENCED_SYMBOL(x) (x)

//==============================================================================
// Global Methods
//==============================================================================
template<typename T> void Swap(T& a, T& b);
uint16_t ShortSwap(const uint16_t i);
uint32_t LongSwap(const uint32_t i);
uint64_t LongLongSwap(const uint64_t i);
float FloatSwap(const float f);
double DoubleSwap(const double d);

//==============================================================================
// Inline Implementation
//==============================================================================
#include "Common.inl"

#endif
