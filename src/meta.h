/********************************************************
 * Meta type definitions and macros.
 ********************************************************/
#ifndef META_H_
#define META_H_ 1

#include <cstdint>
//#include <cstdbool>
typedef int cstdbool;
#define false 0
#define true 1


/**
 * Returns the number of bits within the given type.
 */
#define bitSize(type) ((uint8_t) sizeof(type) << 3)


/**
 * Returns an array structure definition for the given element type and length.
 */
#define array_struct(type, size) struct { const uint64_t length = (size); const type elements[length]; }


/**
 * Fills the given element count times into the array defined by the given address.
 */
#define array_fill(array, count, element) for (int64_t __index = (count); __index >= 0; --__index) (array)[__index] = (element);


#endif
