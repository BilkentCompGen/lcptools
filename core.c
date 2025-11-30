/**
 * @file core.c
 * @brief Implementation of the `core` struct and its associated functions.
 *
 * This file contains the implementation of the `core` struct, which is used to
 * represent a sequence of encoded bits for string data. The stuct supports
 * operations such as compression, comparison, and writing/reading to files.
 *
 * Key operations include:
 * - Encoding strings into bit arrays using coefficient-based encoding.
 * - Constructing `core` objects from strings or other `core` objects.
 * - Compressing bit representations to optimize memory usage.
 * - Writing and reading `core` objects to and from files.
 * - Comparing `core` objects with overloaded operators.
 * - Efficiently handling block-wise bit manipulations.
 *
 * @note The `STATS` macro is used to conditionally compile sections of the code
 * that track additional metadata such as `start` and `end` indices for
 * performance analysis.
 */

#include "core.h"

/**
 * @brief Computes the 32-bit MurmurHash3 hash for a given key.
 *
 * This function computes a 32-bit hash of the input data 'key' with the
 * specified length 'len' and an optional seed value. It processes the
 * input in blocks and handles any remaining bytes.
 *
 * @param key Pointer to the data to be hashed.
 * @param len The length of the data in bytes.
 * @param seed An initial seed value for the hash computation.
 * @return The resulting 32-bit hash value.
 */
uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed) {
    const uint8_t *data = (const uint8_t *)key;
    const int nblocks = len / 4;

    uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    // Body: Process blocks of 4 bytes at a time
    const uint32_t *blocks = (const uint32_t *)(data + nblocks * 4);

    for (int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;

        h1 ^= k1;
        h1 = (h1 << 15) | (h1 >> (32 - 15));
        h1 = h1 * 5 + 0xe6546b64;
    }

    // Tail: Process remaining bytes
    const uint8_t *tail = (const uint8_t *)(data + nblocks * 4);

    uint32_t k1 = 0;

    switch (len & 3) {
    case 3:
        k1 ^= tail[2] << 16;
        break;
    case 2:
        k1 ^= tail[1] << 8;
        break;
    case 1:
        k1 ^= tail[0];
        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;
        h1 ^= k1;
    }

    // Finalization: Mix the hash to ensure the last few bits are fully mixed
    h1 ^= len;

    /* fmix32 */
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    return h1;
}

void init_core1(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index) {
    cr->start = start_index;
    cr->end = end_index;
    cr->label = 0;
    cr->label |= ((distance-2) << 6);
    cr->label |= (alphabet[(int)(*begin)] << 4);
    cr->label |= (alphabet[(int)(*(begin+distance-2))] << 2);
    cr->label |= (alphabet[(int)(*(begin+distance-1))]);
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * distance;
}

void init_core2(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index) {
    cr->start = start_index;
    cr->end = end_index;
    cr->label = 0;
    cr->label |= ((distance-2) << 6);
    cr->label |= (rc_alphabet[(int)(*(begin))] << 4);
    cr->label |= (rc_alphabet[(int)(*(begin-distance+2))] << 2);
    cr->label |= (rc_alphabet[(int)(*(begin-distance+1))]);
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * distance;
}

void init_core3(struct core *cr, struct core *begin, uint64_t distance) {
    cr->start = begin->start;
    cr->end = (begin+distance-1)->end;
    cr->bit_rep = 0;
    cr->bit_size = 0;

    for (struct core *it=begin; it<begin+distance; it++) {
        cr->bit_size += it->bit_size;
    }

    int index = 0;
    for (struct core *it = begin+distance-1; begin <= it && index + it->bit_size <= 64; it--) {
        cr->bit_rep |= (it->bit_rep << index);
        index += it->bit_size;
    }

    cr->bit_rep = 0x7FFFFFFFFFFFFFFF & cr->bit_rep;
    cr->bit_size = minimum(cr->bit_size, 63);

    ulabel data[4];
    data[0] = (begin)->label;
    data[1] = (begin+distance-2)->label;
    data[2] = (begin+distance-1)->label;
    data[3] = distance-2;
    cr->label = MurmurHash3_32((void*)data, 4 * sizeof(ulabel), 42);
}

void init_core4(struct core *cr, ubit_size bit_size, uint64_t bit_rep, ulabel label, uint64_t start, uint64_t end) {
    cr->bit_size = bit_size;
    cr->bit_rep = bit_rep;
    cr->label = label;
    cr->start = start;
    cr->end = end;
}

void core_compress(const struct core *left_core, struct core *right_core) {
    if (left_core->bit_rep & 0x8000000000000000) { // if compressing 1-level cores
        uint64_t left_core_3 = (left_core->bit_rep) & 3;
        uint64_t left_core_2 = (left_core->bit_rep >> 2) & 3;
        uint64_t left_core_middle_count = (left_core->bit_rep & 0x7FFFFFFFFFFFFFFF) >> 6;
        uint64_t left_core_1 = (left_core->bit_rep >> 4) & 3;
        
        uint64_t right_core_3 = (right_core->bit_rep) & 3;
        uint64_t right_core_2 = (right_core->bit_rep >> 2) & 3;
        uint64_t right_core_middle_count = (right_core->bit_rep & 0x7FFFFFFFFFFFFFFF) >> 6;
        uint64_t right_core_1 = (right_core->bit_rep >> 4) & 3;
        
        if (left_core_3 != right_core_3) { // if right characters mismatches
            if ((left_core_3 & 1) != (right_core_3 & 1)) {
                right_core->bit_rep = (right_core_3 & 1); // 0b00 + r % 2
            } else {
                right_core->bit_rep = 2 + ((right_core_3 >> 1) & 1); // 0b10 + r % 2
            }
            right_core->bit_size = 2;
        } 
        else if (left_core_2 != right_core_2) { // if middle characters mismatches
            if ((left_core_2 & 1) != (right_core_2 & 1)) {
                right_core->bit_rep = 4 + (right_core_2 & 1); // 0b100 + r % 2
            } else {
                right_core->bit_rep = 6 + ((right_core_2 >> 1) & 1); // 0b110 + r % 2
            }
            right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
        } 
        else if (left_core_middle_count != right_core_middle_count) { // middle character counts mismatches
            if (left_core_middle_count < right_core_middle_count) {
                // compare left_core_1 with right_core_2
                if ((left_core_1 & 1) != (right_core_2 & 1)) {
                    right_core->bit_rep = 4 * (left_core_middle_count + 1) + (right_core_2 & 1); // 2 * 2 * (mid + 1) + r % 2
                } else {
                    right_core->bit_rep = 2 * (2 * (left_core_middle_count + 1) + 1) + ((right_core_2 >> 1) & 1); // 2 * (2 * (mid + 1) + 1) + r % 2
                }
                right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
            } else {
                // compare left_core_2 with right_core_1
                if ((left_core_2 & 1) != (right_core_1 & 1)) {
                    right_core->bit_rep = 4 * (right_core_middle_count + 1) + (right_core_1 & 1);
                } else {
                    right_core->bit_rep = 2 * (2 * (right_core_middle_count + 1) + 1) + ((right_core_1 >> 1) & 1);
                }
                right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
            }
        } 
        else if (left_core_1 != right_core_1) { // left characters mismatches
            if ((left_core_1 & 1) != (right_core_1 & 1)) {
                right_core->bit_rep = 4 * (left_core_middle_count + 1) + (right_core_1 & 1);
            } else {
                right_core->bit_rep = 2 * (2 * (left_core_middle_count + 1) + 1) + ((right_core_1 >> 1) & 1);
            }
            right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
        } 
        else { // they are same
            right_core->bit_rep = 2 * right_core->bit_size;
            right_core->bit_size = (64 - __builtin_clzll(right_core->bit_rep));
        }
    } else { // if compressing upper level (>1) cores
        ubit_size first_differing_index = 64;
        if (left_core->bit_rep != right_core->bit_rep) {
            first_differing_index = __builtin_ctzll(left_core->bit_rep ^ right_core->bit_rep); // trailing zero count (0-index)
        } else {
            first_differing_index = right_core->bit_size;
        }
        first_differing_index = minimum(first_differing_index, minimum(left_core->bit_size, right_core->bit_size));
        right_core->bit_rep = 2 * first_differing_index + ((right_core->bit_rep >> first_differing_index) & 1);
        right_core->bit_size = right_core->bit_rep == 0 ? 2 : (64 - __builtin_clzll(right_core->bit_rep));
        right_core->bit_size = right_core->bit_size < 2 ? 2 : right_core->bit_size;
    }

    // now, the right core is dependent on the left; hence, its coverage spans towards the left
    right_core->start = left_core->start;
}

void print_core(const struct core *cr) {
    if (cr->bit_rep & 0x8000000000000000) { // if printing 1-level cores
        uint64_t middle_count = (0x7FFFFFFFFFFFFFFF & cr->bit_rep) >> 6;
        uint64_t middle_val = (cr->bit_rep >> 2) & 3;
        printf("%ld", ((cr->bit_rep >> 5) & 1));
        printf("%ld", ((cr->bit_rep >> 4) & 1));
        for (uint64_t i=0; i<middle_count; i++) {
            printf("%ld", ((middle_val >> 1) & 1));
            printf("%ld", (middle_val & 1));           
        }
        printf("%ld", ((cr->bit_rep >> 1) & 1));
        printf("%ld", (cr->bit_rep & 1));
    } else {
        for (ubit_size index = cr->bit_size - 1; 0 < index; index--) {
            printf("%ld", ((cr->bit_rep >> index) & 1));
        }
        printf("%ld", (cr->bit_rep & 1));
    }
}

// core comparison operator implementation

int core_eq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep == rhs->bit_rep;
}

int core_neq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep != rhs->bit_rep;
}

int core_gt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep > rhs->bit_rep;
}

int core_lt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep < rhs->bit_rep;
}

int core_geq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep >= rhs->bit_rep;
}

int core_leq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep <= rhs->bit_rep;
}