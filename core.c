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
        h1 = (h1 << 15) | (h1 >> (32 - 15)); // it should be (h1 << 13) | (h1 >> (32 - 13)) but left at it is, let it be a legacy :)
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
