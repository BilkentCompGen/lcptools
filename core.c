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
 * - Comparing `core` objects with overloaded operators.
 * - Efficiently handling block-wise bit manipulations.
 *
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
static inline uint32_t murmurhash32(const void *key, int len, uint32_t seed) {
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
    int char1 = alphabet[(unsigned char)(*begin)];
    int char2 = alphabet[(unsigned char)(*(begin + distance - 2))];
    int char3 = alphabet[(unsigned char)(*(begin + distance - 1))];
    cr->label = ((distance - 2) << 6) | (char1 << 4) | (char2 << 2) | char3;
}

void init_core2(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index) {
    cr->start = start_index;
    cr->end = end_index;
    int char1 = rc_alphabet[(unsigned char)(*begin)];
    int char2 = rc_alphabet[(unsigned char)(*(begin + distance - 2))];
    int char3 = rc_alphabet[(unsigned char)(*(begin + distance - 1))];
    cr->label = ((distance - 2) << 6) | (char1 << 4) | (char2 << 2) | char3;
}

void init_core3(struct core *cr, struct core *begin, uint64_t distance) {
    cr->start = begin->start;
    cr->end = (begin+distance-1)->end;
    ulabel data[4] = {(begin)->label, (begin+distance-2)->label, (begin+distance-1)->label, distance-2};
    cr->label = murmurhash32((void*)data, 4 * sizeof(ulabel), 42);
}

void init_core4(struct core *cr, ulabel label, uint64_t start, uint64_t end) {
    cr->label = label;
    cr->start = start;
    cr->end = end;
}

void print_core(const struct core *cr) {
    printf("%u", cr->label);
}