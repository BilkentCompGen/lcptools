#include "core.h"

#define minimum(a, b) ((a) < (b) ? (a) : (b))

// Seeds for the label hash
#define LCP_LABEL_SEED_LO 42u
#define LCP_LABEL_SEED_HI 0x9747b28cu


#if LCP_LABEL_BITS != 0

static inline uint32_t rotl32(uint32_t x, int r) {
    return (x << r) | (x >> (32 - r));
}

/**
 * @brief MurmurHash3 x86_32 over four 32-bit words, with a configurable seed.
 *
 * The seed is what lets a 64-bit label be built from two independent runs of
 * this same function.
 */
static inline uint32_t hash4_label32(uint32_t seed, uint32_t w0, uint32_t w1, uint32_t w2, uint32_t w3) {
    const uint32_t c1 = 0xcc9e2d51, c2 = 0x1b873593;
    uint32_t h1 = seed, k1;

    k1 = w0; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;
    k1 = w1; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;
    k1 = w2; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;
    k1 = w3; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;

    h1 ^= 16u;                 /* len */
    h1 ^= h1 >> 16; h1 *= 0x85ebca6b; h1 ^= h1 >> 13; h1 *= 0xc2b2ae35; h1 ^= h1 >> 16;
    return h1;
}

/**
 * @brief Computes a core's label from its three sub-labels and its distance.
 *
 * At LABEL=32 this is one MurmurHash3 x86_32 with the historical seed, so the
 * values are exactly what earlier versions produced.
 *
 * At LABEL=64 the label is two independent 32-bit chains laid side by side:
 *
 *   bits  0..31   fed by the low halves of the sub-labels   (LCP_LABEL_SEED_LO)
 *   bits 32..63   fed by the high halves                    (LCP_LABEL_SEED_HI)
 */
static inline lcp_label hash_label(lcp_label w0, lcp_label w1, lcp_label w2, lcp_label w3) {
#if LCP_LABEL_BITS == 64
    uint32_t lo = hash4_label32(LCP_LABEL_SEED_LO, (uint32_t)w0, (uint32_t)w1, (uint32_t)w2, (uint32_t)w3);
    uint32_t hi = hash4_label32(LCP_LABEL_SEED_HI, (uint32_t)(w0 >> 32), (uint32_t)(w1 >> 32), (uint32_t)(w2 >> 32), (uint32_t)w3);
    return ((lcp_label)hi << 32) | (lcp_label)lo;
#else
    return hash4_label32(LCP_LABEL_SEED_LO, (uint32_t)w0, (uint32_t)w1, (uint32_t)w2, (uint32_t)w3);
#endif
}

/**
 * @brief Builds a level 1 label from the packed symbol encoding.
 */
static inline lcp_label seed_label(lcp_label packed) {
#if LCP_LABEL_BITS == 64
    uint32_t hi = hash4_label32(LCP_LABEL_SEED_HI, (uint32_t)packed, 0u, 0u, 0u);
    return ((lcp_label)hi << 32) | (lcp_label)(uint32_t)packed;
#else
    return packed;
#endif
}

#endif

void init_core1(struct core *cr, const char *begin, lcp_pos distance, lcp_pos start_index, lcp_pos end_index) {
#if LCP_POS_BITS != 0    
    cr->start = start_index;
    cr->end = end_index;
#else
    (void)start_index;
    (void)end_index;
#endif
#if LCP_LABEL_BITS != 0
    lcp_label packed = 0;
    packed |= ((distance-2) << 6);
    packed |= (alphabet[(int)(*begin)] << 4);
    packed |= (alphabet[(int)(*(begin+distance-2))] << 2);
    packed |= (alphabet[(int)(*(begin+distance-1))]);
    // bit_rep carries the packed encoding only; the label may be wider
    cr->bit_rep = 0x8000000000000000 | (uint64_t)packed;
    cr->label = seed_label(packed);
#else
    lcp_label label = 0;
    label |= ((distance-2) << 6);
    label |= (alphabet[(int)(*begin)] << 4);
    label |= (alphabet[(int)(*(begin+distance-2))] << 2);
    label |= (alphabet[(int)(*(begin+distance-1))]);
    cr->bit_rep = 0x8000000000000000 | label;
#endif
    cr->bit_size = 2 * distance;
}

void init_core2(struct core *cr, const char *begin, lcp_pos distance, lcp_pos start_index, lcp_pos end_index) {
#if LCP_POS_BITS != 0
    cr->start = start_index;
    cr->end = end_index;
#else
    (void)start_index;
    (void)end_index;
#endif
#if LCP_LABEL_BITS != 0
    lcp_label packed = 0;
    packed |= ((distance-2) << 6);
    packed |= (rc_alphabet[(int)(*(begin))] << 4);
    packed |= (rc_alphabet[(int)(*(begin-distance+2))] << 2);
    packed |= (rc_alphabet[(int)(*(begin-distance+1))]);
    // bit_rep carries the packed encoding only; the label may be wider
    cr->bit_rep = 0x8000000000000000 | (uint64_t)packed;
    cr->label = seed_label(packed);
#else
    lcp_label label = 0;
    label |= ((distance-2) << 6);
    label |= (rc_alphabet[(int)(*(begin))] << 4);
    label |= (rc_alphabet[(int)(*(begin-distance+2))] << 2);
    label |= (rc_alphabet[(int)(*(begin-distance+1))]);
    cr->bit_rep = 0x8000000000000000 | label;
#endif
    cr->bit_size = 2 * distance;
}

void init_core3(struct core *cr, struct core *begin, lcp_pos distance) {
#if LCP_POS_BITS != 0 
    cr->start = begin->start;
    cr->end = (begin+distance-1)->end;
#endif
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
#if LCP_LABEL_BITS != 0
    cr->label = hash_label(begin->label, (begin + distance - 2)->label, (begin + distance - 1)->label, (lcp_label)(distance - 2));
#endif
}

void init_core4(struct core *cr, ubit_size bit_size, uint64_t bit_rep, lcp_label label, lcp_pos start, lcp_pos end) {
    cr->bit_size = bit_size;
    cr->bit_rep = bit_rep;
#if LCP_LABEL_BITS != 0
    cr->label = label;
#else
    (void)label;
#endif
#if LCP_POS_BITS != 0
    cr->start = start;
    cr->end = end;
#else   
    (void)start;
    (void)end;
#endif
}

void print_core(const struct core *cr) {
    if (cr->bit_rep & 0x8000000000000000) { // if printing 1-level cores
        uint64_t middle_count = (0x7FFFFFFFFFFFFFFF & cr->bit_rep) >> 6;
        uint64_t middle_val = (cr->bit_rep >> 2) & 3;
        printf("%" PRIu64 "\n", ((cr->bit_rep >> 5) & 1));
        printf("%" PRIu64 "\n", ((cr->bit_rep >> 4) & 1));
        for (uint64_t i=0; i<middle_count; i++) {
            printf("%" PRIu64 "\n", ((middle_val >> 1) & 1));
            printf("%" PRIu64 "\n", (middle_val & 1));           
        }
        printf("%" PRIu64 "\n", ((cr->bit_rep >> 1) & 1));
        printf("%" PRIu64 "\n", (cr->bit_rep & 1));
    } else {
        for (ubit_size index = cr->bit_size - 1; 0 < index; index--) {
            printf("%" PRIu64 "\n", ((cr->bit_rep >> index) & 1));
        }
        printf("%" PRIu64 "\n", (cr->bit_rep & 1));
    }
}
