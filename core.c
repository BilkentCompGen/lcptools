#include "core.h"


#if LCP_COMPUTE_LABEL == 1

static inline uint32_t rotl32(uint32_t x, int r) {
    return (x << r) | (x >> (32 - r));
}

static inline lcp_label hash4_label(lcp_label w0, lcp_label w1, lcp_label w2, lcp_label w3) {
    const uint32_t c1 = 0xcc9e2d51, c2 = 0x1b873593;
    uint32_t h1 = 42u, k1;

    k1 = w0; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;
    k1 = w1; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;
    k1 = w2; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;
    k1 = w3; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;

    h1 ^= 16u;                 /* len */
    h1 ^= h1 >> 16; h1 *= 0x85ebca6b; h1 ^= h1 >> 13; h1 *= 0xc2b2ae35; h1 ^= h1 >> 16;
    return h1;
}

#endif

void init_core1(struct core *cr, const char *begin, lcp_pos distance, lcp_pos start_index, lcp_pos end_index) {
    cr->start = start_index;
    cr->end = end_index;
#if LCP_COMPUTE_LABEL == 1
    cr->label = 0;
    cr->label |= ((distance-2) << 6);
    cr->label |= (alphabet[(int)(*begin)] << 4);
    cr->label |= (alphabet[(int)(*(begin+distance-2))] << 2);
    cr->label |= (alphabet[(int)(*(begin+distance-1))]);
    cr->bit_rep = 0x8000000000000000 | cr->label;
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
    cr->start = start_index;
    cr->end = end_index;
#if LCP_COMPUTE_LABEL == 1
    cr->label = 0;
    cr->label |= ((distance-2) << 6);
    cr->label |= (rc_alphabet[(int)(*(begin))] << 4);
    cr->label |= (rc_alphabet[(int)(*(begin-distance+2))] << 2);
    cr->label |= (rc_alphabet[(int)(*(begin-distance+1))]);
    cr->bit_rep = 0x8000000000000000 | cr->label;
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
#if LCP_COMPUTE_LABEL == 1
    cr->label = hash4_label(begin->label, (begin + distance - 2)->label, (begin + distance - 1)->label, (lcp_label)(distance - 2));
#endif
}

void init_core4(struct core *cr, ubit_size bit_size, uint64_t bit_rep, lcp_label label, lcp_pos start, lcp_pos end) {
    cr->bit_size = bit_size;
    cr->bit_rep = bit_rep;
#if LCP_COMPUTE_LABEL == 1
    cr->label = label;
#else
    (void)label;
#endif
    cr->start = start;
    cr->end = end;
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
