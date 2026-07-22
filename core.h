/**
 * @file core.h
 * @brief Header file for the `core` struct, which represents bit-encoded
 * sequences and provides various utilities for encoding, compression, and file
 * I/O operations on bit sequences. (I/O operations will be supported later)
 *
 * This file contains the declaration of the `core` struct along with its member
 * functions, constructors, destructors, and operator overloads. It supports
 * operations for encoding strings into bit sequences, combining multiple
 * sequences, compressing data, and reading/writing from/to files.
 *
 * Key functionalities include:
 * - Encoding sequences of characters (e.g., ACGT) into a compact bit-encoded
 * form.
 * - Supporting reverse complement encoding of DNA sequences.
 * - Saving and loading core from files.
 * - Calculating memory usage of the constructed core structure.
 *
 * Dependencies:
 * - Requires constant.h and encoding.h for auxiliary data structures and
 * utilities.
 *
 * @see constant.h
 * @see encoding.h
 *
 * @namespace lcp
 * @struct core
 *
 *
 * @author Akmuhammet Ashyralyyev
 * @version 1.0
 * @date 2024-09-14
 *
 */

#ifndef CORE_H
#define CORE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "encoding.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> 
#include <inttypes.h>
#include <string.h>

#define DCT_ITERATION_COUNT 1

#define minimum(a, b) ((a) < (b) ? (a) : (b))

typedef uint32_t ubit_size;

#ifdef LABEL64
typedef uint64_t ulabel;
#else
typedef uint32_t ulabel;
#endif

#ifdef LARGE_GENOME
struct core {
    ubit_size bit_size;
    ulabel label;
    uint64_t bit_rep;
    uint64_t start;
    uint64_t end;
};
#else 
struct core {
    ubit_size bit_size;
    ulabel label;
    uint64_t bit_rep;
    uint32_t start;
    uint32_t end;
};
#endif

/**
 * @brief Initializes a core structure with the provided string data and index range.
 * 
 * This function processes a given substring starting at `begin` with a specified 
 * distance (length of the substring) and assigns start and end indices for tracking.
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param begin Pointer to the start of the string data.
 * @param distance Length of the substring to process.
 * @param start_index Start index of the substring within the data.
 * @param end_index End index of the substring within the data.
 */
void init_core1(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index);

/**
 * @brief Initializes a core structure with the provided string data and index range.
 * 
 * Similar to `init_core1`, but initalizes the core structure with reverse complement
 * alphabet encoding.
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param begin Pointer to the start of the string data.
 * @param distance Length of the substring to process.
 * @param start_index Start index of the substring within the data.
 * @param end_index End index of the substring within the data.
 */
void init_core2(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index);

/**
 * @brief Initializes a core structure by combining data from other core structures.
 * 
 * This function initializes a new core structure (`cr`) using a sequence of 
 * `core` objects starting from `begin` with the specified `distance` (number 
 * of `core` objects to process).
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param begin Pointer to the start of the sequence of core structures.
 * @param distance Number of core structures to process in the sequence.
 */
void init_core3(struct core *cr, struct core *begin, uint64_t distance);

/**
 * @brief Directly initializes a core structure with precomputed representation and metadata.
 * 
 * This function allows the initialization of a core structure when the bit 
 * representation, bit size, and label are already computed. Useful for 
 * deserializing or cloning a core structure.
 * 
 * @param cr Pointer to the core structure to initialize.
 * @param bit_size Size of the bit representation in bits.
 * @param bit_rep Pointer to the precomputed bit representation array.
 * @param label Unique label assigned to the core structure.
 * @param start Start index of the substring or sequence represented by the core.
 * @param end End index of the substring or sequence represented by the core.
 */
void init_core4(struct core *cr, ubit_size bit_size, uint64_t bit_rep, ulabel label, uint64_t start, uint64_t end);

/**
 * @brief Output the bit representation of a `core` object.
 *
 * This function iterates over the bits in the `core` object and outputs
 * them to the given output stream. It prints each bit as either 0 or 1,
 * starting from the most significant bit to the least significant bit.
 *
 * @param cr The `core` object to be output.
 */
void print_core(const struct core *cr);

// core comparison operator implementation

/**
 * @brief Operator overload for equality comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the two objects are equal, 0 otherwise.
 */
static inline int core_eq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep == rhs->bit_rep;
}

/**
 * @brief Operator overload for greater-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater, 0 otherwise.
 */
static inline int core_neq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep != rhs->bit_rep;
}

/**
 * @brief Operator overload for smaller-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller, 0 otherwise.
 */
static inline int core_gt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep > rhs->bit_rep;
}

/**
 * @brief Operator overload for not-equal-to comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the two objects are not equal, 0 otherwise.
 */
static inline int core_lt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep < rhs->bit_rep;
}

/**
 * @brief Operator overload for greater-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater than or equal, 0 otherwise.
 */
static inline int core_geq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep >= rhs->bit_rep;
}

/**
 * @brief Operator overload for smaller-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller than or equal, 0 otherwise.
 */
static inline int core_leq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep <= rhs->bit_rep;
}

/**
 * @brief Return the minimum of two ubit_size values.
 *
 * @param a First value.
 * @param b Second value.
 * @return The smaller of a and b.
 */
static inline ubit_size umin(ubit_size a, ubit_size b) { 
    return a < b ? a : b; 
}

/**
 * @brief Compute the bit-length of a 64-bit unsigned integer.
 *
 * Bit-length is defined as the position of the most significant set bit
 * (1-based). Returns 0 if x == 0.
 *
 * Uses __builtin_clzll for efficient leading-zero counting.
 *
 * @param x Input 64-bit unsigned integer.
 * @return Number of significant bits required to represent x.
 */
static inline ubit_size bitlen_u64(uint64_t x) {
    return x ? (ubit_size)(64u - (ubit_size)__builtin_clzll(x)) : 0u;
}

/**
 * @brief Compute the bit-length of x with a minimum result of 2.
 *
 * Ensures the returned bit-length is at least 2, even if x has fewer
 * significant bits (or is zero).
 *
 * @param x Input 64-bit unsigned integer.
 * @return max(bitlen_u64(x), 2).
 */
static inline ubit_size bitlen_min2(uint64_t x) {
    ubit_size bl = bitlen_u64(x);
    return bl < 2u ? 2u : bl;
}

/**
 * @brief Extract a 2-bit symbol from the low 6 bits of a packed value.
 *
 * The low 6 bits are interpreted as three 2-bit slots:
 *   k = 0 → rightmost 2 bits
 *   k = 1 → middle 2 bits
 *   k = 2 → leftmost 2 bits
 *
 * @param rep Packed representation.
 * @param k   Slot index (0–2).
 * @return The 2-bit symbol at slot k.
 */
static inline uint64_t sym2(uint64_t rep, unsigned k) {
    return k ? (rep >> k) & 3ull : rep & 3ull;
}

/**
 * @brief Extract the middle repetition count stored above the low 6 bits.
 *
 * Bits 0–5 are reserved for 2-bit symbols. Bits above 6 contain the
 * repetition count. The most significant bit is masked off before shifting.
 *
 * @param rep Packed representation.
 * @return Middle repetition count.
 */
static inline uint64_t mid_count(uint64_t rep) {
    return (rep & 0x7FFFFFFFFFFFFFFFull) >> 6;
}

/**
 * @brief Emit an encoded index-bit value based on symbol comparison.
 *
 * Computes:
 *   result = i + selected_bit or 2 + i + selected_bit, depending on 
 *  from where the bit is selected
 *
 * If the least significant bits of a2 and b2 differ, select (b2 & 1).
 * Otherwise, select ((b2 >> 1) & 1).
 *
 * @param a2 First 2-bit symbol.
 * @param b2 Second 2-bit symbol.
 * @param i  Base index multiplier.
 * @return Encoded value 2 + i plus chosen bit from b2.
 */
static inline uint64_t emit_idx_bit(uint64_t a2, uint64_t b2, uint64_t i) {
    if ((a2 & 1) != (b2 & 1)) {
        return i + (b2 & 1);
    }
    return 2 + i + ((b2 >> 1) & 1);
}

/**
 * @brief Perform level-1 compression of two adjacent cores.
 *
 * Implements the 3-symbol + middle-count encoding scheme.
 * The caller guarantees that both `left` and `right` are already
 * level-1 encoded.
 *
 * The function compares corresponding components of the packed
 * representations in the following priority order:
 *
 *   1. Rightmost 2-bit symbol (L3 vs R3)
 *   2. Leftmost 2-bit symbol (L2 vs R2)
 *   3. Middle repetition count (Lm vs Rm)
 *      - If counts differ, performs a boundary comparison:
 *          • If Lm < Rm: compare L1 with R2
 *          • Otherwise:  compare L2 with R1
 *        The emitted index depends on the smaller count.
 *   4. Final 2-bit symbol (L1 vs R1)
 *   5. If all components match, emit a “same” encoding.
 *
 * In each mismatch case, the emitted value is:
 *
 *     out = 2 * index + selected_bit
 *
 * where the selected bit is derived from the compared 2-bit symbols.
 *
 * Side effects:
 *   - Overwrites `right->bit_rep` with the compressed result.
 *   - Updates `right->bit_size` to bitlen_min2(out).
 *   - Propagates `left->start` into `right->start`.
 *
 * @param left   Pointer to left core (read-only, level-1 encoded).
 * @param right  Pointer to right core (level-1 encoded, updated in place).
 */
static inline void core_compress_level1(const struct core *left, struct core *right) {
    uint64_t L = left->bit_rep;
    uint64_t R = right->bit_rep;

    uint64_t L3 = sym2(L, 0), L2 = sym2(L, 2), L1 = sym2(L, 4);
    uint64_t R3 = sym2(R, 0), R2 = sym2(R, 2), R1 = sym2(R, 4);

    uint64_t Lm = mid_count(L);
    uint64_t Rm = mid_count(R);

    uint64_t out;

    if (L3 != R3) {
        // base = 0 (2*i with i=0)
        out = emit_idx_bit(L3, R3, 0);
        right->bit_rep = out;
        right->bit_size = 2;        
    }
    else if (L2 != R2) {
        // Your original used base=4/6 -> i=2 (since 2*i = 4)
        out = emit_idx_bit(L2, R2, 4);
        right->bit_rep = out;
        right->bit_size = bitlen_min2(out);
    }
    else if (Lm != Rm) {
        // Preserve your “compare across boundary depending on which count is smaller”
        if (Lm < Rm) {
            // compare left L1 with right R2; index = Lm + 1
            out = emit_idx_bit(L1, R2, 4 * Lm + 4);
        } else {
            // compare left L2 with right R1; index = Rm + 1
            out = emit_idx_bit(L2, R1, 4 * Rm + 4);
        }
        right->bit_rep = out;
        right->bit_size = bitlen_min2(out);
    }
    else if (L1 != R1) {
        // left mismatch: index = Lm + 1
        out = emit_idx_bit(L1, R1, 4 * Lm + 4);
        right->bit_rep = out;
        right->bit_size = bitlen_min2(out);
    }
    else {
        // same
        out = 2ull * (uint64_t)right->bit_size;
        right->bit_rep = out;
        right->bit_size = bitlen_min2(out);
    }

    right->start = left->start;
}


/**
 * @brief Perform upper-level compression via bounded first-difference encoding.
 *
 * The caller guarantees that both `left` and `right` are already
 * upper-level encoded.
 *
 * The algorithm:
 *
 *   1. Compute a comparison bound:
 *        bound = min(left->bit_size, right->bit_size, 64).
 *
 *   2. If the bit representations are identical:
 *        idx = bound.
 *      Otherwise:
 *        - Compute x = left->bit_rep ^ right->bit_rep.
 *        - Find the least significant differing bit:
 *              idx = ctz(x).
 *        - Clamp idx to `bound`.
 *
 *   3. Emit:
 *        out = 2 * idx + bit_at_idx(right)
 *
 *      where bit_at_idx(right) is the bit of right->bit_rep at position idx.
 *
 * Side effects:
 *   - Overwrites `right->bit_rep` with the compressed result.
 *   - Updates `right->bit_size` to bitlen_min2(out).
 *   - Propagates `left->start` into `right->start`.
 *
 * This encoding effectively captures the first differing bit
 * (from the least significant side), bounded by the smaller size.
 *
 * @param left   Pointer to left core (read-only, upper-level encoded).
 * @param right  Pointer to right core (upper-level encoded, updated in place).
 */
static inline void core_compress_upper(const struct core *left, struct core *right) {
    ubit_size bound = umin(right->bit_size, umin(left->bit_size, 64u));

    ubit_size idx;
    if (left->bit_rep == right->bit_rep) {
        idx = bound;
    } else {
        uint64_t x = left->bit_rep ^ right->bit_rep;   // nonzero here
        idx = (ubit_size)__builtin_ctzll(x);
        idx = umin(idx, bound);
    }

    uint64_t out = 2ull * (uint64_t)idx + ((right->bit_rep >> idx) & 1ull);
    right->bit_rep = out;
    right->bit_size = bitlen_min2(out);
    right->start = left->start;
}

#ifdef __cplusplus
}
#endif

#endif
