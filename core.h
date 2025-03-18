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
#include "string.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> 

#define minimum(a, b) ((a) < (b) ? (a) : (b))

typedef uint32_t ulabel;

struct core {
    ulabel label;
    uint64_t start;
    uint64_t end;
};

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
 * @param label Unique label assigned to the core structure.
 * @param start Start index of the substring or sequence represented by the core.
 * @param end End index of the substring or sequence represented by the core.
 */
void init_core4(struct core *cr, ulabel label, uint64_t start, uint64_t end);

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
#define core_eq(lhs, rhs) ((lhs)->label & 3) == ((rhs)->label & 3)

/**
 * @brief Operator overload for not-equal-to comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the two objects are not equal, 0 otherwise.
 */
#define core_neq(lhs, rhs) ((lhs)->label & 3) != ((rhs)->label & 3)

/**
 * @brief Operator overload for greater-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater, 0 otherwise.
 */
#define core_gt(lhs, rhs) ((lhs)->label & 3) > ((rhs)->label & 3)

/**
 * @brief Operator overload for smaller-than comparison between two `core`
 * objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller, 0 otherwise.
 */
#define core_lt(lhs, rhs) ((lhs)->label & 3) < ((rhs)->label & 3)

/**
 * @brief Operator overload for greater-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is greater than or equal, 0 otherwise.
 */
#define core_geq(lhs, rhs) ((lhs)->label & 3) >= ((rhs)->label & 3)

/**
 * @brief Operator overload for smaller-than-or-equal-to comparison between
 * two `core` objects.
 *
 * @param lhs The left-hand side `core` object.
 * @param rhs The right-hand side `core` object.
 * @return 1 if the left-hand object is smaller than or equal, 0 otherwise.
 */
#define core_leq(lhs, rhs) ((lhs)->label & 3) <= ((rhs)->label & 3)

#ifdef __cplusplus
}
#endif

#endif