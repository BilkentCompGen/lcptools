/**
 * @file lps.h
 * @brief Defines the lps struct and its associated methods for handling locally consistent
 * parsing (LCP) of strings.
 *
 * The lps struct is responsible for performing LCP operations on input strings, constructing
 * cores, and supporting functionalities like parsing, compression (deepening), and memory
 * usage calculations. It includes methods for reading and writing the data to files as well
 * as deepening the LCP to higher levels of compression.
 *
 * The lps struct leverages various helper classes like core, encoding, and hash to manage the
 * string data and its compressed forms. Additionally, it supports both standard and
 * reverse-complement parsing for specialized string handling in bioinformatics and other fields.
 *
 * Key functionalities include:
 * - Parsing an input string or file to extract LCP cores.
 * - Performing multi-level compression of LCP cores (deepening).
 * - Saving and loading LCP cores from files.
 * - Calculating memory usage of the constructed LCP structure.
 *
 * Dependencies:
 * - Requires core.h, encoding.h, hash.h, and constant.h for auxiliary data structures and utilities.
 *
 * Example usage:
 * @code
 *   std::string sequence = "AGCTAGCTAG";
 *   lcp::lps parser(sequence);
 *   parser.deepen();
 *   parser.write("output.lps");
 * @endcode
 *
 * @see core.h
 * @see encoding.h
 * @see hash.h
 * @see constant.h
 *
 * @namespace lcp
 * @struct lps
 *
 * @note Destructor handles clean-up of allocated memory for cores.
 *
 * @author Akmuhammet Ashyralyyev
 * @version 1.0
 * @date 2024-09-14
 *
 */

#ifndef LPS_H
#define LPS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "core.h"
#include "encoding.h"
#include <pthread.h>
#include <stdio.h>
#include <math.h>

#ifndef CONSTANT_FACTOR
#define CONSTANT_FACTOR     1.5
#endif

struct lps {
    int level;
    int size;
    struct core *cores;
};

/**
 * @brief Constructs an lps object from a string.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 */
void init_lps(struct lps *lps_ptr, const char *str, int len);

/**
 * @brief Constructs an lps object from a string.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 * @param offset The length of the offset in which each index will be shifted.
 */
void init_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset);

/**
 * @brief Constructs an lps object from a string, with reverse complement
 * transformation.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 */
void init_lps2(struct lps *lps_ptr, const char *str, int len);
/**
 * @brief Initializes an lps object by reading its contents from a binary file.
 *
 * This function reads the level and size of the `lps` object from the provided file, 
 * allocates memory for the cores array if necessary, and then reads each core's data. 
 * If any error occurs during reading, it will print an error message and terminate 
 * the program or return control, depending on error-handling strategy.
 *
 * @param lps_ptr The `lps` object that will be initialized
 * @param in File pointer to the binary file containing the serialized lps data.
 * 
 * @note The caller must ensure that the `FILE *in` is a valid and open binary file 
 *       for reading. The function will allocate memory for `lps_ptr->cores` if 
 *       `lps_ptr->size > 0`. The caller is responsible for freeing this memory later.
 */
void init_lps3(struct lps *lps_ptr, FILE *in);

/**
 * @brief Constructs an lps object from a string, using split and merge paradigm.
 * The give string will be roughly split into the length of `chunk_size`
 * and processed independently. The final cores will be merged into a single array.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 * @param chunk_size The length of the chunks to be processed.
 */
void init_lps4(struct lps *lps_ptr, const char *str, int len, int lcp_level, int chunk_size);

/**
 * @brief Destructor for the lps object. Frees dynamically allocated memory for cores.
 */
void free_lps(struct lps *lps_ptr);

/**
 * @brief Serializes and writes an lps object to a binary file.
 *
 * This function writes the level, size, and all core objects of the `lps` object 
 * to the specified binary file. Each core's data, including its bit representation, 
 * is written sequentially to the file. The resulting file can later be read to 
 * reconstruct the `lps` object.
 *
 * @param lps_ptr The `lps` object that will be initialized
 * @param out File pointer to the binary file where the lps data will be written.
 * 
 * @note The caller must ensure that the `FILE *out` is a valid and open binary file 
 *       for writing. If the file cannot be written to, the function behavior is undefined.
 */
void write_lps(struct lps *lps_ptr, FILE *out);

/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
int parse1(const char *begin, const char *end, struct core *cores, uint64_t offset);

/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores using complement alphabet.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order based on complement), 
 * it builds and stores these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
int parse2(const char *begin, const char *end, struct core *cores, uint64_t offset);

/**
 * @brief Parses a array of cores to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a array of `core` structures defined by iterators `begin` and `end` and 
 * identifies key segments, called "cores," that represent the (LCP) regions. By analyzing
 * `core` structure relationships in the array (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the `core` array to parse.
 * @param end Iterator pointing to the end of the `core` array to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @return Size of the cores identified in the given string.
 */
int parse3(struct core *begin, struct core *end, struct core *cores);

/**
 * @brief Calculates and returns the memory size used by the `lps` structure.
 *
 * @return The memory size (in bytes) used by the `lps` structure.
 */
int64_t lps_memsize(const struct lps *lps_ptr);

/**
 * @brief Deepens the compression level of the LCP structure. This method compresses the
 * existing cores and finds new cores.
 *
 * @param lps_ptr The `lps` object that will be parsed over.
 * @return 1 if successful in deepening the structure, 0 otherwise.
 */
int lps_deepen1(struct lps *lps_ptr);

/**
 * @brief Deepens the compression level of the LCP structure to a specific level.
 *
 * @param lps_ptr The `lps` object that will be parsed over.
 * @param lcp_level The target compression level to deepen to.
 * @return 1 if deepening was successful, 0 otherwise.
 */
int lps_deepen(struct lps *lps_ptr, int lcp_level);

/**
 * @brief Deepens the compression level of the LCP structure in parallel. This method 
 * compresses the existing cores and finds new cores.
 *
 * @param lps_ptr The `lps` object that will be parsed over.
 * @param thread_number The number of theads to run the DCT.
 * @return 1 if successful in deepening the structure, 0 otherwise.
 */
int lps_deepen1_parallel(struct lps *lps_ptr, int thread_number);

/**
 * @brief Deepens the compression level of the LCP structure to a specific level in threads.
 *
 * @param lps_ptr The `lps` object that will be parsed over.
 * @param lcp_level The target compression level to deepen to.
 * @param thread_number The number of theads to run the DCT.
 * @return 1 if deepening was successful, 0 otherwise.
 */
int lps_deepen_parallel(struct lps *lps_ptr, int lcp_level, int thread_number);

/**
 * @brief Outputs the representation of a `lcp` pointer.
 *
 * This function iterates over the cores in the `cores` array and outputs them.
 * It prints each cores as bit representation, after initially printing the LCP level.
 *
 * @param lps_ptr A pointer to the `lps` object to be output.
 */
void print_lps(const struct lps *lps_ptr);

/**
 * @brief Equality operator for comparing two lcp::lps objects.
 *
 * This function compares the sizes of the core arrays of both objects, and then
 * compares each core element by dereferencing the core pointers. If all elements match,
 * the function returns true; otherwise, it returns false.
 *
 * @param lhs The left-hand side lps object to compare.
 * @param rhs The right-hand side lps object to compare.
 * @return 1 if both lps objects are equal, 0 otherwise.
 */
int lps_eq(const struct lps *lhs, const struct lps *rhs);

/**
 * @brief Inequality operator for comparing two lcp::lps objects.
 *
 * This function first checks if the sizes of the core arrays are different.
 * If they are, the objects are not equal. It then checks each core element.
 * If any element is different, the function returns true; otherwise, it returns false.
 *
 * @param lhs The left-hand side lps object to compare.
 * @param rhs The right-hand side lps object to compare.
 * @return 1 if the lps objects are not equal, 0 otherwise.
 */
int lps_neq(const struct lps *lhs, const struct lps *rhs);


#ifdef __cplusplus
}
#endif

#endif
