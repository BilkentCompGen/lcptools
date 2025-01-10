#include <string>
#include <stdint.h>
#include <map>
#include <vector>

// to represent 15 chars = uint32_t, 31 chars = uint64_t, 63 chars = unsigned __int128
typedef uint64_t kmer_type;
typedef std::string::iterator Iter;
typedef std::vector<struct minimizer> Vec;

/**
 * @brief Struct to represent a minimizer
 *
 * This struct represents a minimizer, with members kmer (the encoded version of the kmer) and
 * position, which is the position of the kmer in the sequence. The less-than operator is
 * overloaded to facilitate the process of comparing the encoded kmer value of this kmer with
 * that of another kmer.
 */
struct minimizer {
    kmer_type kmer;
    uint64_t position;
    
    minimizer(kmer_type kmer, uint64_t position) : kmer(kmer), position(position) {}

    bool operator<(const struct minimizer& other) const {
        return kmer < other.kmer;
    }
};

/**
 * @brief Encodes a substring of DNA characters into a unique unsigned integer.
 *
 * This function encodes a sequence of characters (typically representing nucleotides
 * in a DNA string) using a predefined mapping into a 2-bit encoded integer for each
 * character. The final encoded value is returned as an unsigned integer.
 *
 * @param map    Pointer to an array where each index represents the mapping of a character.
 * @param begin  Iterator to the beginning of the substring to encode.
 * @param end    Iterator to the end of the substring to encode.
 * @return A unique unsigned integer representing the encoded DNA sequence.
 */
inline kmer_type encode(int* map, Iter begin, Iter end) {
    kmer_type res = 0;
    for (Iter it = begin; it < end; it++) {
        res *= 4;
        res |= map[static_cast<size_t>(*it)];
    }
    return res;
};

/**
 * @brief Encodes a substring of DNA characters into a unique unsigned integer using
 * reverse complement.
 *
 * This function encodes a sequence of characters (typically representing nucleotides
 * in a DNA string) using a predefined mapping into a 2-bit encoded integer for each
 * character in reversed order. The final encoded value is returned as an unsigned integer.
 *
 * @param rc_map Pointer to an array where each index represents the mapping of a character.
 * @param begin  Iterator to the beginning of the substring to encode.
 * @param end    Iterator to the end of the substring to encode.
 * @return A unique unsigned integer representing the encoded DNA sequence.
 */
inline kmer_type rc_encode(int* rc_map, Iter begin, Iter end) {
    kmer_type res = 0;
    for (Iter it = end-1; begin <= it; it--) {
        res *= 4;
        res |= rc_map[static_cast<size_t>(*it)];
    }
    return res;
};

/**
 * @brief Initializes a character map for nucleotide encoding.
 *
 * This function initializes a mapping of DNA characters (A, T, G, C) and their lowercase
 * counterparts (a, t, g, c) into integers (0, 3, 2, 1, respectively). All other characters
 * are mapped to 0 by default.
 *
 * @param map A 128-element integer array where each index corresponds to a character
 *            in the ASCII table.
 */
void init_map(int map[128]) {
    for (int i = 0; i < 128; i++)
        map[i] = 0;

    map['A'] = 0; map['a'] = 0;
    map['T'] = 3; map['t'] = 3;
    map['G'] = 2; map['g'] = 2;
    map['C'] = 1; map['c'] = 1;
};

/**
 * @brief Initializes a character map for complement nucleotide encoding.
 *
 * This function initializes a mapping of DNA characters (A, T, G, C) and their lowercase
 * counterparts (a, t, g, c) into integers (1, 0, 1, 2, respectively). All other characters
 * are mapped to 0 by default.
 *
 * @param map A 128-element integer array where each index corresponds to a character
 *            in the ASCII table.
 */
void init_rc_map(int map[128]) {
    for (int i = 0; i < 128; i++)
        map[i] = 0;

    map['A'] = 3; map['a'] = 3;
    map['T'] = 0; map['t'] = 0;
    map['G'] = 1; map['g'] = 1;
    map['C'] = 2; map['c'] = 2;
};


/**
 * @brief   Finds and inserts the lexicographically smallest k-mer into a minimizer vector.
 *
 * This function searches through a sequence of k-mers within the given range and
 * identifies the lexicographically smallest k-mer. It then encodes this k-mer and appends it
 * as a minimizer to the `minimizers` vector if it hasn't been added already.
 *
 * @param begin Iterator pointing to the beginning of the sequence.
 * @param end Iterator pointing to the end of the sequence.
 * @param index The index of the first k-mer in the sequence.
 * @param kmerSize The size of the k-mer (number of characters).
 * @param minimizers A vector to store the resulting minimizers.
 * @param map A pointer to an array used to encode k-mers.
 */
void process(Iter begin, Iter end, uint64_t index, int kmerSize, Vec &minimizers, int* map) {
    Iter min_kmer = begin; // default minimal kmer is selected to be first one
    uint64_t min_kmer_index = index;
    Iter cur_kmer = begin+1;  // processing is started with second kmer
    uint64_t cur_kmer_index = index+1;
    
    while (cur_kmer < end) {
        int temp_index = 0;

        while (temp_index < kmerSize && 
                ((*(cur_kmer+temp_index))|0x20) == ((*(min_kmer+temp_index))|0x20))
            temp_index++;

        if (temp_index != kmerSize && 
            ((*(cur_kmer+temp_index))|0x20) < ((*(min_kmer+temp_index))|0x20)) {
            min_kmer = cur_kmer;
            min_kmer_index = cur_kmer_index;
        }

        cur_kmer++;
        cur_kmer_index++;
    }
    
    if (minimizers.empty() || minimizers.back().position != min_kmer_index) {
        minimizers.emplace_back(encode(map, min_kmer, min_kmer+kmerSize), min_kmer_index);
    }
};

/**
 * @brief   Finds and inserts the lexicographically smallest canonical k-mer into a minimizer vector.
 *
 * This function searches through a sequence of k-mers within the given range and
 * identifies the lexicographically smallest k-mer. It then encodes this k-mer and appends it
 * as a minimizer to the `minimizers` vector if it hasn't been added already.
 *
 * @param begin Iterator pointing to the beginning of the sequence.
 * @param end Iterator pointing to the end of the sequence.
 * @param index The index of the first k-mer in the sequence.
 * @param kmerSize The size of the k-mer (number of characters).
 * @param minimizers A vector to store the resulting minimizers.
 * @param map A pointer to an array used to encode k-mers.
 */
void process2(Iter begin, Iter end, uint64_t index, int kmerSize, Vec &minimizers, int* map, int* rc_map) {
    kmer_type min_kmer = std::min(encode(map, begin, begin+kmerSize), rc_encode(rc_map, begin, begin+kmerSize));
    uint64_t min_kmer_index = index;

    Iter cur_kmer = begin+1;  // processing is started with second kmer
    uint64_t cur_kmer_index = index+1;
    
    while (cur_kmer<end) {
        kmer_type fwd_kmer = encode(map, cur_kmer, cur_kmer + kmerSize);
        kmer_type rc_kmer = rc_encode(rc_map, cur_kmer, cur_kmer + kmerSize);

        if (fwd_kmer < min_kmer) {
            min_kmer = fwd_kmer;
            min_kmer_index = cur_kmer_index;
        }

        if (rc_kmer < min_kmer) {
            min_kmer = rc_kmer;
            min_kmer_index = cur_kmer_index;
        }

        cur_kmer++;
        cur_kmer_index++;
    }
    
    if (minimizers.empty() || minimizers.back().position != min_kmer_index) {
        minimizers.emplace_back(min_kmer, min_kmer_index);
    }
};

/**
 * @brief   Finds and inserts the lexicographically smallest k-mer into a minimizer map.
 *
 * This function searches through a sequence of k-mers within the given range and
 * identifies the lexicographically smallest k-mer. It then encodes this k-mer and appends it
 * as a minimizer to the map, where key is the minimizer and value is the vector of integers where
 * the same minimizer appeared in the sequence, and to the vector of minimizers
 *
 * @param begin          Iterator pointing to the beginning of the sequence.
 * @param end            Iterator pointing to the end of the sequence.
 * @param current_index  The index of the first k-mer in the sequence.
 * @param kmerSize       The size of the k-mer (number of characters).
 * @param map            A pointer to an array used to encode k-mers.
 * @param minimizerMap   A map that stores kmers and their position in given sequence (in order).
 */
uint64_t process3(Iter begin, Iter end, int64_t previous_index, uint64_t current_index, int kmerSize, int *map, std::map<kmer_type, std::vector<uint64_t>>& minimizerMap) {
    Iter min_kmer = begin;
    Iter cur_kmer = begin+1;
    uint64_t min_kmer_index = current_index;
    uint64_t cur_kmer_index = current_index + 1;
   
    while (cur_kmer < end) {
        int temp_index = 0;

        while (temp_index < kmerSize && 
                ((*(cur_kmer+temp_index))|0x20) == ((*(min_kmer+temp_index))|0x20))
            temp_index++;

        if (temp_index != kmerSize && 
            ((*(cur_kmer+temp_index))|0x20) < ((*(min_kmer+temp_index))|0x20)) {
            min_kmer = cur_kmer;
            min_kmer_index = cur_kmer_index;
        }

        cur_kmer++;
        cur_kmer_index++;
    }

    if (previous_index == -1 || (uint64_t)previous_index != min_kmer_index) {
        kmer_type kmer = encode(map, min_kmer, min_kmer+kmerSize);
        minimizerMap[kmer].push_back(min_kmer_index);
    }

    return (uint64_t)min_kmer_index;
};

uint64_t process4(Iter begin, Iter end, int64_t previous_index, uint64_t current_index, int kmerSize, int *map, int *rc_map, std::map<kmer_type, std::vector<uint64_t>>& minimizerMap) {
    kmer_type min_kmer = std::min(encode(map, begin, begin+kmerSize), rc_encode(rc_map, begin, begin+kmerSize));
    uint64_t min_kmer_index = current_index;

    Iter cur_kmer = begin+1;  // processing is started with second kmer
    uint64_t cur_kmer_index = current_index+1;
    
    while (cur_kmer<end) {
        kmer_type fwd_kmer = encode(map, cur_kmer, cur_kmer + kmerSize);
        kmer_type rc_kmer = rc_encode(rc_map, cur_kmer, cur_kmer + kmerSize);

        if (fwd_kmer < min_kmer) {
            min_kmer = fwd_kmer;
            min_kmer_index = cur_kmer_index;
        }

        if (rc_kmer < min_kmer) {
            min_kmer = rc_kmer;
            min_kmer_index = cur_kmer_index;
        }

        cur_kmer++;
        cur_kmer_index++;
    }
    
    if (previous_index == -1 || (uint64_t)previous_index != min_kmer_index) {
        kmer_type kmer = encode(map, begin+(min_kmer_index-current_index), begin+(min_kmer_index-current_index+kmerSize));
        minimizerMap[kmer].push_back(min_kmer_index);
    }

    return (uint64_t)min_kmer_index;
};