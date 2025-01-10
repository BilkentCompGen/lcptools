#include <string>
#include <vector>
#include <map>
#include <stdint.h>

typedef uint32_t kmer_type;
typedef std::string::iterator Iter;
typedef std::vector<struct syncmer> Vec;

/**
 * @brief Struct to represent a syncmer
 *
 * This struct represents a syncmer, with members kmer (the encoded version of the kmer)
 * and position, which is the position of the kmer in the sequence. The less-than operator is
 * overloaded to facilitate the process of comparing the encoded kmer value of this kmer with
 * that of another kmer.
 */
struct syncmer {
    kmer_type kmer;
    uint64_t position;

    syncmer(kmer_type kmer, uint64_t position) : kmer(kmer), position(position) {}

    bool operator<(const struct syncmer& other) const {
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
kmer_type encode(int* map, Iter begin, Iter end) {
    kmer_type res = 0;
    for (Iter it = begin; it < end; it++) {
        res *= 4;
        res |= map[static_cast<size_t>(*it)];
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
 * @brief   Identifies and inserts a syncmer into the syncmer vector if it meets certain conditions.
 *
 * This function scans through a k-mer to identify the lexicographically smallest s-mer (substring of size `smerSize`)
 * and compares its index with predefined indicex (`beginIndex`). If the smallest s-mer is found
 * at either of these positions, the corresponding syncmer is encoded and appended to the `syncmers` vector.
 *
 * @param smerIndex      The index to compare with the smallest s-mer's position.
 * @param begin          Iterator pointing to the beginning of the current k-mer.
 * @param end            Iterator pointing to the end of the current k-mer.
 * @param current_index  The index of the current k-mer in the sequence.
 * @param kmerSize       The size of the k-mer (number of characters). 
 * @param smerSize       The size of the s-mer (number of characters).
 * @param syncmers       A vector to store the resulting syncmers.
 * @param map            A pointer to an array used to encode s-mers.
 */
void process(Iter begin, Iter end, uint64_t current_index, int kmerSize, int smerSize, int smerIndex, Vec &syncmers, int* map) {
    Iter min_smer = begin;
    int min_smer_index = smerIndex;
    Iter cur_smer = begin+1;
    
    while (cur_smer+smerSize <= end) {   
        int temp_index = 0;

        while (temp_index < smerSize && 
                ((*(min_smer+temp_index))|0x20) == ((*(cur_smer+temp_index))|0x20))
            temp_index++;

        if (temp_index != smerSize && 
            ((*(min_smer+temp_index))|0x20) > (*((cur_smer+temp_index))|0x20)) {   
            min_smer = cur_smer;
            min_smer_index = cur_smer - begin;

            if (smerIndex < min_smer_index)
                return;
        }

        cur_smer++;
    }
    
    if (min_smer_index == smerIndex) {
        syncmers.emplace_back(encode(map, begin, begin+kmerSize), current_index);
    }
};

/**
 * @brief   Identifies and inserts a syncmer into the syncmer map if it meets certain conditions.
 *
 * This function scans through a k-mer to identify the lexicographically smallest s-mer (substring of size `smerSize`)
 * and compares its index with predefined indicex (`beginIndex`). If the smallest s-mer is found
 * at either of these positions, the corresponding syncmer is encoded and appended to the `syncmers` map where
 * key is the encoded value and value is a vector of positions of syncmers that appear in string.
 *
 * @param smerIndex      The index to compare with the smallest s-mer's position.
 * @param begin          Iterator pointing to the beginning of the current k-mer.
 * @param end            Iterator pointing to the end of the current k-mer.
 * @param current_index  The index of the current k-mer in the sequence.
 * @param kmerSize       The size of the k-mer (number of characters). 
 * @param smerSize       The size of the s-mer (number of characters).
 * @param map            A pointer to an array used to encode s-mers.
 * @param syncmerMap     A map tha stores syncmers and their locations in the sequence.
 */
void process3(Iter begin, Iter end, uint64_t current_index, int kmerSize, int smerSize, int smerIndex, int* map, std::map<kmer_type, std::vector<uint64_t>>& syncmerMap) {
    Iter min_smer = begin;
    int min_smer_index = smerIndex; // 0
    Iter cur_smer = begin+1;
    
    while (cur_smer+smerSize <= end) {   
        int temp_index = 0;

        while (temp_index < smerSize && 
                ((*(min_smer+temp_index))|0x20) == ((*(cur_smer+temp_index))|0x20))
            temp_index++;

        if (temp_index != smerSize && 
            ((*(min_smer+temp_index))|0x20) > (*((cur_smer+temp_index))|0x20)) {   
            min_smer = cur_smer;
            min_smer_index = cur_smer - begin;

            if (smerIndex < min_smer_index)
                return;
        }

        cur_smer++;
    }
    
    if (min_smer_index == smerIndex) {
        kmer_type kmer = encode(map, begin, begin+kmerSize);
        syncmerMap[kmer].push_back(current_index);
    }
};
