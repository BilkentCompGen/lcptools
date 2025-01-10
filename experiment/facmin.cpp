/**
 * @file    facmin.cpp
 * @brief   Locating Canonical Minimizers in Genome Sequences
 *
 * This program reads the genome sequence and locates the minimizers in each window, using the
 * predefined k-mer and window sizes. Statistical data such as the total number of minimizers in the genome,
 * the number of unique minimizers, the mean of the distances between consecutive minimizers and their
 * standard deviation are calculated.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm> // sort, to find distinct minimizers
#include <cassert>
#include "helper.cpp"
#include "sketches/minimizer.cpp"  // minimizer struct, map init, process2

#define CAPACITY 250000000

/**
 * @brief   Finds minimizers for each window in a genomic sequence and calculates processing time.
 *
 * This function scans through a given genomic sequence with windows of size `windowSize`.
 * For each window, it identifies the lexicographically smallest k-mer (minimizer) using the
 * `emplaceMinimizer` function and stores it in the `minimizers` vector. It also calculates the distances
 * between consecutive minimizers and records the time taken to process the sequence.
 *
 * @param gapSize        The suffix and prefix gaps size.
 * @param intraGapSize   The itra-gaps size.
 * @param sequence       The input sequence of nucleotides/characters to be processed.
 * @param minimizers     A vector to store the resulting minimizers.
 * @param kmerSize       The size of the k-mer (number of characters).
 * @param windowSize     The size of the sliding window in which the minimizers are searched.
 * @param map            A pointer to an array used for encoding k-mers.
 * @param distances      An array to store the frequencies of distances between consecutive minimizers.
 */
void findMinimizers(int &gapSize, int &intraGapSize, std::string &sequence, Vec &minimizers, int kmerSize, int windowSize, int *map, int *rc_map, int *distances) {
    int current_index = 0;
    
    for (Iter it = sequence.begin(); it < sequence.end() - windowSize - kmerSize; it++) {
        process2(it, it+windowSize, current_index, kmerSize, minimizers, map, rc_map);
        current_index++;
    }
    
    gapSize += minimizers.begin()->position;
    gapSize += sequence.size() - ((minimizers.end() - 1)->position+kmerSize); // size - (last.start+kmer)

    for (Vec::iterator it = minimizers.begin()+1; it < minimizers.end(); it++) {
        distances[ it->position - (it - 1)->position ]++;
        if ((it-1)->position + kmerSize < it->position) {
            intraGapSize += (it->position - ((it-1)->position + kmerSize));
        }
        assert(it->position - (it-1)->position <= windowSize); // two subsequent minimizers cannot have distance more than window size
    }

    std::cout << "Length of the processed sequence: " << format_int(sequence.size()) << 
    " minimizer count: " << format_int(minimizers.size()) << std::endl;
};

// Read and process the genome sequence and print the results
int main(int argc, char** argv) {
    
    // check if the correct number of arguments is provided
    if (argc < 4) {
        std::cerr << "Wrong format: " << argv[0] << " [infile] [kmer-size] [window-size]" << std::endl;
        return -1;
    }

    // open the FASTA file using the filename from the command-line argument
    std::ifstream input(argv[1]);
    if (!input.good()) {
        std::cerr << "Error opening: " << argv[1] << std::endl;
        return -1;
    }

    int kmer_size = atoi(argv[2]);
    int window_size = atoi(argv[3]);

    assert(window_size <= kmer_size);

    // variables
    std::string gen, line, id;
    int distances[window_size+1] = { 0 };
    int map[128];
    int rc_map[128];
    int gapSize = 0;
    int intraGapSize = 0;

    gen.reserve(CAPACITY);
    init_map(map);
    init_rc_map(rc_map);

    std::vector<Vec> minimizers;

    // Open genome file
    std::fstream genome;
    genome.open(argv[1], std::ios::in);
    
    // Read the file line by line
    if (genome.is_open()) {  
        
        std::cout << "Program begins" << std::endl;
        std::cout << "K-mer size: " << kmer_size << ", Window size: " << window_size << std::endl;

        while (getline(genome, line)) {

            if (line[0] == '>') {

                // Process previous chromosome before moving into new one
                if (gen.size() != 0) {
                    Vec sequence_minimizers;
                    sequence_minimizers.reserve(3 * gen.size() / window_size);
                    findMinimizers(gapSize, intraGapSize, gen, sequence_minimizers, kmer_size, window_size, map, rc_map, distances);
                    minimizers.push_back(sequence_minimizers);
                }

                id = line.substr(1);
                std::cout << "Processing started for " << id << std::endl;

                gen.clear();

                continue;    
            }
            else if (line[0] != '>'){
                gen += line;
            }
        }
        
        // Process last chromosome before calculating stats
        if (gen.size() != 0) {
            Vec sequence_minimizers;
            sequence_minimizers.reserve(3 * gen.size() / window_size);
            findMinimizers(gapSize, intraGapSize, gen, sequence_minimizers, kmer_size, window_size, map, rc_map, distances);
            minimizers.push_back(sequence_minimizers);
        }

        genome.close();
    }
    
    // Count total number of minimizers
    int numOfMinimizers = 0;
    for (size_t i = 0; i < minimizers.size(); i++) {
        numOfMinimizers += minimizers[i].size();   
    }

    // Count distinct minimizers
    int numOfDistintMinimizers = 0;
    std::cout << "Counting distinct minimizers..." << std::endl;
    std::vector<kmer_type> flattened_minimizers;
    flattened_minimizers.reserve(numOfMinimizers);

    for (const Vec &sequence_minimizers : minimizers) {
        for (const minimizer &min : sequence_minimizers) {
            flattened_minimizers.push_back(min.kmer);
        }
    }

    std::sort(flattened_minimizers.begin(), flattened_minimizers.end());

    numOfDistintMinimizers = 1;
    for (std::vector<kmer_type>::iterator it = flattened_minimizers.begin() + 1; it < flattened_minimizers.end(); it++) {
        if (*(it-1) != *(it)) {
            numOfDistintMinimizers++;
        }
    }

    std::cout << "Calculating stats..." << std::endl;

    double total_size = 0;
    for (int i=0; i<=window_size; i++) {
        total_size += (i * distances[i]);
    }
    total_size += gapSize;
    total_size += minimizers.size() * kmer_size;
    
    // Calculate stats
    double average = mean(distances, window_size+1);
    double std_dev = stdev(distances, window_size+1, average);
    
    // Output the results
    std::cout << "Total k-mers: " << format_int(numOfMinimizers) << std::endl;
    std::cout << "Unique k-mers: " << format_int(numOfDistintMinimizers) << std::endl;
    std::cout << "Avg Dist. : " << format_double(average) << std::endl;
    std::cout << "StdDev Dist. : " << format_double(std_dev) << std::endl;
    std::cout << "Gap size: " << format_int(gapSize) << std::endl;
    std::cout << "Total size: " << format_double(total_size) << std::endl;
};
