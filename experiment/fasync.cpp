/**
 * @file    fasync.cpp
 * @brief   Locating Syncmers in Genome Sequences
 *
 * This program reads the genome sequences and locates the syncmers in it by
 * distinguishing the kmers that have their lexicographically smallest smers starting/
 * ending in the chosen indices. Statistical data such as the total number of syncmers,
 * the number of unique syncmers, the mean of the distances between consecutive
 * syncmers and their standard deviation are calculated.
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "helper.cpp"
#include "sketches/syncmer.cpp"

#define CAPACITY            250000000
#define DISTANCE_ARRAY_SIZE 9999

/**
 * @brief   Finds syncmers in a genomic sequence and calculates processing time.
 *
 * This function processes a genomic sequence by sliding a window of size `kmerSize` across it, calling
 * `emplaceSyncmer` to identify syncmers based on predefined index (`beginIndex`). The
 * function also calculates distances between consecutive syncmers and measures the time taken to process
 * the entire sequence.
 *
 * @param gapSize        The suffix and prefix gaps size.
 * @param intraGapSize   The itra-gaps size.
 * @param sequence       The input sequence of nucleotides/characters to be processed.
 * @param syncmers       A vector to store the resulting syncmers.
 * @param kmerSize       The size of the k-mer (number of characters).
 * @param smerSize       The size of the s-mer (number of characters).
 * @param smerIndex      The index to compare with the smallest s-mer's position.
 * @param map            A pointer to an array used for encoding k-mers.
 * @param distances      An array to store the frequencies of distances between consecutive syncmers.
 */
void findSyncmers(int &gapSize, int &intraGapSize, std::string &sequence, Vec &syncmers, int kmerSize, int smerSize, int smerIndex, int* map, int* distances) {
    uint64_t current_index = 0;

    for (Iter it = sequence.begin(); it <= sequence.end() - kmerSize; ++it) {
        process(it, it + kmerSize, current_index, kmerSize, smerSize, smerIndex, syncmers, map);
        current_index++;
    }

    gapSize += syncmers.begin()->position;
    gapSize += sequence.size() - ((syncmers.end() - 1)->position + kmerSize);

    for (Vec::iterator it = syncmers.begin()+1; it < syncmers.end(); it++) {
        distances[ it->position - (it - 1)->position ]++;
        if ((it-1)->position + kmerSize < it->position) {
            intraGapSize += (it->position - ((it-1)->position + kmerSize));
        }
    }

    std::cout << "Length of the processed sequence: " << format_int(sequence.size()) << 
    " syncmer count: " << format_int(syncmers.size()) << std::endl;
};

// Read and process the sequence and print the results
int main(int argc, char** argv) {
    
    // Check if the correct number of arguments is provided
    if (argc < 5) {
        std::cerr << "Wrong format: " << argv[0] << " [infile] [kmer-size] [smer-size] [smer-index]" << std::endl;
        return -1;
    }

    // Open the FASTA file using the filename from the command-line argument
    std::ifstream input(argv[1]);
    if (!input.good()) {
        std::cerr << "Error opening: " << argv[1] << ". You have failed." << std::endl;
        return -1;
    }

    int kmer_size = atoi(argv[2]);
    int smer_size = atoi(argv[3]);
    int smer_index = atoi(argv[4]);

    // Variables
    std::string gen, line, id;
    int distances[DISTANCE_ARRAY_SIZE];
    for (int i=0; i<DISTANCE_ARRAY_SIZE; i++) {
        distances[i] = 0;
    }
    int map[128];
    int gapSize = 0;
    int intraGapSize = 0;

    gen.reserve(CAPACITY);
    init_map(map);
    
    std::vector<Vec> syncmers;

    // Open genome file
    std::fstream genome;
    genome.open(argv[1], std::ios::in);

    // Read the file line by line
    if (genome.is_open()) {  
        
        std::cout << "Program begins" << std::endl;
        std::cout << "K-mer size: " << kmer_size << " S-mer size: " << smer_size << " S-mer index: " << smer_index << std::endl;

        while (getline(genome, line)) {

            if (line[0] == '>') {

                // Process previous chromosome before moving into new one
                if (gen.size() != 0) {
                    Vec sequence_syncmers;
                    sequence_syncmers.reserve(gen.size());
                    findSyncmers(gapSize, intraGapSize, gen, sequence_syncmers, kmer_size, smer_size, smer_index, map, distances);
                    syncmers.push_back(sequence_syncmers);
                }

                id = line.substr(1);
                std::cout << "Processing started for " << id << std::endl;

                // OMG!!
                gen.clear();

                continue;    
            }
            else if (line[0] != '>'){
                gen += line;
            }
        }

        // Process last chromosome before calculating stats
        if (gen.size() != 0) {
            Vec sequence_syncmers;
            sequence_syncmers.reserve(gen.size());
            findSyncmers(gapSize, intraGapSize, gen, sequence_syncmers, kmer_size, smer_size, smer_index, map, distances);
            syncmers.push_back(sequence_syncmers);
        }

        genome.close();
    }
    
    // Count total number of Sycnmers
    int numOfSyncmers = 0;
    for (size_t i = 0; i < syncmers.size(); i++) {
        numOfSyncmers += syncmers[i].size();   
    }

    // Count distinct syncmers
    int numOfDistintSyncmers = 0;
    std::cout << "Counting distinct syncmers..." << std::endl;
    std::vector<kmer_type> flattened_syncmers;
    flattened_syncmers.reserve(numOfSyncmers);

    for (const Vec &sequence_syncmers : syncmers) {
        for (const struct syncmer &syncmer : sequence_syncmers) {
            flattened_syncmers.push_back(syncmer.kmer);
        }
    }

    std::sort(flattened_syncmers.begin(), flattened_syncmers.end());

    numOfDistintSyncmers = 1;
    for (std::vector<kmer_type>::iterator it = flattened_syncmers.begin() + 1; it < flattened_syncmers.end(); it++) {
        if (*(it-1) != *(it)) {
            numOfDistintSyncmers++;
        }
    }
    
    std::cout << "Calculating stats..." << std::endl;

    double total_size = 0;
    for (int i=0; i<DISTANCE_ARRAY_SIZE; i++) {
        total_size += (i * distances[i]);
    }
    total_size += gapSize;
    total_size += syncmers.size() * kmer_size;

    // Calculate stats
    double average = mean(distances, DISTANCE_ARRAY_SIZE);
    double std_dev = stdev(distances, DISTANCE_ARRAY_SIZE, average);

    // Output the results
    std::cout << "Total k-mers: " << format_int(numOfSyncmers) << std::endl;
    std::cout << "Unique k-mers: " << format_int(numOfDistintSyncmers) << std::endl;
    std::cout << "Avg Dist. : " << format_double(average) << std::endl;
    std::cout << "StdDev Dist. : " << format_double(std_dev) << std::endl;
    std::cout << "Gap size: " << format_double(gapSize+intraGapSize) << std::endl;
    std::cout << "Total size: " << format_double(total_size) << std::endl;
};