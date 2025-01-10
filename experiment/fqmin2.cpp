/**
 * @file    fqmin2.cpp
 * @brief   Locating Minimizers in Genome Sequences
 * 
 */

#include <thread>
#include <mutex>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cassert>
#include "helper.cpp"
#include "sketches/minimizer.cpp"

#define PASS_NUMBER 7
#define KMER_VALUES_SIZE 22
#define WINDOW_VALUES_SIZE 4

std::mutex mtx;

void calculate_metrics(int kmer_size_values[KMER_VALUES_SIZE], int window_size_values[WINDOW_VALUES_SIZE], stats_type results[KMER_VALUES_SIZE][WINDOW_VALUES_SIZE][3]) {
    for (int i = 0; i < KMER_VALUES_SIZE; i++) {
        std::cout << " & " << kmer_size_values[i];
    }
    std::cout << " \\\\" << std::endl;

    for (int j = 0; j < WINDOW_VALUES_SIZE; j++) {
        std::cout << window_size_values[j];

        for (int i = 0; i < KMER_VALUES_SIZE; i++) {
            stats_type TP = results[i][j][0]; // True Positive
            stats_type FP = results[i][j][1]; // False Positive
            stats_type FN = results[i][j][2]; // False Negative

            // calculate Precision
            double precision = (TP + FP == 0) ? 0.0 : static_cast<double>(TP) / (TP + FP);

            // calculate Sensitivity
            double sensitivity = (TP + FN == 0) ? 0.0 : static_cast<double>(TP) / (TP + FN);

            std::cout << " & " << std::fixed << std::setprecision(5) << precision;
            std::cout << "," << std::fixed << std::setprecision(5) << sensitivity;
        }
        std::cout << " \\\\" << std::endl; 
    }
};

void findMinimizers(std::string &sequence, int kmerSize, int windowSize, int *map, std::map<kmer_type, std::vector<uint64_t>>& minimizerMap) {
    uint64_t current_index = 0;
    int64_t previous_index = -1;

    for (Iter it = sequence.begin(); it < sequence.end() - windowSize - kmerSize; it++) {
        previous_index = process3(it, it+windowSize, previous_index, current_index, kmerSize, map, minimizerMap);
        current_index++;
    }
};

void t_process(int thread_index, const char *fa, const char *mf, const char *ff, int kmer_size, int window_size, stats_type (&results)[3], stats_type (&stats)[5]) {

    {
        std::lock_guard<std::mutex> lock(mtx);
        std::cout << "Thread " << thread_index << " started processing from k=" << kmer_size << " w= " << window_size << std::endl;
    }

    // open the fa file
    std::ifstream fastaFile(fa);
    if (!fastaFile.good()) {
        std::cerr << "Error opening: " << fa << std::endl;
        return;
    }

    // open the maf file
    std::ifstream mafFile(mf);
    if (!mafFile.good()) {
        std::cerr << "Error opening: " << mf << std::endl;
        return;
    }

    // open the fastq file
    std::ifstream fastqFile(ff);
    if (!fastqFile.good()) {
        std::cerr << "Error opening: " << ff << std::endl;
        return;
    }

    int map[128] = { 0 };
    init_map(map);

    std::string fa_line;
    std::string fq_id, maf_id;
    std::string fq_line, maf_line;
    std::string maf_sign;

    std::string sequence;
    sequence.reserve(250000000);

    while (getline(fastaFile, fa_line)) {
        if (fa_line[0] != '>') {
            sequence += fa_line;
        }
    }

    std::map<kmer_type, std::vector<uint64_t>> mapReference;
    findMinimizers(sequence, kmer_size, window_size, map, mapReference);

    fastaFile.close();

    stats[2] = 0;
    stats[3] = 0;
    stats[4] = 0;
    for (const std::pair<const kmer_type, std::vector<uint64_t>> &sim_pair : mapReference) {
        if (stats[2] < sim_pair.second.size()) 
            stats[2] = sim_pair.second.size();
        if (MAXIMUM_FREQ_THRESHOLD <= sim_pair.second.size())
            stats[3]++;
        stats[4] += sim_pair.second.size();
    }

    stats_type true_positive = 0; // A minimizer that is present in both the original read and the simulated read.
    stats_type false_positive = 0; // A minimizer that is present in the simulated read but not in the original read.
    stats_type false_negative = 0; // A minimizer that is present in the original read but not in the simulated read.

    while (true) {
        // read fastq read. if there is a read in there, then there must be simulated reads in maf
        if (! getline(fastqFile, fq_line)) { // reads first line of read: @ID
            break;
        }

        fq_id = fq_line.substr(1, fq_line.rfind('/') - 1);
        getline(fastqFile, fq_line); // fq_line contains high quality read now

        // process maf file to find the reference sequence
        while (true) {
            getline(mafFile, maf_line); // reads first line of read: a
            getline(mafFile, maf_line); // read ref             

            // parse the line
            std::istringstream issRef(maf_line);
            issRef >> sequence >> sequence >> sequence >> sequence >> sequence >> sequence >> sequence;

            getline(mafFile, maf_line); // read simulated read line

            // parse the line to get id
            std::istringstream issId(maf_line);
            issId >> maf_id >> maf_id >> maf_sign >> maf_sign >> maf_sign;
            maf_id = maf_id.substr(0, maf_id.rfind('/'));

            getline(mafFile, maf_line); // skip empty line (last line)

            if (maf_id == fq_id) { // found id match
                break;
            }
            
            for(int i=1; i<PASS_NUMBER; i++) { // skip other simulated reads
                getline(mafFile, maf_line);
                getline(mafFile, maf_line);
                getline(mafFile, maf_line);
                getline(mafFile, maf_line);
            }
        }

        // remove all alignment information ('-') from ref
        Iter it_prev = sequence.begin();
        Iter it_cur = sequence.begin();
        for (; it_cur < sequence.end(); it_cur++) {
            if (*it_cur != '-') {  
                *it_prev = *it_cur;
                it_prev++;
            }
        }
        sequence.erase(it_prev, sequence.end()); 
        
        // now we have sequence, fq_line, and sign
        
        if (maf_sign == "-") {
            reverse_complement(fq_line);
            stats[0]++;
        } else {
            stats[1]++;
        }

        // this is where the fun begins
        std::map<kmer_type, std::vector<uint64_t>> mapGTRead;
        findMinimizers(sequence, kmer_size, window_size, map, mapGTRead);

        std::map<kmer_type, std::vector<uint64_t>> mapSimRead;
        findMinimizers(fq_line, kmer_size, window_size, map, mapSimRead);

        // Algorithm to find how many minimizers match
        alignment_global_stats<kmer_type, uint64_t>(mapReference, mapGTRead, mapSimRead, true_positive, false_positive, false_negative);

        // DONE

        // MOVE FASTQ FILE
        getline(fastqFile, fq_line); // move +
        getline(fastqFile, fq_line); // move quality score line

        // MOVE MAF FILE
        for(int i=1; i<PASS_NUMBER; i++) { // skip other simulated reads
            getline(mafFile, maf_line);
            getline(mafFile, maf_line);
            getline(mafFile, maf_line);
            getline(mafFile, maf_line);
        }
    }

    results[0] = true_positive;
    results[1] = false_positive;
    results[2] = false_negative;

    mafFile.close();
    fastqFile.close();
};

// Read and process the genome sequence and print the results
int main(int argc, char **argv) {
    // check if the correct number of arguments is provided
    if (argc < 4) {
        std::cerr << "Wrong format: " << argv[0] << " [fa-file] [maf-file] [fq-file]" << std::endl;
        return -1;
    }

    int kmer_size_values[KMER_VALUES_SIZE] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
    int window_size_values[WINDOW_VALUES_SIZE] = {10, 11, 15, 19};

    stats_type results[KMER_VALUES_SIZE][WINDOW_VALUES_SIZE][3];
    stats_type stats[KMER_VALUES_SIZE][WINDOW_VALUES_SIZE][5];

    for (int i=0; i<KMER_VALUES_SIZE; i++) {
        for (int j=0; j<WINDOW_VALUES_SIZE; j++) {
            for (int k=0; k<3; k++) {
                results[i][j][k] = 0;
            }
            for (int k=0; k<5; k++) {
                stats[i][j][k] = 0;
            }
        }
    }

    std::cout << "Program begins..." << std::endl;

    std::thread threads[KMER_VALUES_SIZE][WINDOW_VALUES_SIZE];
    int thread_id = 0;

    for (int i=0; i<KMER_VALUES_SIZE; i++) {
        for (int j=0; j<WINDOW_VALUES_SIZE; j++) {
            threads[i][j] = std::thread(t_process, thread_id, argv[1], argv[2], argv[3], kmer_size_values[i], window_size_values[j], std::ref(results[i][j]), std::ref(stats[i][j]));
            thread_id++;
        }
    }

    for (int i=0; i<KMER_VALUES_SIZE; i++) {
        for (int j=0; j<WINDOW_VALUES_SIZE; j++) {
            threads[i][j].join();
        }
    }

    for (int i=0; i<KMER_VALUES_SIZE; i++) {
        for (int j=0; j<WINDOW_VALUES_SIZE; j++) {
            int kmer_size = kmer_size_values[i];
            int window_size = window_size_values[j];

            std::cout << "k=" << kmer_size 
                      << " W=" << window_size 
                      << " TP: " << results[i][j][0] 
                      << " FP: " << results[i][j][1]
                      << " FN: " << results[i][j][2]
                      << " fwd_reads: " << stats[i][j][0]
                      << " rc_read: " << stats[i][j][1]
                      << " max_count: " << stats[i][j][2]
                      << " exceeding #: " << stats[i][j][3]
                      << " total sketches: " << stats[i][j][4] << std::endl;
        }
    }

    calculate_metrics(kmer_size_values, window_size_values, results);

    return 0;
};
