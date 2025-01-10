/**
 * @file    fasync.cpp
 * @brief   Locating Syncmers in Genome Sequences
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
#include "sketches/syncmer.cpp"

#define PASS_NUMBER 7
#define KMER_VALUES_SIZE 4
#define MAX_SMER_SIZE 19
#define SMER_BEGIN_INDEX 0

std::mutex mtx;

void calculate_metrics(int kmer_size_values[KMER_VALUES_SIZE], stats_type results[KMER_VALUES_SIZE][MAX_SMER_SIZE][4]) {
    
    for (int i = 0; i < KMER_VALUES_SIZE; i++) {
        if (kmer_size_values[i] > MAX_SMER_SIZE) {
            std::cerr << "Error: kmer_size_values[" << i << "] exceeds MAX_SMER_SIZE." << std::endl;
            std::cout << "i=" << i << " val=" << kmer_size_values[i] << std::endl;
            continue;
        }

        std::cout << kmer_size_values[i];

        for (int j = 2; j < kmer_size_values[i]; j++) {
            stats_type TP = results[i][j][0]; // True Positive
            stats_type FP = results[i][j][1]; // False Positive
            stats_type FN = results[i][j][2]; // False Negative

            // calculate Precision
            double precision = (TP + FP == 0) ? 0.0 : static_cast<double>(TP) / (TP + FP);

            // calculate Sensitivity
            double sensitivity = (TP + FN == 0) ? 0.0 : static_cast<double>(TP) / (TP + FN);
                   
            std::cout << " & " << std::fixed << std::setprecision(5) << precision;
            std::cout << "," << std::fixed << std::setprecision(5) << sensitivity; 

            assert(TP+FP+FN == results[i][j][3]);
        }
        std::cout << " \\\\" << std::endl;
    }
};

void findSyncmers(std::string &sequence, int kmerSize, int smerSize, int smerIndex, int* map, std::map<kmer_type, std::vector<uint64_t>>& syncmerMap) {
    uint64_t current_index = 0;

    for (Iter it = sequence.begin(); it <= sequence.end() - kmerSize; ++it) {
        process3(it, it + kmerSize, current_index, kmerSize, smerSize, smerIndex, map, syncmerMap);
        current_index++;
    }
};

void t_process(int thread_index, const char *mf, const char *ff, int kmer_size, int smer_size, stats_type &tp, stats_type &fp, stats_type &fn, stats_type &total, stats_type &fwd_reads, stats_type &rc_reads) {

    {
        std::lock_guard<std::mutex> lock(mtx);
        std::cout << "Thread " << thread_index << " started processing from k=" << kmer_size << " s= " << smer_size << std::endl;
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

    std::string fq_id, maf_id;
    std::string fq_line, maf_line;
    std::string maf_sign;
    std::string sequence;

    stats_type true_positive = 0; // A syncmer that is present in both the original read and the simulated read.
    stats_type false_positive = 0; // A syncmer that is present in the simulated read but not in the original read.
    stats_type false_negative = 0; // A syncmer that is present in the original read but not in the simulated read.
    stats_type total_syncmers = 0;

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
            rc_reads++;
        } else {
            fwd_reads++;
        }

        // this is where the fun begins
        std::map<kmer_type, std::vector<uint64_t>> mapGTRead;
        findSyncmers(sequence, kmer_size, smer_size, SMER_BEGIN_INDEX, map, mapGTRead);

        std::map<kmer_type, std::vector<uint64_t>> mapSimRead;
        findSyncmers(fq_line, kmer_size, smer_size, SMER_BEGIN_INDEX, map, mapSimRead);

        // Algorithm to find how many syncmers match
        alignment_pairwise_stats<kmer_type, uint64_t>(mapGTRead, mapSimRead, true_positive, false_positive, false_negative, total_syncmers);

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

    tp = true_positive;
    fp = false_positive;
    fn = false_negative;
    total = total_syncmers;
};

// Read and process the genome sequence and print the results
int main(int argc, char **argv) {

    if (argc < 3) {
        std::cerr << "Wrong format: " << argv[0] << " [maf-file] [fq-file]" << std::endl;
        return -1;
    }

    int kmer_size_values[KMER_VALUES_SIZE] = {10, 11, 15, 19};

    stats_type results[KMER_VALUES_SIZE][MAX_SMER_SIZE][4];
    stats_type reads[KMER_VALUES_SIZE][MAX_SMER_SIZE][2];

    for (int i=0; i<KMER_VALUES_SIZE; i++) {
        for (int j=0; j<MAX_SMER_SIZE; j++) {
            for (int k=0; k<4; k++) {
                results[i][j][k] = 0;
            }
            for (int k=0; k<2; k++) {
                reads[i][j][k] = 0;
            }
        }
    }

    std::cout << "Program begins..." << std::endl;

    std::thread threads[KMER_VALUES_SIZE][MAX_SMER_SIZE];
    int thread_id = 0;

    for (int i=0; i<KMER_VALUES_SIZE; i++) {
        for (int j=2; j<kmer_size_values[i]; j++) {
            threads[i][j] = std::thread(t_process, thread_id, argv[1], argv[2], kmer_size_values[i], j, std::ref(results[i][j][0]), std::ref(results[i][j][1]), std::ref(results[i][j][2]), std::ref(results[i][j][3]), std::ref(reads[i][j][0]), std::ref(reads[i][j][1]));
            thread_id++;
        }
    }

    for (int i=0; i<KMER_VALUES_SIZE; i++) {
        for (int j=2; j<kmer_size_values[i]; j++) {
            threads[i][j].join();
        }
    }

    for (int i=0; i<KMER_VALUES_SIZE; i++) {
        for (int j=2; j<kmer_size_values[i]; j++) {
            int kmer_size = kmer_size_values[i];
            int smer_size = j;

            std::cout << "k=" << kmer_size 
                      << " s=" << smer_size 
                      << " TP: " << results[i][j][0] 
                      << " FP: " << results[i][j][1]
                      << " FN: " << results[i][j][2] 
                      << " Total: " << results[i][j][3]
                      << " fwd_reads: " << reads[i][j][0]
                      << " rc_read: " << reads[i][j][1] << std::endl;
        }
    }
    
    calculate_metrics(kmer_size_values, results);

    return 0;
    
};
