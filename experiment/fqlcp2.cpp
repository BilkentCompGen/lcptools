/**
 * @file    fqlcp2.cpp
 * @brief   Locating LCP cores in Genome Sequences
 * 
 */

#include <thread>
#include <mutex>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iomanip>
#include <cmath>
#include <cassert>
#include "lps.h"
#include "helper.cpp"

#define PASS_NUMBER 7
#define LCP_LEVEL_MIN 2
#define LCP_LEVEL_MAX 6
#define LCP_LEVEL_COUNT 5

std::mutex mtx;

typedef uint32_t kmer_type;
typedef std::string::iterator Iter;

void calculate_metrics(stats_type results[LCP_LEVEL_COUNT][3]) {
    std::cout << "LCP level" << std::endl;

    for (int i = LCP_LEVEL_MIN; i <= LCP_LEVEL_MAX; i++) {
        std::cout << i;

        stats_type TP = results[i-LCP_LEVEL_MIN][0]; // True Positive
        stats_type FP = results[i-LCP_LEVEL_MIN][1]; // False Positive
        stats_type FN = results[i-LCP_LEVEL_MIN][2]; // False Negative

        // calculate Precision
        double precision = (TP + FP == 0) ? 0.0 : static_cast<double>(TP) / (TP + FP);

        // calculate Sensitivity
        double sensitivity = (TP + FN == 0) ? 0.0 : static_cast<double>(TP) / (TP + FN);
        
        std::cout << " & " << std::fixed << std::setprecision(5) << precision;
        std::cout << "," << std::fixed << std::setprecision(5) << sensitivity << std::endl;
    }
};

void findLcpCores(std::string &sequence, bool rc, int lcpLevel, std::map<kmer_type, std::vector<uint64_t>>& lcpCoresMap) {
    struct lps str;

    if (rc) {
        init_lps2(&str, sequence.c_str(), sequence.size());
    } else {
        init_lps(&str, sequence.c_str(), sequence.size());
    }  
    lps_deepen(&str, lcpLevel);

    for(int i=0; i<str.size; i++) {
        lcpCoresMap[str.cores[i].label].push_back(str.cores[i].start);
    }
    free_lps(&str);
};

void t_process(int thread_index, const char *fa, const char *mf, const char *ff, int lcpLevel, stats_type (&results)[3], stats_type (&stats)[5]) {

    {
        std::lock_guard<std::mutex> lock(mtx);
        std::cout << "Thread " << thread_index << " started processing from level=" << lcpLevel << std::endl;
    }

    // open the fasta file
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

    LCP_INIT();

    std::string fa_line;
    std::string fq_id, maf_id;
    std::string fq_line, maf_line;
    std::string maf_sign;

    std::map<kmer_type, std::vector<uint64_t>> mapReference;
    std::string sequence;
    sequence.reserve(250000000);

    while (getline(fastaFile, fa_line)) {
        if (fa_line[0] != '>') {
            sequence += fa_line;
        }
    }

    findLcpCores(sequence, false, lcpLevel, mapReference);
    fastaFile.close();
    sequence.clear();
    
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

    stats_type true_positive = 0; // A lcp core that is present in both the original read and the simulated read.
    stats_type false_positive = 0; // A lcp core that is present in the simulated read but not in the original read.
    stats_type false_negative = 0; // A lcp core that is present in the original read but not in the simulated read.

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
        
        // now we have fq_line, index, size and sign of read
        Iter it_prev = sequence.begin();
        for (Iter it_cur = sequence.begin(); it_cur < sequence.end(); it_cur++) {
            if (*it_cur != '-') {  
                *it_prev = *it_cur;
                it_prev++;
            }
        }
        sequence.erase(it_prev, sequence.end());

        // this is where the fun begins
        std::map<kmer_type, std::vector<uint64_t>> mapGTRead;
        findLcpCores(sequence, false, lcpLevel, mapGTRead);

        std::map<kmer_type, std::vector<uint64_t>> mapSimRead;
        findLcpCores(fq_line, maf_sign == "-", lcpLevel, mapSimRead);

        if (maf_sign == "-") {
            stats[0]++;
        } else {
            stats[1]++;
        }

        // Algorithm to find how many lcp cores match
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

    stats_type results[LCP_LEVEL_COUNT][3];
    stats_type stats[LCP_LEVEL_COUNT][5];

    for (int i=0; i<LCP_LEVEL_COUNT; i++) {
        for (int k=0; k<3; k++) {
            results[i][k] = 0;
        }
        for (int k=0; k<5; k++) {
            stats[i][k] = 0;
        }
    }

    std::cout << "Program begins..." << std::endl;

    std::thread threads[LCP_LEVEL_COUNT];
    int thread_id = 0;

    for (int i=0; i<LCP_LEVEL_COUNT; i++) {
        threads[i] = std::thread(t_process, thread_id, argv[1], argv[2], argv[3], i+LCP_LEVEL_MIN, std::ref(results[i]), std::ref(stats[i]));
        thread_id++;
    }

    for (int i=0; i<LCP_LEVEL_COUNT; i++) {
        threads[i].join();
    }

    for (int i=0; i<LCP_LEVEL_COUNT; i++) {
        std::cout << "l=" << i+LCP_LEVEL_MIN 
                    << " TP: " << results[i][0] 
                    << " FP: " << results[i][1]
                    << " FN: " << results[i][2] 
                    << " fwd_reads: " << stats[i][0]
                    << " rc_read: " << stats[i][1]
                    << " max_count: " << stats[i][2]
                    << " exceeding #: " << stats[i][3]
                    << " total sketches: " << stats[i][4] << std::endl;
    }

    calculate_metrics(results);

    return 0;
};
