/**
 * @file    falcp-mem.cpp
 * @brief   Analysis and Processing of Genomic Data
 *
 * This program is designed for in-depth analysis of genomic sequences. It reads
 * genomic data, processes it through multiple levels of analysis (defined by
 * LCP_LEVEL), and computes execution times of the processing of different LCP levels.
 *
 * The program leverages a series of custom functions and data structures, operating
 * on large volumes of data with efficiency and accuracy in focus.
 */

#include "helper.cpp"
#include "lps.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#define STRING_SIZE 250000000

/**
 * @brief Processes genomic sequence data, tracks execution time, and analyzes
 * various statistics.
 *
 * This function processes a given genomic sequence by performing multiple levels
 * of LCP (Locally Consistent Parsing) analysis. It tracks the execution time for
 * each level of processing.
 *
 * @param sequence The genomic sequence (string) to be analyzed.
 * @param total_core_counts An array storing the number of LCP cores found at each level.
 * @param durations A vector storing the durations (in milliseconds) of each level's
 *                  processing time.
 */
void process(std::string &sequence, 
             size_t (&total_core_counts)[LCP_LEVEL], 
             std::vector<std::chrono::milliseconds> &durations) {

	auto start = std::chrono::high_resolution_clock::now();
    struct lps str;
    init_lps(&str, sequence.c_str(), sequence.size());

	auto extraction_end = std::chrono::high_resolution_clock::now();
	total_core_counts[0] += str.size;

	durations[0] += std::chrono::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(extraction_end - start).count());

	for (int i = 1; i < LCP_LEVEL; i++) {

		auto start_level = std::chrono::high_resolution_clock::now();

		lps_deepen1_parallel(&str, 4);

		auto stop_level = std::chrono::high_resolution_clock::now();
		total_core_counts[i] += str.size;

		durations[i] += std::chrono::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(stop_level - start_level).count());
	}

	std::cout << "Length of the processed sequence: " << format_int(sequence.size()) << std::endl;

    free_lps(&str);
	sequence.clear();
}

/**
 * @brief The entry point of the program.
 *
 * The main function coordinates the entire genomic data analysis process. It
 * initializes necessary data structures, reads input genomic sequences.
 *
 */
int main(int argc, char **argv) {

	if (argc < 2) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] " << std::endl;
		return -1;
	}

	std::ifstream input(argv[1]);
	if (!input.good()) {
		std::cerr << "Error opening: " << argv[1] << " . You have failed." << std::endl;
		return -1;
	}

	// variables
	std::string line;

	std::fstream genome;
	genome.open(argv[1], std::ios::in);

	std::vector<std::chrono::milliseconds> durations(LCP_LEVEL);
	size_t total_core_counts[LCP_LEVEL] = {0, 0, 0, 0, 0, 0, 0, 0};

	// read file
	if (genome.is_open()) {

		std::string sequence, id;
		sequence.reserve(STRING_SIZE);

		// initializing coefficients of the alphabet and hash tables
        LCP_INIT();

		std::cout << "Program begins" << std::endl;

		while (getline(genome, line)) {

			if (line[0] == '>') {

				// process previous chromosome before moving into new one
				if (sequence.size() != 0) {
					process(sequence, total_core_counts, durations);
				}

				id = line.substr(1);
				std::cout << "Processing started for " << id << std::endl;
				continue;

			} else if (line[0] != '>') {
				sequence += line;
			}
		}

		if (sequence.size() != 0) {
			process(sequence, total_core_counts, durations);
		}

		genome.close();
	}

	std::cout << std::endl;

	std::string sep = " & ";

	std::cout << "LCP level";
	for (int i = 0; i < LCP_LEVEL; i++) {
		std::cout << sep << i + 1;
	}
	std::cout << std::endl;

	// Total Cores
	std::cout << "Total \\# of Cores";
	for (int i = 0; i < LCP_LEVEL; i++) {
		std::cout << sep << format_int(total_core_counts[i]);
	}
	std::cout << " \\\\" << std::endl;

	// Execution Time
	std::cout << "Exec. Time (sec)";
	for (int i = 0; i < LCP_LEVEL; i++) {
		std::cout << sep << format_double(((double)durations[i].count()) / 1000);
	}
	std::cout << " \\\\" << std::endl << std::endl;

	return 0;
};
