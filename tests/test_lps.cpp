#include "core.h"
#include "lps.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

void log(const std::string &message) {
	std::cout << message << std::endl;
}

void test_lps_split_init() {

    LCP_INIT();

	std::ifstream genome("data/test.fasta");
    std::string sequence, line;

    getline(genome, line); // skip first header line

    while (getline(genome, line)) {
        if (line[0] != '>') {
            sequence += line;
        } else {
            break;
        }
    }
    genome.close();

    struct lps lps_obj1;
    init_lps(&lps_obj1, sequence.c_str(), sequence.size());
    lps_deepen(&lps_obj1, 7);

    struct lps lps_obj2;
    init_lps4(&lps_obj2, sequence.c_str(), sequence.size(), 7, 100000);

    assert(lps_eq(&lps_obj1, &lps_obj2) && "LCP split and merge result should be same as processing linearly");

    free_lps(&lps_obj1);
    free_lps(&lps_obj2);

    log("...  test_lps_split_init passed!");
}

void test_lps_file_io() {

	LCP_INIT();

	std::string test_string = "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT";
	struct lps lps_obj;
    init_lps(&lps_obj, test_string.c_str(), test_string.size());

	// write to file
	FILE *out = fopen("lps_test.bin", "wb");
    if (!out) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    write_lps(&lps_obj, out);

    fclose(out);

	// read from file
	FILE *in = fopen("lps_test.bin", "rb");
    struct lps lps_obj_from_file;
	init_lps3(&lps_obj_from_file, in);
	fclose(in);

	// compare the read object with the original
	assert(lps_obj.level == lps_obj_from_file.level && "Level should match after reading from file");
	assert(lps_obj.size == lps_obj_from_file.size && "Core size should match after reading from file");
	for (int i = 0; i < lps_obj.size; ++i) {
		assert(core_eq(&(lps_obj.cores[i]), &(lps_obj_from_file.cores[i])) && "Cores should match after reading from file");
	}

	// clean up the test file
	std::remove("lps_test.bin");

    free_lps(&lps_obj);
    free_lps(&lps_obj_from_file);

	log("...  test_lps_file_io passed!");
}

void test_lps_consistency() {

    std::ifstream genome("data/test.fasta");
    
    LCP_INIT();

    std::string sequence;
    std::string line;

    getline(genome, line); // skip first header line

    while (getline(genome, line)) {
        if (line[0] != '>') {
            sequence += line;
        } else {
            break;
        }
    }
    genome.close();

    struct lps lps_obj;
    init_lps(&lps_obj, sequence.c_str(), sequence.size());
    lps_deepen(&lps_obj, 5);

    int start = lps_obj.cores[500].start;
    int end = lps_obj.cores[500].end;

    // check if core is identified in the given intervals
    std::string subsequence1(sequence.begin()+start, sequence.begin()+end);

    struct lps lps_obj1;
    init_lps_offset(&lps_obj1, subsequence1.c_str(), subsequence1.size(), start);
    lps_deepen(&lps_obj1, 5);

    assert(core_eq(&(lps_obj.cores[500]), &(lps_obj1.cores[0])) && "Core should be identified in the original subsequence");

    // check if core will not be identified in the given refined intervals
    std::string subsequence2(sequence.begin()+start+1, sequence.begin()+end);

    struct lps lps_obj2;
    init_lps_offset(&lps_obj2, subsequence2.c_str(), subsequence2.size(), start+1);
    lps_deepen(&lps_obj2, 5);

    assert(lps_obj2.size == 0 && "Core should not be identified in the original subsequence");

    // check if core will not be identified in the given refined intervals
    std::string subsequence3(sequence.begin()+start, sequence.begin()+end-1);

    struct lps lps_obj3;
    init_lps_offset(&lps_obj3, subsequence3.c_str(), subsequence3.size(), start);
    lps_deepen(&lps_obj3, 5);

    assert(lps_obj3.size == 0 && "Core should not be identified in the original subsequence");

    free_lps(&lps_obj);
    free_lps(&lps_obj1);
    free_lps(&lps_obj2);
    free_lps(&lps_obj3);

    log("...  test_lps_consistency passed!");
}

int main() {

	log("Running test_lps...");

    test_lps_split_init();
    test_lps_file_io();
    test_lps_consistency();

	log("All tests in test_lps completed successfully!");

	return 0;
}
