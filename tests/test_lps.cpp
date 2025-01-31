#include "core.h"
#include "lps.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

void log(const std::string &message) {
	std::cout << message << std::endl;
}

void test_lps_constructor() {

	LCP_INIT();

	std::string test_string = "GGGACCTgGTGACCCCAGCcCACGaCAGCCAAGCGCCAGCTGAGCtCAGGTGTGAGGAGATCacaGTCCT";
	// create an lps object
    struct lps lps_obj;
    init_lps(&lps_obj, test_string.c_str(), test_string.size());

	// cores to compare
	std::vector<struct core *> cores;
	cores.resize(31);

	cores[0] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[0], 6, 0x8000000000000000 | 0b01100001, 1, 0, 0);

	cores[1] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[1], 8, 0x8000000000000000 | 0b10000111, 2, 0, 0);

	cores[2] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[2], 6, 0x8000000000000000 | 0b01011110, 3, 0, 0);

	cores[3] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[3], 8, 0x8000000000000000 | 0b10111011, 4, 0, 0);

	cores[4] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[4], 6, 0x8000000000000000 | 0b01101110, 5, 0, 0);

	cores[5] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[5], 6, 0x8000000000000000 | 0b01100001, 6, 0, 0);

	cores[6] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[6], 12, 0x8000000000000000 | 0b100000100, 7, 0, 0);

	cores[7] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[7], 6, 0x8000000000000000 | 0b01010010, 8, 0, 0);

	cores[8] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[8], 10, 0x8000000000000000 | 0b11100100, 9, 0, 0);

	cores[9] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[9], 6, 0x8000000000000000 | 0b01010001, 10, 0, 0);

	cores[10] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[10], 6, 0x8000000000000000 | 0b01100001, 11, 0, 0);

	cores[11] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[11], 6, 0x8000000000000000 | 0b01010010, 12, 0, 0);

	cores[12] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[12], 8, 0x8000000000000000 | 0b10100100, 13, 0, 0);

	cores[13] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[13], 8, 0x8000000000000000 | 0b10010010, 14, 0, 0);

	cores[14] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[14], 6, 0x8000000000000000 | 0b01100110, 15, 0, 0);

	cores[15] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[15], 8, 0x8000000000000000 | 0b10100100, 16, 0, 0);

	cores[16] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[16], 6, 0x8000000000000000 | 0b01010010, 17, 0, 0);

	cores[17] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[17], 6, 0x8000000000000000 | 0b01100111, 18, 0, 0);

	cores[18] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[18], 6, 0x8000000000000000 | 0b01100010, 19, 0, 0);

	cores[19] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[19], 6, 0x8000000000000000 | 0b01100111, 20, 0, 0);

	cores[20] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[20], 6, 0x8000000000000000 | 0b01010010, 21, 0, 0);

	cores[21] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[21], 8, 0x8000000000000000 | 0b10001011, 22, 0, 0);

	cores[22] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[22], 6, 0x8000000000000000 | 0b01111011, 23, 0, 0);

	cores[23] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[23], 6, 0x8000000000000000 | 0b01100010, 24, 0, 0);

	cores[24] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[24], 8, 0x8000000000000000 | 0b10001000, 25, 0, 0);

	cores[25] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[25], 6, 0x8000000000000000 | 0b01100010, 26, 0, 0);

	cores[26] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[26], 6, 0x8000000000000000 | 0b01100011, 27, 0, 0);

	cores[27] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[27], 6, 0x8000000000000000 | 0b01010001, 28, 0, 0);

	cores[28] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[28], 6, 0x8000000000000000 | 0b01010010, 29, 0, 0);

	cores[29] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[29], 6, 0x8000000000000000 | 0b01101101, 30, 0, 0);

	cores[30] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[30], 8, 0x8000000000000000 | 0b10110111, 31, 0, 0);

	// compare the resulting cores at level 1
	assert(lps_obj.size == static_cast<int>(cores.size()) && "Core size at level 1 should match");
	for (int i = 0; i < lps_obj.size; i++) {
		assert(core_eq(&(lps_obj.cores[i]), cores[i]) && "Cores at level 1 should match");
	}

	for (size_t i = 0; i < cores.size(); i++) {
        free(cores[i]);
	}

    free_lps(&lps_obj);

	log("...  test_lps_constructor passed!");
}

void test_lps_reverse_complement() {

    LCP_INIT();

	std::string test_string = "AGGACTgtgatCTCCTCACACCTGAGCTCAGCTGGCGCTTGGCTGTCGtGggCTGGGGTCAccAGGTCCC";

	// create an lps object
    struct lps lps_obj;
    init_lps2(&lps_obj, test_string.c_str(), test_string.size());

	// cores to compare
	std::vector<struct core *> cores;
	cores.resize(31);

	cores[0] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[0], 6, 0x8000000000000000 | 0b01100001, 1, 0, 0);

	cores[1] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[1], 8, 0x8000000000000000 | 0b10000111, 2, 0, 0);

	cores[2] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[2], 6, 0x8000000000000000 | 0b01011110, 3, 0, 0);

	cores[3] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[3], 8, 0x8000000000000000 | 0b10111011, 4, 0, 0);

	cores[4] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[4], 6, 0x8000000000000000 | 0b01101110, 5, 0, 0);

	cores[5] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[5], 6, 0x8000000000000000 | 0b01100001, 6, 0, 0);

	cores[6] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[6], 12, 0x8000000000000000 | 0b100000100, 7, 0, 0);

	cores[7] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[7], 6, 0x8000000000000000 | 0b01010010, 8, 0, 0);

	cores[8] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[8], 10, 0x8000000000000000 | 0b11100100, 9, 0, 0);

	cores[9] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[9], 6, 0x8000000000000000 | 0b01010001, 10, 0, 0);

	cores[10] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[10], 6, 0x8000000000000000 | 0b01100001, 11, 0, 0);

	cores[11] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[11], 6, 0x8000000000000000 | 0b01010010, 12, 0, 0);

	cores[12] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[12], 8, 0x8000000000000000 | 0b10100100, 13, 0, 0);

	cores[13] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[13], 8, 0x8000000000000000 | 0b10010010, 14, 0, 0);

	cores[14] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[14], 6, 0x8000000000000000 | 0b01100110, 15, 0, 0);

	cores[15] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[15], 8, 0x8000000000000000 | 0b10100100, 16, 0, 0);

	cores[16] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[16], 6, 0x8000000000000000 | 0b01010010, 17, 0, 0);

	cores[17] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[17], 6, 0x8000000000000000 | 0b01100111, 18, 0, 0);

	cores[18] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[18], 6, 0x8000000000000000 | 0b01100010, 19, 0, 0);

	cores[19] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[19], 6, 0x8000000000000000 | 0b01100111, 20, 0, 0);

	cores[20] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[20], 6, 0x8000000000000000 | 0b01010010, 21, 0, 0);

	cores[21] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[21], 8, 0x8000000000000000 | 0b10001011, 22, 0, 0);

	cores[22] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[22], 6, 0x8000000000000000 | 0b01111011, 23, 0, 0);

	cores[23] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[23], 6, 0x8000000000000000 | 0b01100010, 24, 0, 0);

	cores[24] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[24], 8, 0x8000000000000000 | 0b10001000, 25, 0, 0);

	cores[25] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[25], 6, 0x8000000000000000 | 0b01100010, 26, 0, 0);

	cores[26] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[26], 6, 0x8000000000000000 | 0b01100011, 27, 0, 0);

	cores[27] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[27], 6, 0x8000000000000000 | 0b01010001, 28, 0, 0);

	cores[28] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[28], 6, 0x8000000000000000 | 0b01010010, 29, 0, 0);

	cores[29] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[29], 6, 0x8000000000000000 | 0b01101101, 30, 0, 0);

	cores[30] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[30], 8, 0x8000000000000000 | 0b10110111, 31, 0, 0);

	// compare the resulting cores at level 1
	assert(lps_obj.size == static_cast<int>(cores.size()) && "Core size at level 1 should match");
	for (int i = 0; i < lps_obj.size; i++) {
		assert(core_eq(&(lps_obj.cores[i]), cores[i]) && "Cores at level 1 should match");
	}

	for (size_t i = 0; i < cores.size(); i++) {
        free(cores[i]);
	}

    free_lps(&lps_obj);

	log("...  test_lps_reverse_complement passed!");
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

void test_lps_deepen() {

	LCP_INIT();

	std::string test_string = "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT";
    struct lps lps_obj;
    init_lps(&lps_obj, test_string.c_str(), test_string.size());

	// deepen to level 2
	int success = lps_deepen(&lps_obj, 2);
	assert(success && "Deepening to level 2 should be successful");

	std::vector<struct core *> level2_cores;
	level2_cores.resize(12);

	level2_cores[0] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[0], 6, 0b110001, 1, 0, 0);

	level2_cores[1] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[1], 6, 0b010001, 1, 0, 0);

	level2_cores[2] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[2], 6, 0b010011, 2, 0, 0);

	level2_cores[3] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[3], 8, 0b10011000, 3, 0, 0);

	level2_cores[4] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[4], 8, 0b10000010, 4, 0, 0);

	level2_cores[5] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[5], 7, 0b1011101, 4, 0, 0);

	level2_cores[6] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[6], 7, 0b1011011, 5, 0, 0);

	level2_cores[7] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[7], 6, 0b010001, 6, 0, 0);

	level2_cores[8] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[8], 6, 0b010001, 7, 0, 0);

	level2_cores[9] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[9], 8, 0b10010010, 8, 0, 0);

	level2_cores[10] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[10], 6, 0b110110, 9, 0, 0);

	level2_cores[11] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[11], 6, 0b100001, 10, 0, 0);

	assert(lps_obj.size == static_cast<int>(level2_cores.size()) && "Core size at level 2 should match");
	for (int i = 0; i < lps_obj.size; i++) {
		assert(core_eq(&(lps_obj.cores[i]), level2_cores[i]) && "Cores at level 2 should match");
	}

	// test level 3 cores
	std::vector<struct core *> level3_cores;
	level3_cores.resize(4);

	level3_cores[0] = (struct core*)malloc(sizeof(struct core));
    init_core4(level3_cores[0], 6, 0b110011, 1, 0, 0);

	level3_cores[1] = (struct core*)malloc(sizeof(struct core));
    init_core4(level3_cores[1], 6, 0b110111, 1, 0, 0);

	level3_cores[2] = (struct core*)malloc(sizeof(struct core));
    init_core4(level3_cores[2], 8, 0b11101100, 1, 0, 0);

	level3_cores[3] = (struct core*)malloc(sizeof(struct core));
    init_core4(level3_cores[3], 9, 0b110000101, 1, 0, 0);

	// deepen to level 3
	success = lps_deepen(&lps_obj, 3);
	assert(success && "Deepening to level 3 should be successful");

	assert(lps_obj.size == static_cast<int>(level3_cores.size()) && "Core size at level 3 should match");
	for (int i = 0; i < lps_obj.size; i++) {
		assert(core_eq(&(lps_obj.cores[i]), level3_cores[i]) && "Cores at level 3 should match");
	}

	// attempt to deepen to a lower level (should not do anything)
	success = lps_deepen(&lps_obj, 3);
	assert(success == 0 && "Deepening to a lower level should be unsuccessful");

    free_lps(&lps_obj);

    for (size_t i = 0; i < level2_cores.size(); i++) {
        free(level2_cores[i]);
	}

    for (size_t i = 0; i < level3_cores.size(); i++) {
        free(level3_cores[i]);
	}

	log("...  test_lps_deepen passed!");
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

    int start = lps_obj.cores[5000].start;
    int end = lps_obj.cores[5000].end;

    // check if core is identified in the given intervals
    std::string subsequence1(sequence.begin()+start, sequence.begin()+end);

    struct lps lps_obj1;
    init_lps_offset(&lps_obj1, subsequence1.c_str(), subsequence1.size(), start);
    lps_deepen(&lps_obj1, 5);

    assert(core_eq(&(lps_obj.cores[5000]), &(lps_obj1.cores[0])) && "Core should be identified in the original subsequence");

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

	test_lps_constructor();
    test_lps_reverse_complement();
    test_lps_split_init();
    test_lps_file_io();
	test_lps_deepen();
    test_lps_consistency();

	log("All tests in test_lps completed successfully!");

	return 0;
}
