#include "core.h"
#include "lps.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>

#define CORE_COUNT(arr) (static_cast<int>(sizeof(arr) / sizeof((arr)[0])))


inline void log(const std::string &message) {
	std::cout << message << std::endl;
}

// reads the first sequence from a FASTA file
inline std::string read_first_sequence(const char *path) {
	std::ifstream genome(path);
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

	return sequence;
}

// expected core, only the bit representation matters for core_eq
struct expected_core {
	ubit_size bit_size;
	ublock bits;
};

// level 1 cores of "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT"
// (same result as the reverse complement of "AGGACTGTGATCTCCTCACACCTGAGCTCAGCTGGCGCTTGGCTGTCGTGGGCTGGGGTCACCAGGTCCC")
static const expected_core LEVEL1_CORES[] = {
	{6, 0b100001},
	{8, 0b00010111},
	{6, 0b011110},
	{8, 0b11101011},
	{6, 0b101110},
	{6, 0b100001},
	{12, 0b000101010100},
	{6, 0b010010},
	{10, 0b1001010100},
	{6, 0b010001},
	{6, 0b100001},
	{6, 0b010010},
	{8, 0b10010100},
	{8, 0b01000010},
	{6, 0b100110},
	{8, 0b10010100},
	{6, 0b010010},
	{6, 0b100111},
	{6, 0b100010},
	{6, 0b100111},
	{6, 0b010010},
	{8, 0b00101011},
	{6, 0b111011},
	{6, 0b100010},
	{8, 0b00101000},
	{6, 0b100010},
	{6, 0b100011},
	{6, 0b010001},
	{6, 0b010010},
	{6, 0b101101},
	{8, 0b11010111},
};

static const expected_core LEVEL2_CORES[] = {
	{6, 0b110001},
	{6, 0b010001},
	{6, 0b00010011},
	{8, 0b10011000},
	{8, 0b10000010},
	{7, 0b1011101},
	{7, 0b1011011},
	{6, 0b010001},
	{6, 0b010001},
	{8, 0b10010010},
	{6, 0b110110},
	{6, 0b100001},
};

static const expected_core LEVEL3_CORES[] = {
	{6, 0b110011},
	{6, 0b110111},
	{8, 0b11101100},
	{9, 0b110000101},
};

// compares the cores of an lps object against an expected table
static void assert_cores_match(const struct lps *lps_ptr, const expected_core *expected, int count, const char *what) {
	assert(lps_ptr->size == count && "Core size should match");
	for (int i = 0; i < count; i++) {
		if (lps_ptr->cores[i].bit_size != expected[i].bit_size || lps_ptr->cores[i].bit_rep[0] != expected[i].bits) {
			std::cerr << what << ": core mismatch at index " << i << std::endl;
			assert(false && "Cores should match");
		}
	}
}

void test_lps_constructor() {

	LCP_INIT();

	std::string test_string = "GGGACCTgGTGACCCCAGCcCACGaCAGCCAAGCGCCAGCTGAGCtCAGGTGTGAGGAGATCacaGTCCT";
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	assert_cores_match(&lps_obj, LEVEL1_CORES, CORE_COUNT(LEVEL1_CORES), "test_lps_constructor");

	free_lps(&lps_obj);

	log("...  test_lps_constructor passed!");
};

void test_lps_reverse_complement() {

	LCP_INIT();

	std::string test_string = "AGGACTgtgatCTCCTCACACCTGAGCTCAGCTGGCGCTTGGCTGTCGtGggCTGGGGTCAccAGGTCCC";
	struct lps lps_obj;
	init_lps2(&lps_obj, test_string.c_str(), test_string.size());

	assert_cores_match(&lps_obj, LEVEL1_CORES, CORE_COUNT(LEVEL1_CORES), "test_lps_reverse_complement");

	free_lps(&lps_obj);

	log("...  test_lps_reverse_complement passed!");
}

void test_lps_split_init() {

	LCP_INIT();

	std::string sequence = read_first_sequence("data/test.fasta");

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
};

void test_lps_deepen() {

	LCP_INIT();

	std::string test_string = "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT";
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	// deepen to level 2
	int success = lps_deepen(&lps_obj, 2);
	assert(success && "Deepening to level 2 should be successful");

	assert_cores_match(&lps_obj, LEVEL2_CORES, CORE_COUNT(LEVEL2_CORES), "test_lps_deepen (level 2)");

	// deepen to level 3
	success = lps_deepen(&lps_obj, 3);
	assert(success && "Deepening to level 3 should be successful");

	assert_cores_match(&lps_obj, LEVEL3_CORES, CORE_COUNT(LEVEL3_CORES), "test_lps_deepen (level 3)");

	// attempt to deepen to a lower level (should not do anything)
	success = lps_deepen(&lps_obj, 3);
	assert(success == 0 && "Deepening to a lower level should be unsuccessful");

	free_lps(&lps_obj);

	log("...  test_lps_deepen passed!");
}

void test_lps_consistency() {

	LCP_INIT();

	std::string sequence = read_first_sequence("data/test.fasta");

	struct lps lps_obj;
	init_lps(&lps_obj, sequence.c_str(), sequence.size());
	lps_deepen(&lps_obj, 5);

	int start = lps_obj.cores[5000].start;
	int end = lps_obj.cores[5000].end;

	// check if core is identified in the given intervals
	std::string subsequence1(sequence.begin() + start, sequence.begin() + end);

	struct lps lps_obj1;
	init_lps_offset(&lps_obj1, subsequence1.c_str(), subsequence1.size(), start);
	lps_deepen(&lps_obj1, 5);

	assert(core_eq(&(lps_obj.cores[5000]), &(lps_obj1.cores[0])) && "Core should be identified in the original subsequence");

	// check if core will not be identified in the given refined intervals
	std::string subsequence2(sequence.begin() + start + 1, sequence.begin() + end);

	struct lps lps_obj2;
	init_lps_offset(&lps_obj2, subsequence2.c_str(), subsequence2.size(), start + 1);
	lps_deepen(&lps_obj2, 5);

	assert(lps_obj2.size == 0 && "Core should not be identified in the original subsequence");

	// check if core will not be identified in the given refined intervals
	std::string subsequence3(sequence.begin() + start, sequence.begin() + end - 1);

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
