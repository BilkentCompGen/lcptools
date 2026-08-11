#include "encoding.h"
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>


inline void log(const std::string &message) {
	std::cout << message << std::endl;
}

// expected encoding and reverse complement encoding for a symbol
struct encoding_case {
	char symbol;
	int encoding;
	int rc_encoding;
};

static void assert_encodings(const encoding_case *cases, int count) {
	for (int i = 0; i < count; i++) {
		assert(alphabet[(int)cases[i].symbol] == cases[i].encoding && "Symbol encoding should match");
		assert(rc_alphabet[(int)cases[i].symbol] == cases[i].rc_encoding && "Reverse complement encoding should match");
	}
}

void test_encoding_initialization_default() {

	LCP_INIT();

	static const encoding_case cases[] = {
		{'A', 0, 3}, {'C', 1, 2}, {'G', 2, 1}, {'T', 3, 0},
		{'a', 0, 3}, {'c', 1, 2}, {'g', 2, 1}, {'t', 3, 0},
	};
	assert_encodings(cases, sizeof(cases) / sizeof(cases[0]));

	// check dictionary bit size
	assert(alphabet_bit_size == 2 && "Alphabet bit size should be 2");

	log("...  test_encoding_initialization_default passed!");
};

void test_encoding_initialization_from_file() {

	static const encoding_case cases[] = {
		{'A', 5, 2}, {'C', 3, 3}, {'G', 7, 0}, {'T', 8, 1},
	};
	const int count = sizeof(cases) / sizeof(cases[0]);

	// create a temporary encoding file
	std::ofstream encoding_file("encoding_test.txt");
	for (int i = 0; i < count; i++) {
		encoding_file << cases[i].symbol << " " << cases[i].encoding << " " << cases[i].rc_encoding << "\n";
	}
	encoding_file.close();

	LCP_INIT_FILE("encoding_test.txt", 0);

	assert_encodings(cases, count);

	// check dictionary bit size
	assert(alphabet_bit_size == 4 && "Alphabet bit size should be 4");

	// clean up the temporary file
	std::remove("encoding_test.txt");

	log("...  test_encoding_initialization_from_file passed!");
};

int main() {
	log("Running test_encoding...");

	test_encoding_initialization_default();
	test_encoding_initialization_from_file();

	log("All tests in test_encoding completed successfully!");

	return 0;
};
