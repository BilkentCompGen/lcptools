#include "encoding.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>

void log(const std::string &message) {
	std::cout << message << std::endl;
}

void test_encoding_initialization_default() {

	LCP_INIT();

	// check default alphabet
	assert(alphabet['A'] == 0 && "A should be encoded as 0");
	assert(alphabet['C'] == 1 && "C should be encoded as 1");
	assert(alphabet['G'] == 2 && "G should be encoded as 2");
	assert(alphabet['T'] == 3 && "T should be encoded as 3");
	assert(alphabet['a'] == 0 && "a should be encoded as 0");
	assert(alphabet['c'] == 1 && "c should be encoded as 1");
	assert(alphabet['g'] == 2 && "g should be encoded as 2");
	assert(alphabet['t'] == 3 && "t should be encoded as 3");

	// check reverse complement alphabet
	assert(rc_alphabet['A'] == 3 && "Reverse complement of A should be 3");
	assert(rc_alphabet['C'] == 2 && "Reverse complement of C should be 2");
	assert(rc_alphabet['G'] == 1 && "Reverse complement of G should be 1");
	assert(rc_alphabet['T'] == 0 && "Reverse complement of T should be 0");
	assert(rc_alphabet['a'] == 3 && "Reverse complement of a should be 3");
	assert(rc_alphabet['c'] == 2 && "Reverse complement of c should be 2");
	assert(rc_alphabet['g'] == 1 && "Reverse complement of g should be 1");
	assert(rc_alphabet['t'] == 0 && "Reverse complement of t should be 0");

	log("...  test_encoding_initialization_default passed!");
}

void test_encoding_initialization_from_file() {

	// create a temporary encoding file
	std::ofstream encoding_file("encoding_test.txt");
	encoding_file << "A 2 2\n";
	encoding_file << "C 3 3\n";
	encoding_file << "G 1 0\n";
	encoding_file << "T 0 1\n";
	encoding_file.close();

    LCP_INIT_FILE("encoding_test.txt", 0);

	// check alphabet
	assert(alphabet['A'] == 2 && "A should be encoded as 5");
	assert(alphabet['C'] == 3 && "C should be encoded as 3");
	assert(alphabet['G'] == 1 && "G should be encoded as 7");
	assert(alphabet['T'] == 0 && "T should be encoded as 8");

	// check reverse complement alphabet
	assert(rc_alphabet['A'] == 2 && "Reverse complement of A should be 2");
	assert(rc_alphabet['C'] == 3 && "Reverse complement of C should be 3");
	assert(rc_alphabet['G'] == 0 && "Reverse complement of G should be 0");
	assert(rc_alphabet['T'] == 1 && "Reverse complement of T should be 1");

	// clean up the temporary file
	std::remove("encoding_test.txt");

	log("...  test_encoding_initialization_from_file passed!");
}

int main() {
	log("Running test_encoding...");

	test_encoding_initialization_default();
	test_encoding_initialization_from_file();

	log("All tests in test_encoding completed successfully!");

	return 0;
}
