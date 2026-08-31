#include "core.h"
#include <cassert>
#include <string>
#include <iostream>


inline void log(const std::string &message) {
	std::cout << message << std::endl;
}

void test_core_constructors() {

	struct core core1;
    init_core4(&core1, 4, 0b1111, 2, 0, 10); // 1111 in binary

	assert(core1.bit_size == 4 && "Label length should be 4");
	assert(core1.bit_rep == 0b1111 && "Label should be 0b1111");
#if LCP_LABEL_BITS != 0
	assert(core1.label == 2 && "Core label should be 2");
#endif
#if LCP_POS_BITS != 0
	assert(core1.start == 0 && "Start should be 0");
	assert(core1.end == 10 && "End should be 10");
#endif

	log("...  test_core_constructors passed!");
}

#if !LCP_ALPHABET_PROTEIN

void test_core_compress() {

	struct core core1;
    init_core4(&core1, 3, (1ull << 63) | 0b01110001, 10, 0, 3); // 101 in binary
	struct core core2;
    init_core4(&core2, 3, (1ull << 63) | 0b10100001, 11, 0, 3); // 111 in binary

    core_compress_level1(&core2, &core1);

	// expected result after compressing 101 and 111 is 10 (binary) => 2 in decimal
	assert(core1.bit_rep == 0b1001 && "Compressed core's label should be 0b1001");
	assert(core1.bit_size == 4 && "Compressed core's label length should be 4");
#if LCP_LABEL_BITS != 0
	assert(core1.label == 10 && "Core's label should be 10");
#endif

	struct core core3;
    init_core4(&core3, 3, 0b101, 10, 0, 3); // 101 in binary
	struct core core4;
    init_core4(&core4, 3, 0b111, 11, 0, 3); // 111 in binary

    core_compress_upper(&core4, &core3);

	// expected result after compressing 101 and 111 is 10 (binary) => 2 in decimal
	assert(core3.bit_rep == 0b10 && "Compressed core's label should be 0b10");
	assert(core3.bit_size == 2 && "Compressed core's label length should be 2");

	struct core core5;
    init_core4(&core5, 3, (1ull << 63) | 0b11100001, 10, 0, 3); // 101 in binary
	struct core core6;
    init_core4(&core6, 3, (1ull << 63) | 0b10110001, 11, 0, 3); // 111 in binary

    core_compress_level1(&core6, &core5);

	// expected result after compressing 101 and 111 is 10 (binary) => 2 in decimal
	assert(core5.bit_rep == 0b1100 && "Compressed core's label should be 0b1110");
	assert(core5.bit_size == 4 && "Compressed core's label length should be 4");


	log("...  test_core_compress passed!");
}

#endif /* !LCP_ALPHABET_PROTEIN */

void test_core_operator_overloads() {

	struct core core1;
    init_core4(&core1, 4, 0b1010, 0, 0, 0); // 1010 in binary
	struct core core2;
    init_core4(&core2, 4, 0b1010, 1, 1, 0); // 1010 in binary
	struct core core3;
    init_core4(&core3, 3, 0b101, 2, 2, 0); // 101 in binary

	assert(core_eq(&core1, &core2) && "core1 should be equal to core2");
	assert(core_neq(&core1, &core3) && "core1 should not be equal to core3");
	assert(core_lt(&core3, &core1) && "core3 should be less than core1");
	assert(core_gt(&core1, &core3) && "core1 should be greater than core3");
	assert(core_geq(&core1, &core2) && "core1 should be greater than or equal to core2");
	assert(core_leq(&core3, &core1) && "core3 should be less than or equal to core1");

	log("...  test_core_operator_overloads passed!");
}

int main() {
	log("Running test_core...");

	test_core_constructors();
#if !LCP_ALPHABET_PROTEIN
	test_core_compress();
#endif
	test_core_operator_overloads();

	log("All tests in test_core completed successfully!");

	return 0;
}