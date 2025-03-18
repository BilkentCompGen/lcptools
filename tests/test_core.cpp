#include "core.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

void log(const std::string &message) {
	std::cout << message << std::endl;
}

void test_core_constructors() {

	struct core core1;
    init_core4(&core1, 2, 0, 10); // 1111 in binary

	assert(core1.label == 2 && "Core label should be 1");
	assert(core1.start == 0 && "Start should be 0");
	assert(core1.end == 10 && "End should be 10");

	log("...  test_core_constructors passed!");
}

void test_core_operator_overloads() {

	struct core core1;
    init_core4(&core1, 0b1010, 0, 0); // 1010 in binary
	struct core core2;
    init_core4(&core2, 0b1010, 1, 1); // 1010 in binary
	struct core core3;
    init_core4(&core3, 0b101, 2, 2); // 101 in binary

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
    test_core_operator_overloads();

	log("All tests in test_core completed successfully!");

	return 0;
}