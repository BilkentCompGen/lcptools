#include "core.h"
#include "lps.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

void log(const std::string &message) {
	std::cout << message << std::endl;
};

void test_lps_constructor() {

	LCP_INIT();

	std::string test_string = "GGGACCTgGTGACCCCAGCcCACGaCAGCCAAGCGCCAGCTGAGCtCAGGTGTGAGGAGATCacaGTCCT";
	// create an lps object
    struct lps lps_obj;
    init_lps(&lps_obj, test_string.c_str(), test_string.size());

	// cores to compare
	std::vector<struct core *> cores;
	cores.resize(31);

	ublock *p0 = (ublock*)malloc(sizeof(ublock));
	p0[0] = 0b100001;
	cores[0] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[0], 6, p0, 1, 0, 0);
	ublock *p1 = (ublock*)malloc(sizeof(ublock));
	p1[0] = 0b00010111;
	cores[1] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[1], 8, p1, 2, 0, 0);
	ublock *p2 = (ublock*)malloc(sizeof(ublock));
	p2[0] = 0b011110;
	cores[2] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[2], 6, p2, 3, 0, 0);
	ublock *p3 = (ublock*)malloc(sizeof(ublock));
	p3[0] = 0b11101011;
	cores[3] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[3], 8, p3, 4, 0, 0);
	ublock *p4 = (ublock*)malloc(sizeof(ublock));
	p4[0] = 0b101110;
	cores[4] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[4], 6, p4, 5, 0, 0);
	ublock *p5 = (ublock*)malloc(sizeof(ublock));
	p5[0] = 0b100001;
	cores[5] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[5], 6, p5, 6, 0, 0);
	ublock *p6 = (ublock*)malloc(sizeof(ublock));
	p6[0] = 0b000101010100;
	cores[6] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[6], 12, p6, 7, 0, 0);
	ublock *p7 = (ublock*)malloc(sizeof(ublock));
	p7[0] = 0b010010;
	cores[7] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[7], 6, p7, 8, 0, 0);
	ublock *p8 = (ublock*)malloc(sizeof(ublock));
	p8[0] = 0b1001010100;
	cores[8] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[8], 10, p8, 9, 0, 0);
	ublock *p9 = (ublock*)malloc(sizeof(ublock));
	p9[0] = 0b010001;
	cores[9] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[9], 6, p9, 10, 0, 0);
	ublock *p10 = (ublock*)malloc(sizeof(ublock));
	p10[0] = 0b100001;
	cores[10] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[10], 6, p10, 11, 0, 0);

	ublock *p11 = (ublock*)malloc(sizeof(ublock));
	p11[0] = 0b010010;
	cores[11] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[11], 6, p11, 12, 0, 0);
	ublock *p12 = (ublock*)malloc(sizeof(ublock));
	p12[0] = 0b10010100;
	cores[12] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[12], 8, p12, 13, 0, 0);
	ublock *p13 = (ublock*)malloc(sizeof(ublock));
	p13[0] = 0b01000010;
	cores[13] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[13], 8, p13, 14, 0, 0);
	ublock *p14 = (ublock*)malloc(sizeof(ublock));
	p14[0] = 0b100110;
	cores[14] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[14], 6, p14, 15, 0, 0);
	ublock *p15 = (ublock*)malloc(sizeof(ublock));
	p15[0] = 0b10010100;
	cores[15] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[15], 8, p15, 16, 0, 0);
	ublock *p16 = (ublock*)malloc(sizeof(ublock));
	p16[0] = 0b010010;
	cores[16] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[16], 6, p16, 17, 0, 0);
	ublock *p17 = (ublock*)malloc(sizeof(ublock));
	p17[0] = 0b100111;
	cores[17] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[17], 6, p17, 18, 0, 0);
	ublock *p18 = (ublock*)malloc(sizeof(ublock));
	p18[0] = 0b100010;
	cores[18] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[18], 6, p18, 19, 0, 0);
	ublock *p19 = (ublock*)malloc(sizeof(ublock));
	p19[0] = 0b100111;
	cores[19] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[19], 6, p19, 20, 0, 0);
	ublock *p20 = (ublock*)malloc(sizeof(ublock));
	p20[0] = 0b010010;
	cores[20] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[20], 6, p20, 21, 0, 0);

	ublock *p21 = (ublock*)malloc(sizeof(ublock));
	p21[0] = 0b00101011;
	cores[21] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[21], 8, p21, 22, 0, 0);
	ublock *p22 = (ublock*)malloc(sizeof(ublock));
	p22[0] = 0b111011;
	cores[22] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[22], 6, p22, 23, 0, 0);
	ublock *p23 = (ublock*)malloc(sizeof(ublock));
	p23[0] = 0b100010;
	cores[23] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[23], 6, p23, 24, 0, 0);
	ublock *p24 = (ublock*)malloc(sizeof(ublock));
	p24[0] = 0b00101000;
	cores[24] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[24], 8, p24, 25, 0, 0);
	ublock *p25 = (ublock*)malloc(sizeof(ublock));
	p25[0] = 0b100010;
	cores[25] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[25], 6, p25, 26, 0, 0);
	ublock *p26 = (ublock*)malloc(sizeof(ublock));
	p26[0] = 0b100011;
	cores[26] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[26], 6, p26, 27, 0, 0);
	ublock *p27 = (ublock*)malloc(sizeof(ublock));
	p27[0] = 0b010001;
	cores[27] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[27], 6, p27, 28, 0, 0);
	ublock *p28 = (ublock*)malloc(sizeof(ublock));
	p28[0] = 0b010010;
	cores[28] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[28], 6, p28, 29, 0, 0);
	ublock *p29 = (ublock*)malloc(sizeof(ublock));
	p29[0] = 0b101101;
	cores[29] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[29], 6, p29, 30, 0, 0);
	ublock *p30 = (ublock*)malloc(sizeof(ublock));
	p30[0] = 0b11010111;
	cores[30] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[30], 8, p30, 31, 0, 0);

	// compare the resulting cores at level 1
	assert(lps_obj.size == static_cast<int>(cores.size()) && "Core size at level 1 should match");
	for (int i = 0; i < lps_obj.size; i++) {
		assert(core_eq(&(lps_obj.cores[i]), cores[i]) && "Cores at level 1 should match");
	}

	for (size_t i = 0; i < cores.size(); i++) {
		free_core(cores[i]);
        free(cores[i]);
	}

    free_lps(&lps_obj);

	log("...  test_lps_constructor passed!");
};

void test_lps_reverse_complement() {

    LCP_INIT();

	std::string test_string = "AGGACTgtgatCTCCTCACACCTGAGCTCAGCTGGCGCTTGGCTGTCGtGggCTGGGGTCAccAGGTCCC";

	// create an lps object
    struct lps lps_obj;
    init_lps2(&lps_obj, test_string.c_str(), test_string.size());

	// cores to compare
	std::vector<struct core *> cores;
	cores.resize(31);

	ublock *p0 = (ublock*)malloc(sizeof(ublock));
	p0[0] = 0b100001;
	cores[0] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[0], 6, p0, 1, 0, 0);
	ublock *p1 = (ublock*)malloc(sizeof(ublock));
	p1[0] = 0b00010111;
	cores[1] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[1], 8, p1, 2, 0, 0);
	ublock *p2 = (ublock*)malloc(sizeof(ublock));
	p2[0] = 0b011110;
	cores[2] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[2], 6, p2, 3, 0, 0);
	ublock *p3 = (ublock*)malloc(sizeof(ublock));
	p3[0] = 0b11101011;
	cores[3] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[3], 8, p3, 4, 0, 0);
	ublock *p4 = (ublock*)malloc(sizeof(ublock));
	p4[0] = 0b101110;
	cores[4] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[4], 6, p4, 5, 0, 0);
	ublock *p5 = (ublock*)malloc(sizeof(ublock));
	p5[0] = 0b100001;
	cores[5] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[5], 6, p5, 6, 0, 0);
	ublock *p6 = (ublock*)malloc(sizeof(ublock));
	p6[0] = 0b000101010100;
	cores[6] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[6], 12, p6, 7, 0, 0);
	ublock *p7 = (ublock*)malloc(sizeof(ublock));
	p7[0] = 0b010010;
	cores[7] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[7], 6, p7, 8, 0, 0);
	ublock *p8 = (ublock*)malloc(sizeof(ublock));
	p8[0] = 0b1001010100;
	cores[8] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[8], 10, p8, 9, 0, 0);
	ublock *p9 = (ublock*)malloc(sizeof(ublock));
	p9[0] = 0b010001;
	cores[9] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[9], 6, p9, 10, 0, 0);
	ublock *p10 = (ublock*)malloc(sizeof(ublock));
	p10[0] = 0b100001;
	cores[10] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[10], 6, p10, 11, 0, 0);

	ublock *p11 = (ublock*)malloc(sizeof(ublock));
	p11[0] = 0b010010;
	cores[11] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[11], 6, p11, 12, 0, 0);
	ublock *p12 = (ublock*)malloc(sizeof(ublock));
	p12[0] = 0b10010100;
	cores[12] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[12], 8, p12, 13, 0, 0);
	ublock *p13 = (ublock*)malloc(sizeof(ublock));
	p13[0] = 0b01000010;
	cores[13] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[13], 8, p13, 14, 0, 0);
	ublock *p14 = (ublock*)malloc(sizeof(ublock));
	p14[0] = 0b100110;
	cores[14] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[14], 6, p14, 15, 0, 0);
	ublock *p15 = (ublock*)malloc(sizeof(ublock));
	p15[0] = 0b10010100;
	cores[15] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[15], 8, p15, 16, 0, 0);
	ublock *p16 = (ublock*)malloc(sizeof(ublock));
	p16[0] = 0b010010;
	cores[16] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[16], 6, p16, 17, 0, 0);
	ublock *p17 = (ublock*)malloc(sizeof(ublock));
	p17[0] = 0b100111;
	cores[17] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[17], 6, p17, 18, 0, 0);
	ublock *p18 = (ublock*)malloc(sizeof(ublock));
	p18[0] = 0b100010;
	cores[18] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[18], 6, p18, 19, 0, 0);
	ublock *p19 = (ublock*)malloc(sizeof(ublock));
	p19[0] = 0b100111;
	cores[19] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[19], 6, p19, 20, 0, 0);
	ublock *p20 = (ublock*)malloc(sizeof(ublock));
	p20[0] = 0b010010;
	cores[20] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[20], 6, p20, 21, 0, 0);

	ublock *p21 = (ublock*)malloc(sizeof(ublock));
	p21[0] = 0b00101011;
	cores[21] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[21], 8, p21, 22, 0, 0);
	ublock *p22 = (ublock*)malloc(sizeof(ublock));
	p22[0] = 0b111011;
	cores[22] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[22], 6, p22, 23, 0, 0);
	ublock *p23 = (ublock*)malloc(sizeof(ublock));
	p23[0] = 0b100010;
	cores[23] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[23], 6, p23, 24, 0, 0);
	ublock *p24 = (ublock*)malloc(sizeof(ublock));
	p24[0] = 0b00101000;
	cores[24] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[24], 8, p24, 25, 0, 0);
	ublock *p25 = (ublock*)malloc(sizeof(ublock));
	p25[0] = 0b100010;
	cores[25] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[25], 6, p25, 26, 0, 0);
	ublock *p26 = (ublock*)malloc(sizeof(ublock));
	p26[0] = 0b100011;
	cores[26] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[26], 6, p26, 27, 0, 0);
	ublock *p27 = (ublock*)malloc(sizeof(ublock));
	p27[0] = 0b010001;
	cores[27] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[27], 6, p27, 28, 0, 0);
	ublock *p28 = (ublock*)malloc(sizeof(ublock));
	p28[0] = 0b010010;
	cores[28] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[28], 6, p28, 29, 0, 0);
	ublock *p29 = (ublock*)malloc(sizeof(ublock));
	p29[0] = 0b101101;
	cores[29] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[29], 6, p29, 30, 0, 0);
	ublock *p30 = (ublock*)malloc(sizeof(ublock));
	p30[0] = 0b11010111;
	cores[30] = (struct core*)malloc(sizeof(struct core));
    init_core4(cores[30], 8, p30, 31, 0, 0);

	// compare the resulting cores at level 1
	assert(lps_obj.size == static_cast<int>(cores.size()) && "Core size at level 1 should match");
	for (int i = 0; i < lps_obj.size; i++) {
		assert(core_eq(&(lps_obj.cores[i]), cores[i]) && "Cores at level 1 should match");
	}

	for (size_t i = 0; i < cores.size(); i++) {
		free_core(cores[i]);
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
};

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

    ublock *p0 = (ublock*)malloc(sizeof(ublock));
	p0[0] = 0b110001;
	level2_cores[0] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[0], 6, p0, 1, 0, 0);

	ublock *p1 = (ublock*)malloc(sizeof(ublock));
	p1[0] = 0b010001;
	level2_cores[1] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[1], 6, p1, 1, 0, 0);

	ublock *p2 = (ublock*)malloc(sizeof(ublock));
	p2[0] = 0b00010011;
	level2_cores[2] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[2], 6, p2, 2, 0, 0);

	ublock *p3 = (ublock*)malloc(sizeof(ublock));
	p3[0] = 0b10011000;
	level2_cores[3] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[3], 8, p3, 3, 0, 0);

	ublock *p4 = (ublock*)malloc(sizeof(ublock));
	p4[0] = 0b10000010;
	level2_cores[4] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[4], 8, p4, 4, 0, 0);

	ublock *p5 = (ublock*)malloc(sizeof(ublock));
	p5[0] = 0b1011101;
	level2_cores[5] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[5], 7, p5, 4, 0, 0);

	ublock *p6 = (ublock*)malloc(sizeof(ublock));
	p6[0] = 0b1011011;
	level2_cores[6] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[6], 7, p6, 5, 0, 0);

	ublock *p7 = (ublock*)malloc(sizeof(ublock));
	p7[0] = 0b010001;
	level2_cores[7] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[7], 6, p7, 6, 0, 0);

	ublock *p8 = (ublock*)malloc(sizeof(ublock));
	p8[0] = 0b010001;
	level2_cores[8] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[8], 6, p8, 7, 0, 0);

	ublock *p9 = (ublock*)malloc(sizeof(ublock));
	p9[0] = 0b10010010;
	level2_cores[9] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[9], 8, p9, 8, 0, 0);

	ublock *p10 = (ublock*)malloc(sizeof(ublock));
	p10[0] = 0b110110;
	level2_cores[10] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[10], 6, p10, 9, 0, 0);

	ublock *p11 = (ublock*)malloc(sizeof(ublock));
	p11[0] = 0b100001;
	level2_cores[11] = (struct core*)malloc(sizeof(struct core));
    init_core4(level2_cores[11], 6, p11, 10, 0, 0);

	assert(lps_obj.size == static_cast<int>(level2_cores.size()) && "Core size at level 2 should match");
	for (int i = 0; i < lps_obj.size; i++) {
		assert(core_eq(&(lps_obj.cores[i]), level2_cores[i]) && "Cores at level 2 should match");
	}

	// test level 3 cores
	std::vector<struct core *> level3_cores;
	level3_cores.resize(4);

	ublock *p12 = (ublock*)malloc(sizeof(ublock));
	p12[0] = 0b110011;
	level3_cores[0] = (struct core*)malloc(sizeof(struct core));
    init_core4(level3_cores[0], 6, p12, 1, 0, 0);

	ublock *p13 = (ublock*)malloc(sizeof(ublock));
	p13[0] = 0b110111;
	level3_cores[1] = (struct core*)malloc(sizeof(struct core));
    init_core4(level3_cores[1], 6, p13, 1, 0, 0);

	ublock *p14 = (ublock*)malloc(sizeof(ublock));
	p14[0] = 0b11101100;
	level3_cores[2] = (struct core*)malloc(sizeof(struct core));
    init_core4(level3_cores[2], 8, p14, 1, 0, 0);

	ublock *p15 = (ublock*)malloc(sizeof(ublock));
	p15[0] = 0b110000101;
	level3_cores[3] = (struct core*)malloc(sizeof(struct core));
    init_core4(level3_cores[3], 9, p15, 1, 0, 0);

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
		free_core(level2_cores[i]);
        free(level2_cores[i]);
	}

    for (size_t i = 0; i < level3_cores.size(); i++) {
		free_core(level3_cores[i]);
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
