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

// The long sequence the FASTA-based tests run on. A DNA build reads the test
// genome; a protein build cannot, since a nucleotide file exercises four of the
// 25 residues and nothing else, so it gets a generated peptide over the full
// alphabet instead.
inline std::string test_sequence() {
#if LCP_ALPHABET_PROTEIN
	static const char SYMS[] = "ABCDEFGHIKLMNOPQRSTUVWXYZ";
	std::string peptide;
	unsigned x = 12345u;
	for (int i = 0; i < 200000; i++) {
		x = x * 1103515245u + 12345u;
		peptide += SYMS[(x >> 16) % 25];
	}
	return peptide;
#else
	return read_first_sequence("data/test.fasta");
#endif
}

// Expected core. core_eq compares bit_rep alone, so the tables below pin every
// other field too: bit_size, label, and the [start, end) interval.
struct expected_core {
	ubit_size bit_size;
	uint64_t  bits;
	lcp_label label;
	lcp_pos   start;
	lcp_pos   end;
};

// compares the cores of an lps object against an expected table
static void assert_cores_match(const struct lps *lps_ptr, const expected_core *expected, int count, const char *what) {
	assert(lps_ptr->size == count && "Core size should match");
	for (int i = 0; i < count; i++) {
		const struct core *actual = &(lps_ptr->cores[i]);
		bool ok = actual->bit_size == expected[i].bit_size && actual->bit_rep == expected[i].bits;
#if LCP_LABEL_BITS != 0
		ok = ok && (uint32_t)actual->label == expected[i].label;
#endif
#if LCP_POS_BITS != 0
		ok = ok && actual->start == expected[i].start && actual->end == expected[i].end;
#endif
		if (!ok) {
			std::cerr << what << ": core mismatch at index " << i << "\n"
			          << "    expected bit_size " << expected[i].bit_size
			          << " bits 0x" << std::hex << expected[i].bits << std::dec
			          << " label " << expected[i].label
			          << " [" << expected[i].start << ", " << expected[i].end << ")\n"
			          << "    actual   bit_size " << actual->bit_size
			          << " bits 0x" << std::hex << actual->bit_rep << std::dec
#if LCP_LABEL_BITS != 0
			          << " label " << actual->label
#endif
#if LCP_POS_BITS != 0
			          << " [" << actual->start << ", " << actual->end << ")"
#endif
			          << std::endl;
			assert(false && "Cores should match");
		}
	}
}

#if LCP_ALPHABET_PROTEIN

// Golden tables for ALPHABET=protein. Residues are 8 bits, so a level 1 core
// packs first/middle/last into bits 0..23 and the run length into bits 24..62,
// below the sentinel at bit 63, and bit_size counts 8 bits per residue.
//
// LEVEL 1 is derived from the LCP rules rather than copied from this library.
// Scanning the peptide gives 25 LMIN, 4 LMAX, 2 RINT and 1 SSEQ segment, three
// of them longer than three residues; re-running that derivation reproduces all
// 32 cores below including their intervals, which is what validates the table.
// The peptide was chosen because it exercises all four segment kinds -- the DNA
// test string produces no SSEQ at all.
//
// LEVELS 2 and 3 are regression snapshots, as with the other alphabets: they
// fall out of the DCT compression and parse3, which would have to be modelled
// bit for bit to predict by hand.

// residues 1..77 of E. coli aspartyl-tRNA synthetase
#define PROTEIN_TEST_STRING "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"

#if LCP_FIXED_CORE
static const expected_core PROTEIN_LEVEL1_CORES[] = {
	{24, LCP_LEVEL1_TAG | 0x10b0912,   17500434,   0,   3},
	{24, LCP_LEVEL1_TAG | 0x1120017,   17956887,   2,   5},
	{24, LCP_LEVEL1_TAG | 0x1080009,   17301513,   5,   8},
	{24, LCP_LEVEL1_TAG | 0x10f100f,   17764367,   8,  11},
	{24, LCP_LEVEL1_TAG | 0x10f0811,   17762321,  10,  13},
	{24, LCP_LEVEL1_TAG | 0x1110514,   17892628,  12,  15},
	{24, LCP_LEVEL1_TAG | 0x1140911,   18090257,  14,  17},
	{24, LCP_LEVEL1_TAG | 0x1070511,   17237265,  17,  20},
	{24, LCP_LEVEL1_TAG | 0x111100f,   17895439,  19,  22},
	{24, LCP_LEVEL1_TAG | 0x1100f0a,   17829642,  20,  23},
	{24, LCP_LEVEL1_TAG | 0x10a0404,   17433604,  22,  25},
	{24, LCP_LEVEL1_TAG | 0x1040410,   17040400,  23,  26},
	{24, LCP_LEVEL1_TAG | 0x104100a,   17043466,  24,  27},
	{24, LCP_LEVEL1_TAG | 0x10a060a,   17434122,  26,  29},
	{24, LCP_LEVEL1_TAG | 0x1080414,   17302548,  29,  32},
	{24, LCP_LEVEL1_TAG | 0x10f000e,   17760270,  32,  35},
	{24, LCP_LEVEL1_TAG | 0x10e080a,   17696778,  34,  37},
	{24, LCP_LEVEL1_TAG | 0x1111014,   17895444,  37,  40},
	{24, LCP_LEVEL1_TAG | 0x1060306,   17171206,  40,  43},
	{24, LCP_LEVEL1_TAG | 0x106120f,   17175055,  42,  45},
	{24, LCP_LEVEL1_TAG | 0x10f030c,   17761036,  44,  47},
	{24, LCP_LEVEL1_TAG | 0x10c0a11,   17566225,  46,  49},
	{24, LCP_LEVEL1_TAG | 0x1060004,   17170436,  49,  52},
	{24, LCP_LEVEL1_TAG | 0x1090014,   17367060,  52,  55},
	{24, LCP_LEVEL1_TAG | 0x1140f14,   18091796,  54,  57},
	{24, LCP_LEVEL1_TAG | 0x1140914,   18090260,  56,  59},
	{24, LCP_LEVEL1_TAG | 0x109000a,   17367050,  59,  62},
	{24, LCP_LEVEL1_TAG | 0x10a0e03,   17436163,  61,  64},
	{24, LCP_LEVEL1_TAG | 0x103000f,   16973839,  63,  66},
	{24, LCP_LEVEL1_TAG | 0x1050414,   17105940,  66,  69},
	{24, LCP_LEVEL1_TAG | 0x1041414,   17044500,  67,  70},
	{24, LCP_LEVEL1_TAG | 0x1141407,   18093063,  68,  71},
	{24, LCP_LEVEL1_TAG | 0x1140711,   18089745,  69,  72},
	{24, LCP_LEVEL1_TAG | 0x10a0009,   17432585,  72,  75},
	{24, LCP_LEVEL1_TAG | 0x1150910,   18155792,  75,  78},
};
#else
static const expected_core PROTEIN_LEVEL1_CORES[] = {
	{24, LCP_LEVEL1_TAG | 0x10b0912,   17500434,   0,   3},
	{24, LCP_LEVEL1_TAG | 0x1120017,   17956887,   2,   5},
	{24, LCP_LEVEL1_TAG | 0x1080009,   17301513,   5,   8},
	{24, LCP_LEVEL1_TAG | 0x10f100f,   17764367,   8,  11},
	{24, LCP_LEVEL1_TAG | 0x10f0811,   17762321,  10,  13},
	{24, LCP_LEVEL1_TAG | 0x1110514,   17892628,  12,  15},
	{24, LCP_LEVEL1_TAG | 0x1140911,   18090257,  14,  17},
	{24, LCP_LEVEL1_TAG | 0x1070511,   17237265,  17,  20},
	{32, LCP_LEVEL1_TAG | 0x2110f0a,   34672394,  19,  23},
	{32, LCP_LEVEL1_TAG | 0x20a0410,   34210832,  22,  26},
	{24, LCP_LEVEL1_TAG | 0x104100a,   17043466,  24,  27},
	{24, LCP_LEVEL1_TAG | 0x10a060a,   17434122,  26,  29},
	{24, LCP_LEVEL1_TAG | 0x1080414,   17302548,  29,  32},
	{24, LCP_LEVEL1_TAG | 0x10f000e,   17760270,  32,  35},
	{24, LCP_LEVEL1_TAG | 0x10e080a,   17696778,  34,  37},
	{24, LCP_LEVEL1_TAG | 0x1111014,   17895444,  37,  40},
	{24, LCP_LEVEL1_TAG | 0x1060306,   17171206,  40,  43},
	{24, LCP_LEVEL1_TAG | 0x106120f,   17175055,  42,  45},
	{24, LCP_LEVEL1_TAG | 0x10f030c,   17761036,  44,  47},
	{24, LCP_LEVEL1_TAG | 0x10c0a11,   17566225,  46,  49},
	{24, LCP_LEVEL1_TAG | 0x1060004,   17170436,  49,  52},
	{24, LCP_LEVEL1_TAG | 0x1090014,   17367060,  52,  55},
	{24, LCP_LEVEL1_TAG | 0x1140f14,   18091796,  54,  57},
	{24, LCP_LEVEL1_TAG | 0x1140914,   18090260,  56,  59},
	{24, LCP_LEVEL1_TAG | 0x109000a,   17367050,  59,  62},
	{24, LCP_LEVEL1_TAG | 0x10a0e03,   17436163,  61,  64},
	{24, LCP_LEVEL1_TAG | 0x103000f,   16973839,  63,  66},
	{24, LCP_LEVEL1_TAG | 0x1050414,   17105940,  66,  69},
	{32, LCP_LEVEL1_TAG | 0x2041407,   33821703,  67,  71},
	{24, LCP_LEVEL1_TAG | 0x1140711,   18089745,  69,  72},
	{24, LCP_LEVEL1_TAG | 0x10a0009,   17432585,  72,  75},
	{24, LCP_LEVEL1_TAG | 0x1150910,   18155792,  75,  78},
};
#endif

#if LCP_FIXED_CORE
static const expected_core PROTEIN_LEVEL2_CORES[] = {
	{ 6,         0x000000000000002e, 3457858291,   2,  13},
	{ 6,         0x0000000000000021, 2169992986,   8,  17},
	{ 9,         0x00000000000000d7, 2796433796,  12,  22},
	{ 6,         0x0000000000000032,   67839844,  17,  25},
	{10,         0x0000000000000273,   46115713,  22,  29},
	{ 9,         0x000000000000013b, 2219236039,  24,  35},
	{ 7,         0x000000000000004b, 3628824713,  32,  43},
	{ 6,         0x0000000000000011,   63918098,  40,  49},
	{ 8,         0x0000000000000049, 3912630473,  44,  55},
	{12,         0x00000000000008cb, 3223361699,  52,  62},
	{ 7,         0x000000000000006d,   73414458,  56,  66},
	{10,         0x0000000000000299, 1582691688,  61,  70},
	{ 9,         0x0000000000000196, 1128430323,  66,  72},
};
#else
static const expected_core PROTEIN_LEVEL2_CORES[] = {
	{ 6,         0x000000000000002e, 3457858291,   2,  13},
	{ 6,         0x0000000000000021, 2169992986,   8,  17},
	{ 9,         0x0000000000000152, 2593600787,  14,  26},
	{ 9,         0x0000000000000173, 1957220616,  19,  29},
	{ 9,         0x000000000000013b, 2219236039,  24,  35},
	{ 7,         0x000000000000004b, 3628824713,  32,  43},
	{ 6,         0x0000000000000011,   63918098,  40,  49},
	{ 8,         0x0000000000000049, 3912630473,  44,  55},
	{12,         0x00000000000008cb, 3223361699,  52,  62},
	{ 7,         0x000000000000006d,   73414458,  56,  66},
	{ 7,         0x0000000000000051,  362123574,  61,  71},
};
#endif

#if LCP_FIXED_CORE
static const expected_core PROTEIN_LEVEL3_CORES[] = {
	{ 6,         0x0000000000000031,  304380204,   8,  29},
	{ 9,         0x00000000000000f8, 2854964615,  17,  43},
	{ 9,         0x0000000000000117, 3701813736,  24,  55},
	{ 7,         0x0000000000000074, 1415562910,  44,  70},
};
#else
static const expected_core PROTEIN_LEVEL3_CORES[] = {
	{ 6,         0x0000000000000011, 3972495471,   2,  29},
	{ 9,         0x00000000000000f8, 1323032698,  14,  43},
	{ 9,         0x0000000000000117, 3701813736,  24,  55},
	{ 7,         0x0000000000000074,  932510592,  44,  71},
};
#endif

void test_lps_constructor() {

	LCP_INIT();

	std::string test_string = PROTEIN_TEST_STRING;
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	assert_cores_match(&lps_obj, PROTEIN_LEVEL1_CORES, CORE_COUNT(PROTEIN_LEVEL1_CORES), "test_lps_constructor");

	free_lps(&lps_obj);

	log("...  test_lps_constructor passed!");
}

void test_lps_deepen() {

	LCP_INIT();

	std::string test_string = PROTEIN_TEST_STRING;
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	int success = lps_deepen(&lps_obj, 2);
	assert(success && "Deepening to level 2 should be successful");
	assert(lps_obj.level == 2 && "Level should be 2 after deepening");
	assert_cores_match(&lps_obj, PROTEIN_LEVEL2_CORES, CORE_COUNT(PROTEIN_LEVEL2_CORES), "test_lps_deepen (level 2)");

	success = lps_deepen(&lps_obj, 3);
	assert(success && "Deepening to level 3 should be successful");
	assert(lps_obj.level == 3 && "Level should be 3 after deepening");
	assert_cores_match(&lps_obj, PROTEIN_LEVEL3_CORES, CORE_COUNT(PROTEIN_LEVEL3_CORES), "test_lps_deepen (level 3)");

	success = lps_deepen(&lps_obj, 3);
	assert(success == 0 && "Deepening to a level already reached should be unsuccessful");

	free_lps(&lps_obj);

	log("...  test_lps_deepen passed!");
}

// Properties of the wider layout that hold by construction, independent of the
// tables above: 8 bits per residue, sentinel set, packed length agreeing with
// bit_size, and every symbol field a valid encoding.
void test_lps_protein_layout() {

	LCP_INIT();

	std::string test_string = PROTEIN_TEST_STRING;
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	assert(lps_obj.size > 0 && "Protein parsing should produce cores");

	for (int i = 0; i < lps_obj.size; i++) {
		const struct core *cr = &(lps_obj.cores[i]);
		assert((cr->bit_rep & LCP_LEVEL1_TAG) && "Level 1 core should carry the sentinel bit");
		assert(cr->bit_size % LCP_SYMBOL_BITS == 0 && "bit_size should be a whole number of residues");
		unsigned residues = cr->bit_size / LCP_SYMBOL_BITS;
		assert(residues >= 3 && "A core spans at least three residues");
		assert(mid_count(cr->bit_rep) + 2 == residues && "Packed length should agree with bit_size");
		assert(((cr->bit_rep >> (2u * LCP_SYMBOL_BITS)) & LCP_SYMBOL_MASK) < 25 && "first residue in range");
		assert(((cr->bit_rep >> LCP_SYMBOL_BITS) & LCP_SYMBOL_MASK) < 25 && "middle residue in range");
		assert((cr->bit_rep & LCP_SYMBOL_MASK) < 25 && "last residue in range");
#if LCP_POS_BITS != 0
		assert(cr->end - cr->start == residues && "Interval should match the residue count");
#endif
	}

	free_lps(&lps_obj);

	log("...  test_lps_protein_layout passed!");
}

#else /* !LCP_ALPHABET_PROTEIN */

#if !LCP_FIXED_CORE

// level 1 cores of "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT"
// (same result as the reverse complement of "AGGACTGTGATCTCCTCACACCTGAGCTCAGCTGGCGCTTGGCTGTCGTGGGCTGGGGTCACCAGGTCCC")
static const expected_core LEVEL1_CORES[] = {
	{ 6, LCP_LEVEL1_TAG | 0b01100001,         97,   2,   5},
	{ 8, LCP_LEVEL1_TAG | 0b10000111,        135,   3,   7},
	{ 6, LCP_LEVEL1_TAG | 0b01011110,         94,   5,   8},
	{ 8, LCP_LEVEL1_TAG | 0b10111011,        187,   6,  10},
	{ 6, LCP_LEVEL1_TAG | 0b01101110,        110,   8,  11},
	{ 6, LCP_LEVEL1_TAG | 0b01100001,         97,  10,  13},
	{12, LCP_LEVEL1_TAG | 0b100000100,       260,  11,  17},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  15,  18},
	{10, LCP_LEVEL1_TAG | 0b11100100,        228,  17,  22},
	{ 6, LCP_LEVEL1_TAG | 0b01010001,         81,  20,  23},
	{ 6, LCP_LEVEL1_TAG | 0b01100001,         97,  23,  26},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  25,  28},
	{ 8, LCP_LEVEL1_TAG | 0b10100100,        164,  27,  31},
	{ 8, LCP_LEVEL1_TAG | 0b10010010,        146,  29,  33},
	{ 6, LCP_LEVEL1_TAG | 0b01100110,        102,  32,  35},
	{ 8, LCP_LEVEL1_TAG | 0b10100100,        164,  34,  38},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  36,  39},
	{ 6, LCP_LEVEL1_TAG | 0b01100111,        103,  38,  41},
	{ 6, LCP_LEVEL1_TAG | 0b01100010,         98,  41,  44},
	{ 6, LCP_LEVEL1_TAG | 0b01100111,        103,  43,  46},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  46,  49},
	{ 8, LCP_LEVEL1_TAG | 0b10001011,        139,  47,  51},
	{ 6, LCP_LEVEL1_TAG | 0b01111011,        123,  50,  53},
	{ 6, LCP_LEVEL1_TAG | 0b01100010,         98,  53,  56},
	{ 8, LCP_LEVEL1_TAG | 0b10001000,        136,  54,  58},
	{ 6, LCP_LEVEL1_TAG | 0b01100010,         98,  56,  59},
	{ 6, LCP_LEVEL1_TAG | 0b01100011,         99,  58,  61},
	{ 6, LCP_LEVEL1_TAG | 0b01010001,         81,  61,  64},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  63,  66},
	{ 6, LCP_LEVEL1_TAG | 0b01101101,        109,  65,  68},
	{ 8, LCP_LEVEL1_TAG | 0b10110111,        183,  66,  70},
};

static const expected_core LEVEL2_CORES[] = {
	{ 6,                  0b00110001, 3555632935,   2,  10},
	{ 6,                  0b00010001, 3424314286,   5,  13},
	{ 6,                  0b00010011, 1785772797,   8,  18},
	{ 8,                  0b10011000,  919618966,  15,  26},
	{ 8,                  0b10000010,  984734144,  20,  31},
	{ 7,                  0b01011101,  641965825,  25,  35},
	{ 7,                  0b01011011, 3374161539,  29,  39},
	{ 6,                  0b00010001, 4234938273,  36,  46},
	{ 6,                  0b00010001,   31892384,  41,  51},
	{ 8,                  0b10010010, 4085211265,  47,  58},
	{ 6,                  0b00110110, 1654454705,  54,  64},
	{ 6,                  0b00100001, 3806595524,  58,  68},
};

static const expected_core LEVEL3_CORES[] = {
	{ 6,                 0b000110011,  223778135,   5,  31},
	{ 6,                 0b000110111, 1340599299,  15,  39},
	{ 8,                 0b011101100, 3371074405,  25,  51},
	{ 9,                 0b110000101,  259839279,  36,  64},
};

void test_lps_constructor() {

	LCP_INIT();

	std::string test_string = "GGGACCTgGTGACCCCAGCcCACGaCAGCCAAGCGCCAGCTGAGCtCAGGTGTGAGGAGATCacaGTCCT";
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	assert_cores_match(&lps_obj, LEVEL1_CORES, CORE_COUNT(LEVEL1_CORES), "test_lps_constructor");

	free_lps(&lps_obj);

	log("...  test_lps_constructor passed!");
}

void test_lps_reverse_complement() {

	LCP_INIT();

	std::string test_string = "AGGACTgtgatCTCCTCACACCTGAGCTCAGCTGGCGCTTGGCTGTCGtGggCTGGGGTCAccAGGTCCC";
	struct lps lps_obj;
	init_lps2(&lps_obj, test_string.c_str(), test_string.size());

	assert_cores_match(&lps_obj, LEVEL1_CORES, CORE_COUNT(LEVEL1_CORES), "test_lps_reverse_complement");

	free_lps(&lps_obj);

	log("...  test_lps_reverse_complement passed!");
}

void test_lps_deepen() {

	LCP_INIT();

	std::string test_string = "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT";
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	// deepen to level 2
	int success = lps_deepen(&lps_obj, 2);
	assert(success && "Deepening to level 2 should be successful");
	assert(lps_obj.level == 2 && "Level should be 2 after deepening");
	assert_cores_match(&lps_obj, LEVEL2_CORES, CORE_COUNT(LEVEL2_CORES), "test_lps_deepen (level 2)");

	// deepen to level 3
	success = lps_deepen(&lps_obj, 3);
	assert(success && "Deepening to level 3 should be successful");
	assert(lps_obj.level == 3 && "Level should be 3 after deepening");
	assert_cores_match(&lps_obj, LEVEL3_CORES, CORE_COUNT(LEVEL3_CORES), "test_lps_deepen (level 3)");

	// attempt to deepen to a level already reached (should not do anything)
	success = lps_deepen(&lps_obj, 3);
	assert(success == 0 && "Deepening to a lower level should be unsuccessful");
	assert(lps_obj.level == 3 && "Level should be unchanged after a no-op deepen");
	assert_cores_match(&lps_obj, LEVEL3_CORES, CORE_COUNT(LEVEL3_CORES), "test_lps_deepen (no-op)");

	free_lps(&lps_obj);

	log("...  test_lps_deepen passed!");
}

#if LCP_POS_BITS != 0

void test_lps_consistency() {

	LCP_INIT();

	std::string sequence = test_sequence();

	struct lps lps_obj;
	init_lps(&lps_obj, sequence.c_str(), sequence.size());
	lps_deepen(&lps_obj, 5);

	assert(lps_obj.size > 5000 && "Test sequence should yield enough cores at level 5");

	int start = lps_obj.cores[5000].start;
	int end = lps_obj.cores[5000].end;

	// the core is identified when its own interval is parsed
	std::string subsequence1(sequence.begin() + start, sequence.begin() + end);
	struct lps lps_obj1;
	init_lps_offset(&lps_obj1, subsequence1.c_str(), subsequence1.size(), start);
	lps_deepen(&lps_obj1, 5);
	assert(lps_obj1.size > 0 && "Parsing the core's own interval should yield a core");
	assert(core_eq(&(lps_obj.cores[5000]), &(lps_obj1.cores[0])) && "Core should be identified in the original subsequence");

	// it is not identified when the interval is trimmed at the front
	std::string subsequence2(sequence.begin() + start + 1, sequence.begin() + end);
	struct lps lps_obj2;
	init_lps_offset(&lps_obj2, subsequence2.c_str(), subsequence2.size(), start + 1);
	lps_deepen(&lps_obj2, 5);
	assert(lps_obj2.size == 0 && "Core should not be identified when trimmed at the front");

	// nor when it is trimmed at the back
	std::string subsequence3(sequence.begin() + start, sequence.begin() + end - 1);
	struct lps lps_obj3;
	init_lps_offset(&lps_obj3, subsequence3.c_str(), subsequence3.size(), start);
	lps_deepen(&lps_obj3, 5);
	assert(lps_obj3.size == 0 && "Core should not be identified when trimmed at the back");

	free_lps(&lps_obj);
	free_lps(&lps_obj1);
	free_lps(&lps_obj2);
	free_lps(&lps_obj3);

	log("...  test_lps_consistency passed!");
}

#endif /* LCP_POS_BITS != 0 */

#else /* LCP_FIXED_CORE */

// Golden tables for CORE=fixed
static const expected_core FIXED_LEVEL1_CORES[] = {
	{ 6, LCP_LEVEL1_TAG | 0b01100001,         97,   2,   5},
	{ 6, LCP_LEVEL1_TAG | 0b01000101,         69,   3,   6},
	{ 6, LCP_LEVEL1_TAG | 0b01010111,         87,   4,   7},
	{ 6, LCP_LEVEL1_TAG | 0b01011110,         94,   5,   8},
	{ 6, LCP_LEVEL1_TAG | 0b01111010,        122,   6,   9},
	{ 6, LCP_LEVEL1_TAG | 0b01101011,        107,   7,  10},
	{ 6, LCP_LEVEL1_TAG | 0b01101110,        110,   8,  11},
	{ 6, LCP_LEVEL1_TAG | 0b01100001,         97,  10,  13},
	{ 6, LCP_LEVEL1_TAG | 0b01000101,         69,  11,  14},
	{ 6, LCP_LEVEL1_TAG | 0b01010101,         85,  12,  15},
	{ 6, LCP_LEVEL1_TAG | 0b01010101,         85,  13,  16},
	{ 6, LCP_LEVEL1_TAG | 0b01010100,         84,  14,  17},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  15,  18},
	{ 6, LCP_LEVEL1_TAG | 0b01100101,        101,  17,  20},
	{ 6, LCP_LEVEL1_TAG | 0b01010101,         85,  18,  21},
	{ 6, LCP_LEVEL1_TAG | 0b01010100,         84,  19,  22},
	{ 6, LCP_LEVEL1_TAG | 0b01010001,         81,  20,  23},
	{ 6, LCP_LEVEL1_TAG | 0b01100001,         97,  23,  26},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  25,  28},
	{ 6, LCP_LEVEL1_TAG | 0b01100101,        101,  27,  30},
	{ 6, LCP_LEVEL1_TAG | 0b01010100,         84,  28,  31},
	{ 6, LCP_LEVEL1_TAG | 0b01010000,         80,  29,  32},
	{ 6, LCP_LEVEL1_TAG | 0b01000010,         66,  30,  33},
	{ 6, LCP_LEVEL1_TAG | 0b01100110,        102,  32,  35},
	{ 6, LCP_LEVEL1_TAG | 0b01100101,        101,  34,  37},
	{ 6, LCP_LEVEL1_TAG | 0b01010100,         84,  35,  38},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  36,  39},
	{ 6, LCP_LEVEL1_TAG | 0b01100111,        103,  38,  41},
	{ 6, LCP_LEVEL1_TAG | 0b01100010,         98,  41,  44},
	{ 6, LCP_LEVEL1_TAG | 0b01100111,        103,  43,  46},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  46,  49},
	{ 6, LCP_LEVEL1_TAG | 0b01001010,         74,  47,  50},
	{ 6, LCP_LEVEL1_TAG | 0b01101011,        107,  48,  51},
	{ 6, LCP_LEVEL1_TAG | 0b01111011,        123,  50,  53},
	{ 6, LCP_LEVEL1_TAG | 0b01100010,         98,  53,  56},
	{ 6, LCP_LEVEL1_TAG | 0b01001010,         74,  54,  57},
	{ 6, LCP_LEVEL1_TAG | 0b01101000,        104,  55,  58},
	{ 6, LCP_LEVEL1_TAG | 0b01100010,         98,  56,  59},
	{ 6, LCP_LEVEL1_TAG | 0b01100011,         99,  58,  61},
	{ 6, LCP_LEVEL1_TAG | 0b01010001,         81,  61,  64},
	{ 6, LCP_LEVEL1_TAG | 0b01010010,         82,  63,  66},
	{ 6, LCP_LEVEL1_TAG | 0b01101101,        109,  65,  68},
	{ 6, LCP_LEVEL1_TAG | 0b01110101,        117,  66,  69},
	{ 6, LCP_LEVEL1_TAG | 0b01010111,         87,  67,  70},
};

static const expected_core FIXED_LEVEL2_CORES[] = {
	{ 7,               0b00001100100,  228420930,   3,   9},
	{ 6,               0b00000010001, 3052224931,   6,  13},
	{ 9,               0b00011011001, 3108203892,   8,  15},
	{11,               0b10110011100, 2696927458,  10,  16},
	{ 8,               0b00011000011, 1485993061,  12,  18},
	{ 8,               0b00011011001,  788648406,  14,  21},
	{ 8,               0b00010010001, 1452379690,  17,  23},
	{ 8,               0b00010000001, 4088632988,  20,  30},
	{ 7,               0b00000100100,  525311635,  25,  32},
	{ 8,               0b00010011101, 3887058755,  28,  35},
	{ 6,               0b00000010011,  405967400,  32,  39},
	{ 6,               0b00000010001, 4234938273,  36,  46},
	{ 7,               0b00000100111, 2905066031,  41,  50},
	{ 9,               0b00111011001, 4122272571,  46,  53},
	{ 9,               0b00100100111, 3190073306,  48,  57},
	{ 7,               0b00001111011,  138128131,  53,  59},
	{ 6,               0b00000110110, 1654454705,  55,  64},
	{ 6,               0b00000100001, 3806595524,  58,  68},
};

static const expected_core FIXED_LEVEL3_CORES[] = {
	{ 7,                 0b001110001, 3545232573,   6,  18},
	{ 7,                 0b000110110,   10338784,  10,  23},
	{ 9,                 0b101101000, 2596542711,  12,  30},
	{ 8,                 0b010000001, 2632502658,  17,  35},
	{ 6,                 0b000111011,  268756671,  28,  50},
	{ 6,                 0b000111011, 1672441688,  36,  57},
	{ 7,                 0b001000001, 4063200252,  48,  68},
};

void test_lps_constructor() {

	LCP_INIT();

	std::string test_string = "GGGACCTgGTGACCCCAGCcCACGaCAGCCAAGCGCCAGCTGAGCtCAGGTGTGAGGAGATCacaGTCCT";
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	assert_cores_match(&lps_obj, FIXED_LEVEL1_CORES, CORE_COUNT(FIXED_LEVEL1_CORES), "test_lps_constructor");

	free_lps(&lps_obj);

	log("...  test_lps_constructor passed!");
}

void test_lps_reverse_complement() {

	LCP_INIT();

	std::string test_string = "AGGACTgtgatCTCCTCACACCTGAGCTCAGCTGGCGCTTGGCTGTCGtGggCTGGGGTCAccAGGTCCC";
	struct lps lps_obj;
	init_lps2(&lps_obj, test_string.c_str(), test_string.size());

	assert_cores_match(&lps_obj, FIXED_LEVEL1_CORES, CORE_COUNT(FIXED_LEVEL1_CORES), "test_lps_reverse_complement");

	free_lps(&lps_obj);

	log("...  test_lps_reverse_complement passed!");
}

void test_lps_deepen() {

	LCP_INIT();

	std::string test_string = "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT";
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	int success = lps_deepen(&lps_obj, 2);
	assert(success && "Deepening to level 2 should be successful");
	assert(lps_obj.level == 2 && "Level should be 2 after deepening");
	assert_cores_match(&lps_obj, FIXED_LEVEL2_CORES, CORE_COUNT(FIXED_LEVEL2_CORES), "test_lps_deepen (level 2)");

	success = lps_deepen(&lps_obj, 3);
	assert(success && "Deepening to level 3 should be successful");
	assert(lps_obj.level == 3 && "Level should be 3 after deepening");
	assert_cores_match(&lps_obj, FIXED_LEVEL3_CORES, CORE_COUNT(FIXED_LEVEL3_CORES), "test_lps_deepen (level 3)");

	success = lps_deepen(&lps_obj, 3);
	assert(success == 0 && "Deepening to a lower level should be unsuccessful");
	assert(lps_obj.level == 3 && "Level should be unchanged after a no-op deepen");
	assert_cores_match(&lps_obj, FIXED_LEVEL3_CORES, CORE_COUNT(FIXED_LEVEL3_CORES), "test_lps_deepen (no-op)");

	free_lps(&lps_obj);

	log("...  test_lps_deepen passed!");
}

// Structural invariants, independent of the tables above: with every emitted
// core a length-3 window, level 1 cores span exactly 3 positions and arrive in
// non-decreasing order of start.
void test_lps_fixed_cores() {

	LCP_INIT();

	std::string test_string = "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT";
	struct lps lps_obj;
	init_lps(&lps_obj, test_string.c_str(), test_string.size());

	assert(lps_obj.size > 0 && "Fixed-length parsing should produce cores");

	for (int i = 0; i < lps_obj.size; i++) {
		const struct core *cr = &(lps_obj.cores[i]);
		assert(cr->bit_size == 6 && "Level 1 fixed core should encode 3 symbols");
		assert((cr->bit_rep & LCP_LEVEL1_TAG) && "Level 1 core should carry the sentinel bit");
#if LCP_POS_BITS != 0
		assert(cr->end - cr->start == 3 && "Fixed core should span 3 positions");
		if (i > 0) {
			assert(lps_obj.cores[i-1].start <= cr->start && "Cores should be ordered by start");
		}
#endif
	}

	free_lps(&lps_obj);

	log("...  test_lps_fixed_cores passed!");
}

#endif /* !LCP_FIXED_CORE */

#endif /* LCP_ALPHABET_PROTEIN */

// ---------------------------------------------------------------------------
// Mode-independent tests: these compare the library against itself rather than
// against a golden table, so they hold for CORE=var and CORE=fixed alike.
// ---------------------------------------------------------------------------

void test_lps_parallel() {

	LCP_INIT();

	std::string test_string = "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT";

	struct lps lps_obj1;
	init_lps(&lps_obj1, test_string.c_str(), test_string.size());
	lps_deepen(&lps_obj1, 3);

	struct lps lps_obj2;
	init_lps(&lps_obj2, test_string.c_str(), test_string.size());
	lps_deepen_parallel(&lps_obj2, 3, 4);

	assert(lps_obj1.size == lps_obj2.size && "Parallel deepen should produce the same core count");
	for (int i = 0; i < lps_obj1.size; i++) {
		assert(core_eq(&(lps_obj1.cores[i]), &(lps_obj2.cores[i])) && "Parallel deepen should produce the same cores");
	}

	free_lps(&lps_obj1);
	free_lps(&lps_obj2);

	log("...  test_lps_parallel passed!");
}

#if LCP_POS_BITS != 0

void test_lps_split_init() {

	LCP_INIT();

	std::string sequence = test_sequence();

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

#endif /* LCP_POS_BITS != 0 */

#if LCP_POS_BITS != 0

void test_lps_core_spans() {

	LCP_INIT();

	std::string sequence = test_sequence();

	struct lps forward;
	init_lps(&forward, sequence.c_str(), sequence.size());
#if !LCP_ALPHABET_PROTEIN
	struct lps revcomp;
	init_lps2(&revcomp, sequence.c_str(), sequence.size());
	struct lps *objs[2] = { &forward, &revcomp };
	const char *names[2] = { "init_lps (parse1)", "init_lps2 (parse2)" };
	const int count = 2;
#else
	struct lps *objs[1] = { &forward };
	const char *names[1] = { "init_lps (parse1)" };
	const int count = 1;
#endif

	for (int k = 0; k < count; k++) {
		assert(objs[k]->size > 0 && "Test sequence should produce cores");
		for (int i = 0; i < objs[k]->size; i++) {
			const struct core *cr = &(objs[k]->cores[i]);
			if (cr->end - cr->start != cr->bit_size / LCP_SYMBOL_BITS) {
				std::cerr << names[k] << ": core " << i << " spans ["
				          << cr->start << ", " << cr->end << ") = "
				          << (cr->end - cr->start) << " positions but encodes "
				          << (cr->bit_size / LCP_SYMBOL_BITS) << " symbols" << std::endl;
				assert(false && "Core interval should match its symbol count");
			}
		}
	}

	free_lps(&forward);
#if !LCP_ALPHABET_PROTEIN
	free_lps(&revcomp);
#endif

	log("...  test_lps_core_spans passed!");
}

#endif /* LCP_POS_BITS != 0 */

#if LCP_LABEL_BITS == 64

void test_lps_label_width() {

	LCP_INIT();

	std::string sequence = test_sequence();

	struct lps lps_obj;
	init_lps(&lps_obj, sequence.c_str(), sequence.size());

	for (int level = 2; level <= 4; level++) {
		lps_deepen1(&lps_obj);
		long high_set = 0;
		for (int i = 0; i < lps_obj.size; i++) {
			if (lps_obj.cores[i].label >> 32) {
				high_set++;
			}
		}
		if (high_set == 0) {
			std::cerr << "level " << level << ": no core uses bits 32..63 of its label" << std::endl;
			assert(false && "64-bit labels should populate their high half");
		}
	}

	free_lps(&lps_obj);

	log("...  test_lps_label_width passed!");
}

#endif /* LCP_LABEL_BITS == 64 */

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
	if (!in) {
		perror("Failed to open file");
		exit(EXIT_FAILURE);
	}
	struct lps lps_obj_from_file;
	init_lps3(&lps_obj_from_file, in);
	fclose(in);

	// compare the read object with the original
	assert(lps_obj.level == lps_obj_from_file.level && "Level should match after reading from file");
	assert(lps_obj.size == lps_obj_from_file.size && "Core size should match after reading from file");
	for (int i = 0; i < lps_obj.size; i++) {
		const struct core *a = &(lps_obj.cores[i]), *b = &(lps_obj_from_file.cores[i]);
		assert(core_eq(a, b) && "Cores should match after reading from file");
		assert(a->bit_size == b->bit_size && "bit_size should survive the round trip");
#if LCP_LABEL_BITS != 0
		assert(a->label == b->label && "label should survive the round trip");
#endif
#if LCP_POS_BITS != 0
		assert(a->start == b->start && a->end == b->end && "positions should survive the round trip");
#endif
	}

	// clean up the test file
	std::remove("lps_test.bin");

	free_lps(&lps_obj);
	free_lps(&lps_obj_from_file);

	log("...  test_lps_file_io passed!");
}

int main() {

	log("Running test_lps...");

	// golden-table tests; both modes have their own tables
	test_lps_constructor();
#if !LCP_ALPHABET_PROTEIN
	test_lps_reverse_complement();
#endif
	test_lps_deepen();
#if LCP_ALPHABET_PROTEIN
	test_lps_protein_layout();
#elif LCP_FIXED_CORE
	test_lps_fixed_cores();
#elif LCP_POS_BITS != 0
	test_lps_consistency();
#endif
	test_lps_parallel();
#if LCP_POS_BITS != 0
	test_lps_split_init();
	test_lps_core_spans();
#endif
#if LCP_LABEL_BITS == 64
	test_lps_label_width();
#endif
	test_lps_file_io();

	log("All tests in test_lps completed successfully!");
	return 0;
}
