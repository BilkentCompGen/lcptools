#include "encoding.h"


int alphabet[128];
#if !LCP_ALPHABET_PROTEIN
int rc_alphabet[128];
#endif
char characters[128];

#if LCP_ALPHABET_PROTEIN
/** 
 * The 20 standard amino acids plus the ambiguity and non-standard codes
 * B (Asx), Z (Glx), X (any), U (selenocysteine) and O (pyrrolysine), in 
 * alphabetical order over the whole set, 0..24.
 */
static const char PROTEIN_SYMBOLS[] = "ABCDEFGHIKLMNOPQRSTUVWXYZ";
#endif

void LCP_SUMMARY(void) {
    printf("# Alphabet encoding summary\n");
    printf("# Coefficients: ");
    for (int i = 0; i < 128; i++) {
        if (alphabet[i] != -1) {
            printf("%c:%d ", i, alphabet[i]);
        }
    }
    printf("\n");
}

void LCP_INIT(void) {
    LCP_INIT2(0);
}

void LCP_INIT2(int verbose) {

    for (int current_index = 0; current_index < 128; current_index++) {
        alphabet[current_index] = -1;
        characters[current_index] = 126;
#if !LCP_ALPHABET_PROTEIN
        rc_alphabet[current_index] = -1;
#endif
    }

#if LCP_ALPHABET_PROTEIN
    for (int i = 0; PROTEIN_SYMBOLS[i] != '\0'; i++) {
        char upper = PROTEIN_SYMBOLS[i];
        char lower = (char)(upper - 'A' + 'a');
        alphabet[(int)upper] = i;
        alphabet[(int)lower] = i;
        characters[i] = upper;
    }
#else
    // init coefficients A/a=0, T/t=3, G/g=2, C/c=1
    alphabet['A'] = 0; alphabet['a'] = 0;
    alphabet['T'] = 3; alphabet['t'] = 3;
    alphabet['G'] = 2; alphabet['g'] = 2;
    alphabet['C'] = 1; alphabet['c'] = 1;

    rc_alphabet['A'] = 3; rc_alphabet['a'] = 3;
    rc_alphabet['T'] = 0; rc_alphabet['t'] = 0;
    rc_alphabet['G'] = 1; rc_alphabet['g'] = 1;
    rc_alphabet['C'] = 2; rc_alphabet['c'] = 2;

    characters[0] = 'A';
    characters[1] = 'C';
    characters[2] = 'G';
    characters[3] = 'T';
#endif

    if (verbose)
        LCP_SUMMARY();
}

int LCP_INIT_FILE(const char *encoding_file, int verbose) {
    
    FILE *encodings = fopen(encoding_file, "r");
    if (!encodings) {
        if (verbose) {
            fprintf(stderr, "Error: Could not open file %s\n", encoding_file);
        }
        return -1;
    }

    // clear arrays
    for (int current_index = 0; current_index < 128; current_index++) {
        alphabet[current_index] = -1;
        characters[current_index] = 126;
#if !LCP_ALPHABET_PROTEIN
        rc_alphabet[current_index] = -1;
#endif
    }

    char character;
    int encoding, mx = -1;

#if LCP_ALPHABET_PROTEIN
    // protein has no complement, so each line is "<symbol> <encoding>"
    while (fscanf(encodings, " %c %d", &character, &encoding) == 2) {
        alphabet[(unsigned char)character] = encoding;
        mx = maximum(encoding, mx);
    }
#else
    int rev_encoding;
    while (fscanf(encodings, " %c %d %d", &character, &encoding, &rev_encoding) == 3) {
        alphabet[(unsigned char)character] = encoding;
        rc_alphabet[(unsigned char)character] = rev_encoding;

        mx = maximum(encoding, mx);
        mx = maximum(rev_encoding, mx);
    }
#endif

    fclose(encodings);

    int bit_count = 0;
    while (mx > 0) {
        bit_count++;
        mx = mx / 2;
    }

    if (bit_count > LCP_SYMBOL_BITS) {
        fprintf(stderr, "Your alphabet needs %d binary digits per symbol, but this build allows at most %d. Rebuild with a wider ALPHABET.\n", bit_count, (int)LCP_SYMBOL_BITS);
        exit(EXIT_FAILURE);
    }

    return 0;
}