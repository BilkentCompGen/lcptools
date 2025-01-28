#include "lps.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1024
#define SEQUENCE_CAPACITY 250000000

void print_usage(const char *lcptools) {
	printf("Usage: %s <command> <filename> <lcp-level> [sequence-size]\n", lcptools);
	printf("Commands:\n");
	printf("  falcpt   Process the fasta file.\n");
	// printf("  fqlcpt   Process the fasta file.\n");
	printf("File extensions:\n");
	printf("  .fasta, .fa, .fastq, .fq\n");
}

int validate_extension(const char *infilename) {
    const char *valid_extensions[] = {".fasta", ".fa", ".fastq", ".fq"};
    uint64_t infilename_len = strlen(infilename);

    for (int i = 0; i < 4; i++) {
        uint64_t ext_len = strlen(valid_extensions[i]);
        if (ext_len < infilename_len &&
            strcmp(infilename + infilename_len - ext_len, valid_extensions[i]) == 0) {
            return 0;
        }
    }
    return -1;
}

int isNumber(const char *str) {
    if (str == NULL || *str == '\0') {
        return -1;
    }

    char *end;
    errno = 0;
    strtol(str, &end, 10);

    // Check if there was an overflow or no valid conversion
    if (errno != 0 || *end != '\0') {
        return 1;
    }

    return 0;
}

void done(FILE *out) {
    char isDone = 0;
    fwrite(&isDone, 1, 1, out);
}

int process_fasta(const char *infilename, const char *outfilename, int lcp_level, long unsigned int sequence_size) {

    FILE *infile = fopen(infilename, "rb");
    FILE *outfile = fopen(outfilename, "wb");

	if (!infile || !outfile) {
        fprintf(stderr, "Error opening file\n");
        if (infile) fclose(infile);
        if (outfile) fclose(outfile);
		return 1;
	}

    // Allocate a buffer for the sequence
	char *sequence = (char *)malloc(sequence_size + 1);
    if (!sequence) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(infile);
        fclose(outfile);
        return 1;
    }
    sequence[0] = '\0';

    // Allocate a buffer for reading lines
    char line[1024];

	// Initialize lcp encoding
    LCP_INIT(0);

	while (fgets(line, sizeof(line), infile)) {

        // Remove newline character at the end of the line
        line[strcspn(line, "\r\n")] = '\0';

		if (line[0] != '>') {
            if (strlen(sequence) + strlen(line) >= sequence_size) {
                fprintf(stderr, "Error: Sequence exceeds buffer size\n");
                free(sequence);
                fclose(infile);
                fclose(outfile);
                return 1;
            }
			strncat(sequence, line, sequence_size - strlen(sequence) - 1);
			continue;
		}

		// Process previous chromosome before moving into new one
		if (strlen(sequence) > 0) {
			struct lps str;
            init_lps(&str, sequence, strlen(sequence));
            lps_deepen(&str, lcp_level);

            write_lps(&str, outfile);

			free_lps(&str);

			sequence[0] = '\0';
		}
	}

	if (strlen(sequence) > 0) {
        struct lps str;
        init_lps(&str, sequence, strlen(sequence));
        lps_deepen(&str, lcp_level);
			
		write_lps(&str, outfile);

		free_lps(&str);
	}

	done(outfile);

    free(sequence);
    fclose(infile);
    fclose(outfile);

	return 0;
}

int main(int argc, char *argv[]) {

	if (argc < 4) {
		print_usage(argv[0]);
		return 1;
	}

	const char *command = argv[1];
	const char *infilename = argv[2];

	if (strcmp(command, "falcpt") != 0) {
		fprintf(stderr, "Error: Unsupported command %s\n", command);
		print_usage(argv[0]);
		return 1;
	}

	if (validate_extension(infilename)) {
		fprintf(stderr, "Error: Invalid file extension. Supported extensions are .fasta, .fa, .fastq, .fq\n");
		return 1;
	}

	if (isNumber(argv[3])) {
		fprintf(stderr, "Error: The lcp level argument must be a positive integer.\n");
		return 1;
	}

	int lcp_level = atol(argv[3]);
	long unsigned int sequence_size = SEQUENCE_CAPACITY;

	if (argc == 5) {
		if (!isNumber(argv[4])) {
			fprintf(stderr, "Error: The sequence size argument must be a positive integer.\n");
			return 1;
		}
		sequence_size = atoi(argv[4]);
	}

	// generate output infilename
    char outfilename[1024];
    snprintf(outfilename, sizeof(outfilename), "%s.lcpt", infilename);

	printf("Output: %s\n", outfilename);

	process_fasta(infilename, outfilename, lcp_level, sequence_size);
	// todo: fastq

	return 0;
}