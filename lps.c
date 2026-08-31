#include "lps.h"

#define minimum(a, b) ((a) < (b) ? (a) : (b))

#if LCP_FIXED_CORE
#include "emit_fixed.h"
#else
#include "emit_var.h"
#endif


/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
int parse1(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = begin;
    const char *it2 = end;
    int core_index = 0;
    int last_invalid_char_index = -1;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (alphabet[(unsigned char)*it1] == -1) {
            last_invalid_char_index = it1 - begin;
            continue;
        }

        if (alphabet[(unsigned char)*it1] == alphabet[(unsigned char)*(it1+1)]) {
            continue;
        }

        // check for RINT core
        if (alphabet[(unsigned char)*(it1+1)] == alphabet[(unsigned char)*(it1+2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 + 2;
            while (temp < end && alphabet[(unsigned char)*(temp-1)] == alphabet[(unsigned char)*temp]) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                    core_index = emit_core1(cores, core_index, begin, offset, it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                core_index = emit_core1(cores, core_index, begin, offset, it1, 2+middle_count, it1-begin+offset, it2-begin+offset);

                continue;
            }
        }

        if (alphabet[(unsigned char)*it1] > alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] < alphabet[(unsigned char)*(it1+2)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                core_index = emit_core1(cores, core_index, begin, offset, it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            alphabet[(unsigned char)*it1] < alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] > alphabet[(unsigned char)*(it1+2)] &&
            alphabet[(unsigned char)*(it1-1)] <= alphabet[(unsigned char)*(it1)] &&
            alphabet[(unsigned char)*(it1+2)] >= alphabet[(unsigned char)*(it1+3)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                core_index = emit_core1(cores, core_index, begin, offset, it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

#if !LCP_ALPHABET_PROTEIN

/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores using complement alphabet.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order based on complement), 
 * it builds and stores these cores for further processing in the LCP framework.
 *
 * Positions are reverse-complement indices: a character at address p sits at
 * index (end - p - 1). Every core here therefore reports
 *
 *     start = end - seg_begin - 1        end = start + seg_len
 *
 * where seg_begin is the segment's highest address.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
int parse2(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = end - 1;
    const char *it2 = begin - 1;
    int core_index = 0;

    // find lcp cores
    for (; begin <= it1 - 2; it1--) {

        // skip invalid character
        if (rc_alphabet[(unsigned char)*it1] == rc_alphabet[(unsigned char)*(it1-1)]) {
            continue;
        }

        // check for RINT core
        if (rc_alphabet[(unsigned char)*(it1-1)] == rc_alphabet[(unsigned char)*(it1-2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 - 2;
            while (begin <= temp && rc_alphabet[(unsigned char)*(temp+1)] == rc_alphabet[(unsigned char)*temp]) {
                temp--;
                middle_count++;
            }
            if (begin <= temp) {
                // check if there is any SSEQ cores left behind
                if (it1 < it2) {
                    core_index = emit_core2(cores, core_index, end, offset, it2+1, it2-it1+2, end-it2-2+offset, end-it1+offset);
                }

                // create RINT core
                it2 = it1 - 2 - middle_count;
                core_index = emit_core2(cores, core_index, end, offset, it1, 2+middle_count, end-it1-1+offset, end-it2-1+offset);

                continue;
            }
        }

        if (rc_alphabet[(unsigned char)*it1] > rc_alphabet[(unsigned char)*(it1-1)] &&
            rc_alphabet[(unsigned char)*(it1-1)] < rc_alphabet[(unsigned char)*(it1-2)]) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                core_index = emit_core2(cores, core_index, end, offset, it2+1, it2-it1+2, end-it2-2+offset, end-it1+offset);
            }

            // create LMIN core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (begin <= it1-3 &&
            rc_alphabet[(unsigned char)*it1] < rc_alphabet[(unsigned char)*(it1-1)] &&
            rc_alphabet[(unsigned char)*(it1-1)] > rc_alphabet[(unsigned char)*(it1-2)] &&
            rc_alphabet[(unsigned char)*(it1+1)] <= rc_alphabet[(unsigned char)*(it1)] &&
            rc_alphabet[(unsigned char)*(it1-2)] >= rc_alphabet[(unsigned char)*(it1-3)]) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                core_index = emit_core2(cores, core_index, end, offset, it2+1, it2-it1+2, end-it2-2+offset, end-it1+offset);
            }

            // create LMAX core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

#endif /* !LCP_ALPHABET_PROTEIN */

/**
 * @brief Parses a array of cores to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a array of `core` structures defined by iterators `begin` and `end` and 
 * identifies key segments, called "cores," that represent the (LCP) regions. By analyzing
 * `core` structure relationships in the array (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the `core` array to parse.
 * @param end Iterator pointing to the end of the `core` array to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @return Size of the cores identified in the given string.
 */
int parse3(struct core *begin, struct core *end, struct core *cores) {

    struct core *it1 = begin;
    struct core *it2 = end;
    int core_index = 0;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (core_eq(it1, it1+1)) {
            continue;
        }

        // check for RINT core
        if (core_eq(it1+1, it1+2)) {

            // count middle characters
            uint32_t middle_count = 1;
            struct core *temp = it1 + 2;
            while (temp < end && core_eq(temp-1, temp)) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    core_index = emit_core3(cores, core_index, it2-1, it1-it2+2);
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                core_index = emit_core3(cores, core_index, it1, it2-it1);

                continue;
            }
        }

        // check for LMIN
        if (core_gt(it1, it1+1) && core_lt(it1+1, it1+2)) {
            
            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                core_index = emit_core3(cores, core_index, it2-1, it1-it2+2);
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            core_lt(it1, it1+1) &&
            core_gt(it1+1, it1+2) &&
            core_leq(it1-1, it1) &&
            core_geq(it1+2, it1+3)) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                core_index = emit_core3(cores, core_index, it2-1, it1-it2+2);
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
            core_index++;

            continue;
        }
    }
    return core_index;
}

void init_lps(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/LCP_CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, 0);
}

void init_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/LCP_CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, offset);
}

#if !LCP_ALPHABET_PROTEIN

void init_lps2(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/LCP_CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse2(str, str+len, lps_ptr->cores, 0);
}

#endif /* !LCP_ALPHABET_PROTEIN */

void init_lps3(struct lps *lps_ptr, FILE *in) {
    // read the level from the binary file
    if (fread(&(lps_ptr->level), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading level from file\n");
        exit(EXIT_FAILURE);
    }

    // read the size (number of cores)
    if(fread(&(lps_ptr->size), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading size from file\n");
        exit(EXIT_FAILURE);
    }

    lps_ptr->cores = NULL;

    if (lps_ptr->size) {
        // allocate memory for the cores array
        lps_ptr->cores = (struct core *)malloc(lps_ptr->size * sizeof(struct core));
        if (fread(lps_ptr->cores, lps_ptr->size * sizeof(struct core), 1, in) != 1) {
            fprintf(stderr, "Error reading cores from file\n");
            exit(EXIT_FAILURE);
        }
    }
}

#if LCP_POS_BITS != 0

void init_lps4(struct lps *lps_ptr, const char *str, int len, int lcp_level, int chunk_size) {

    if (lcp_level < 1)
        return;

    lps_ptr->level = 1;
    lps_ptr->size = 0; 
    int estimated_size = (int)(len / pow((double)LCP_CONSTANT_FACTOR, lcp_level));
    lps_ptr->cores = (struct core *)malloc(estimated_size*sizeof(struct core));

    int str_index = 0, core_index = 0;

    {
        int str_len = minimum(chunk_size, len);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str, str_len, 0);
        lps_deepen(&temp_lps, lcp_level);

        if (temp_lps.size) {
            memcpy(lps_ptr->cores, temp_lps.cores, (temp_lps.size)*sizeof(struct core));
            core_index = (temp_lps.size);
            lps_ptr->size = (temp_lps.size);
            if (temp_lps.size>1)
                str_index = lps_ptr->cores[core_index-2].start;
            else 
                str_index = lps_ptr->cores[core_index-1].start;
        }
        free(temp_lps.cores);
    }

    while (str_index < len) {
        int str_len = minimum(chunk_size, len-str_index);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str+str_index, str_len, str_index);
        lps_deepen(&temp_lps, lcp_level);

        if (1<temp_lps.size) {
            int overlap = 2;
            while (0<overlap) {
                if (lps_ptr->cores[core_index-overlap].start == temp_lps.cores[0].start)
                    break;
                overlap--;
            }
            memcpy(lps_ptr->cores+core_index, temp_lps.cores+overlap, (temp_lps.size-overlap)*sizeof(struct core));
            core_index += (temp_lps.size-overlap);
            lps_ptr->size += (temp_lps.size-overlap);

            if ((uint64_t)str_index < lps_ptr->cores[core_index-2].start) {
                str_index = lps_ptr->cores[core_index-2].start;
                free(temp_lps.cores);
                continue;
            } 
        }
        
        // find next start point
        for(int i=str_index+str_len-1; str_index <= i; i--) {
            if (alphabet[(unsigned char)*(str+i)] == -1) {
                str_index = i+1;
                break;
            }
        }
        if (alphabet[(unsigned char)*(str+str_index)] != -1) { // all of the characters are valid, so not valid cores found
            str_index += str_len;
        }
        
        free(temp_lps.cores);
    }

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));
}

#endif

void free_lps(struct lps *lps_ptr) {
    free(lps_ptr->cores);
    lps_ptr->size = 0;
}

void write_lps(struct lps *lps_ptr, FILE *out) {
    // write the level field
    fwrite(&(lps_ptr->level), sizeof(int), 1, out);

    // write the size (number of cores)
    fwrite(&(lps_ptr->size), sizeof(int), 1, out);

    // write each core object iteratively
    if (lps_ptr->size) {
        fwrite(lps_ptr->cores, lps_ptr->size*sizeof(struct core), 1, out);
    }
}

int64_t lps_memsize(const struct lps *lps_ptr) {
    return sizeof(struct lps) + lps_ptr->size * sizeof(struct core);
}

/**
 * @brief Performs Deterministic Coin Tossing (DCT) compression on binary sequences.
 *
 * This function is a central part of the LCP (Locally Consisted Parsing) algorithm. It identifies differences
 * between consecutive binary strings, compressing the information by focusing on the position and value of
 * the first divergent bit from the right-end of the strings. This difference is used to generate a compact
 * 'core' that encapsulates the unique elements of each sequence in a smaller binary form.
 *
 * This compression significantly reduces redundant information, making further analysis of the sequences
 * within the LCP framework more efficient and manageable.
 *
 * @return 0 if dct is performed, -1 if no enough cores are available for dct.
 */
int lcp_dct(struct lps *lps_ptr) {

    // at least 2 cores are needed for compression
    if (lps_ptr->size < LCP_DCT_ITERATION_COUNT + 1) {
        return -1;
    }

    if (lps_ptr->level == 1) {
        for (uint64_t dct_index = 0; dct_index < LCP_DCT_ITERATION_COUNT; dct_index++) {
            struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

            for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
                core_compress_level1(it_left, it_right);
            }
        }
    } else {
        for (uint64_t dct_index = 0; dct_index < LCP_DCT_ITERATION_COUNT; dct_index++) {
            struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

            for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
                core_compress_upper(it_left, it_right);
            }
        }
    }

    return 0;
}

int lps_deepen1(struct lps *lps_ptr) {

    // compress cores
    if (lcp_dct(lps_ptr) < 0) {
        lps_ptr->size = 0;
        lps_ptr->level++;
        return 0;
    }

    // find new cores
#if !LCP_FIXED_CORE
    int new_size = parse3(lps_ptr->cores + LCP_DCT_ITERATION_COUNT, lps_ptr->cores + lps_ptr->size, lps_ptr->cores);
#else
    struct core *new_cores = (struct core *)malloc((lps_ptr->size - LCP_DCT_ITERATION_COUNT) * sizeof(struct core));
    if (new_cores == NULL) {
        return 0;
    }

    int new_size = parse3(lps_ptr->cores + LCP_DCT_ITERATION_COUNT, lps_ptr->cores + lps_ptr->size, new_cores);

    // replace the previous level; struct core owns no heap memory
    free(lps_ptr->cores);
    lps_ptr->cores = new_cores;
#endif
    lps_ptr->size = new_size;

    lps_ptr->level++;

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));

    return 1;
}

int lps_deepen(struct lps *lps_ptr, int lcp_level) {

    if (lcp_level <= lps_ptr->level)
        return 0;

    while (lps_ptr->level < lcp_level && lps_deepen1(lps_ptr))
        ;

    return 1;
}

typedef struct {
    struct core *cores;     // destination cores (real array)
    struct core dummy_left; // uncompressed core at begin-1
    int offset_begin;       // inclusive i
    int offset_end;         // exclusive i
    int flags;              // 1 bit (process head) 1 bit (level 1 or not)
} dct_worker_args_t;

/**
 * @brief Worker routine for parallel DCT compression.
 *
 * This function is executed by a single POSIX thread and performs
 * compression on a contiguous subrange of cores. The behavior depends
 * on the `flags` field:
 *
 * - flags == -1 : No work is performed (thread exits immediately).
 * - flags % 2   : Use core_compress_level1().
 * - otherwise   : Use core_compress_upper().
 *
 * If flags > 1, the worker also performs an additional boundary
 * compression using `dummy_left` and the first core in its range.
 * This is required when the thread's assigned range does not begin
 * at the global DCT starting index, ensuring correctness across
 * chunk boundaries.
 *
 * The function assumes that:
 * - Each thread operates on a disjoint range of cores.
 * - Any required boundary state is provided through `dummy_left`.
 *
 * @param argp Pointer to a dct_worker_args_t structure containing
 *             the thread parameters.
 *
 * @return Always returns NULL (required by pthread signature).
 */
static void *dct_worker(void *argp) {
    dct_worker_args_t *dct_params = (dct_worker_args_t *)argp;

    if (dct_params->flags == -1) return NULL;

    if (dct_params->flags % 2) {
        for (int i = dct_params->offset_end - 1; dct_params->offset_begin < i; i--) {
            core_compress_level1(dct_params->cores + (i - 1), dct_params->cores + i);
        }
        if (dct_params->flags > 1) {
            core_compress_level1(&(dct_params->dummy_left), dct_params->cores + dct_params->offset_begin);
        }
    } else {
        for (int i = dct_params->offset_end - 1; dct_params->offset_begin < i; i--) {
            core_compress_upper(dct_params->cores + (i - 1), dct_params->cores + i);
        }
        if (dct_params->flags > 1) {
            core_compress_upper(&(dct_params->dummy_left), dct_params->cores + dct_params->offset_begin);
        }
    }

    return NULL;
}

/**
 * @brief Executes one parallel DCT sweep over the cores array.
 *
 * This function divides the compression work starting at `dct_index`
 * into approximately equal chunks and distributes them among
 * `thread_number` threads. Each thread processes a contiguous
 * subrange of cores using the dct_worker() routine.
 *
 * Work partitioning:
 * - total = lps_ptr->size - dct_index
 * - chunk size is computed using ceiling division.
 * - Each thread receives a non-overlapping [begin, end) interval.
 *
 * Boundary handling:
 * - Threads whose range does not start at dct_index receive a
 *   copy of the left boundary core (`dummy_left`) to ensure
 *   correct cross-boundary compression.
 *
 * Memory management:
 * - Dynamically allocates thread handles and argument arrays.
 * - Joins all threads before returning.
 *
 * @param lps_ptr       Pointer to the LPS structure containing cores.
 * @param dct_index     Starting index for this DCT iteration.
 * @param thread_number Number of worker threads to spawn.
 *
 * @return 0 on success.
 * @return -1 on allocation failure.
 */
static int run_parallel_sweep(struct lps *lps_ptr, int dct_index, int thread_number) {

    const int total = (int)lps_ptr->size - dct_index;   // number of pairs (i, i+1)
    if (total <= 0) return 0;                               // nothing to do
    const int chunk = (total + thread_number - 1) / thread_number;

    pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * thread_number);
    dct_worker_args_t *args = (dct_worker_args_t *)malloc(sizeof(dct_worker_args_t) * thread_number);

    if (!threads || !args) {
        perror("LCP: couldn't allocate thread arguments\n");
        return -1;
    }

    for (int t = 0; t < thread_number; t++) {
        int begin = dct_index + t * chunk;
        int end = ((begin + chunk) <= lps_ptr->size ? (begin + chunk) : lps_ptr->size);

        args[t].cores = lps_ptr->cores;
        args[t].offset_begin = begin;
        args[t].offset_end = end;

        if ((int)lps_ptr->size <= begin || end <= begin) {
            args[t].flags = -1;
        } else if (begin != dct_index) {
            args[t].dummy_left = lps_ptr->cores[begin-1];
            args[t].flags = 2 + (lps_ptr->level == 1 ? 1 : 0);
        } else {
            args[t].flags = (lps_ptr->level == 1 ? 1 : 0);
        }
    }

    for (int t = 0; t < thread_number; t++) {
        pthread_create(&threads[t], NULL, dct_worker, &args[t]);
    }

    for (int t = 0; t < thread_number; t++) {
        pthread_join(threads[t], NULL);
    }

    free(threads);
    free(args);
    return 0;
}

/**
 * @brief Performs multi-iteration parallel Deterministic Coin Tossing (DCT).
 *
 * This function executes LCP_DCT_ITERATION_COUNT consecutive compression
 * sweeps over the cores stored in `lps_ptr`, using parallel execution
 * for each sweep.
 *
 * Each iteration reduces redundancy by compressing adjacent core pairs.
 * The process prepares the core sequence for subsequent parsing stages
 * in the LCP framework.
 *
 * Preconditions:
 * - At least LCP_DCT_ITERATION_COUNT + 1 cores must be available.
 *
 * @param lps_ptr       Pointer to the LPS structure.
 * @param thread_number Number of worker threads per sweep.
 *
 * @return 0 on success.
 * @return -1 if there are not enough cores for DCT.
 * @return Propagates non-zero errors from run_parallel_sweep().
 */
int lcp_dct_parallel(struct lps *lps_ptr, int thread_number) {
    // at least 2 cores are needed for compression
    if (lps_ptr->size < LCP_DCT_ITERATION_COUNT + 1) {
        return -1;
    }

    for (int dct_index = 0; dct_index < LCP_DCT_ITERATION_COUNT; dct_index++) {
        int rc = run_parallel_sweep(lps_ptr, dct_index, thread_number);
        if (rc != 0) {
            return rc;
        }
    }

    return 0;
}

/**
 * @brief Performs a single LPS deepening step using parallel DCT.
 *
 * This function advances the LPS structure by one level:
 *
 * 1. Executes parallel DCT compression.
 * 2. Parses the resulting cores using parse3() to identify new cores.
 * 3. Updates the core array and shrinks it via realloc().
 * 4. Increments the LPS level.
 *
 * If DCT cannot be performed (insufficient cores),
 * the structure is reset (size = 0) and the level is still incremented.
 *
 * @param lps_ptr       Pointer to the LPS structure to update.
 * @param thread_number Number of worker threads used during DCT.
 *
 * @return 1 if deepening produced a valid new level.
 * @return 0 if no further deepening was possible.
 */
int lps_deepen1_parallel(struct lps *lps_ptr, int thread_number) {

    // compress cores
    if (lcp_dct_parallel(lps_ptr, thread_number) < 0) {
        lps_ptr->size = 0;
        lps_ptr->level++;
        return 0;
    }

    // find new cores
#if !LCP_FIXED_CORE
    int new_size = parse3(lps_ptr->cores + LCP_DCT_ITERATION_COUNT, lps_ptr->cores + lps_ptr->size, lps_ptr->cores);
#else
    struct core *new_cores = (struct core *)malloc((lps_ptr->size - LCP_DCT_ITERATION_COUNT) * sizeof(struct core));
    if (new_cores == NULL) {
        return 0;
    }

    int new_size = parse3(lps_ptr->cores + LCP_DCT_ITERATION_COUNT, lps_ptr->cores + lps_ptr->size, new_cores);

    // replace the previous level; struct core owns no heap memory
    free(lps_ptr->cores);
    lps_ptr->cores = new_cores;
#endif
    lps_ptr->size = new_size;

    lps_ptr->level++;

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));

    return 1;
}

/**
 * @brief Deepens the LPS structure up to a target LCP level.
 *
 * Repeatedly invokes lps_deepen1_parallel() until either:
 * - The desired `lcp_level` is reached, or
 * - Further deepening is no longer possible.
 *
 * If the current level already satisfies the requested level,
 * no work is performed.
 *
 * @param lps_ptr       Pointer to the LPS structure.
 * @param lcp_level     Target LCP level to reach.
 * @param thread_number Number of worker threads used per deepening step.
 *
 * @return 1 if processing completed (even if stopped early).
 * @return 0 if no deepening was required.
 */
int lps_deepen_parallel(struct lps *lps_ptr, int lcp_level, int thread_number) {

    if (lcp_level <= lps_ptr->level)
        return 0;

    while (lps_ptr->level < lcp_level && lps_deepen1_parallel(lps_ptr, thread_number))
        ;

    return 1;
}

void print_lps(const struct lps *lps_ptr) {
    printf("Level: %d \n", lps_ptr->level);
    for(int i=0; i<lps_ptr->size; i++) {
        print_core(&(lps_ptr->cores[i]));
        printf(" ");
    }
}

int lps_eq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 0;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 0;
        }
    }

    return 1;
}

int lps_neq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 1;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 1;
        }
    }

    return 0;
}