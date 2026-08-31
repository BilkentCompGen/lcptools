/**
 * @file emit_var.h
 * @brief Variable-length core emission (CORE=var): one core per segment.
 *
 * Each helper forwards straight to init_coreN, so routing the parse functions 
 * through the seam cannot change what a CORE=var build produces.
 */

#ifndef LCP_EMIT_VAR_H
#define LCP_EMIT_VAR_H

/**
 * @brief Emits one forward-strand core spanning the whole segment.
 *
 * @param cores      Destination array.
 * @param core_index Current write index.
 * @param origin     Beginning of the whole sequence (unused in this mode).
 * @param offset     Global index shift (unused in this mode).
 * @param seg_begin  First character of the segment.
 * @param seg_len    Length of the segment, in characters.
 * @param start      Start position of the emitted core.
 * @param end        End position of the emitted core.
 * @return Updated write index.
 */
static inline int emit_core1(struct core *cores, int core_index, const char *origin, uint64_t offset, const char *seg_begin, lcp_pos seg_len, lcp_pos start, lcp_pos end) {

    (void)origin;
    (void)offset;

    init_core1(&(cores[core_index]), seg_begin, seg_len, start, end);

    return core_index + 1;
}

#if !LCP_ALPHABET_PROTEIN

/**
 * @brief Emits one reverse-complement core spanning the whole segment.
 *
 * `seg_begin` is the highest-address character of the segment; the segment
 * extends downwards from there.
 *
 * @param cores      Destination array.
 * @param core_index Current write index.
 * @param origin     End of the whole sequence (unused in this mode).
 * @param offset     Global index shift (unused in this mode).
 * @param seg_begin  Highest-address character of the segment.
 * @param seg_len    Length of the segment, in characters.
 * @param start      Start position of the emitted core.
 * @param end        End position of the emitted core.
 * @return Updated write index.
 */
static inline int emit_core2(struct core *cores, int core_index, const char *origin, uint64_t offset, const char *seg_begin, lcp_pos seg_len, lcp_pos start, lcp_pos end) {

    (void)origin;
    (void)offset;

    init_core2(&(cores[core_index]), seg_begin, seg_len, start, end);

    return core_index + 1;
}

#endif /* !LCP_ALPHABET_PROTEIN */

/**
 * @brief Emits one core spanning the whole segment of input cores.
 *
 * @param cores      Destination array.
 * @param core_index Current write index.
 * @param seg_begin  First core of the segment.
 * @param seg_len    Length of the segment, in cores.
 * @return Updated write index.
 */
static inline int emit_core3(struct core *cores, int core_index, struct core *seg_begin, lcp_pos seg_len) {

    init_core3(&(cores[core_index]), seg_begin, seg_len);

    return core_index + 1;
}

#endif
