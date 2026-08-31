/**
 * @file emit_fixed.h
 * @brief Fixed-length core emission (CORE=fixed): every length-3 sliding
 *        window of a segment.
 *
 * Positions are derived per window from `origin` and `offset`, so the `start`
 * and `end` a variable-length core would have carried are unused here.
 *
 * This mode emits up to one core per position rather than one per segment,
 * which is why LCP_CONSTANT_FACTOR drops to 1; see config.h.in.
 */

#ifndef LCP_EMIT_FIXED_H
#define LCP_EMIT_FIXED_H

/**
 * @brief Emits every length-3 window of a forward-strand segment.
 *
 * @param cores      Destination array.
 * @param core_index Current write index.
 * @param origin     Beginning of the whole sequence (index origin).
 * @param offset     Global index shift.
 * @param seg_begin  First character of the segment.
 * @param seg_len    Length of the segment, in characters.
 * @param start      Unused in this mode.
 * @param end        Unused in this mode.
 * @return Updated write index.
 */
static inline int emit_core1(struct core *cores, int core_index, const char *origin, uint64_t offset, const char *seg_begin, lcp_pos seg_len, lcp_pos start, lcp_pos end) {

    (void)start;
    (void)end;

    for (lcp_pos i = 0; i + 3 <= seg_len; i++) {
        const char *p = seg_begin + i;
        init_core1(&(cores[core_index]), p, 3, (uint64_t)(p - origin) + offset, (uint64_t)(p + 3 - origin) + offset);
        core_index++;
    }

    return core_index;
}

#if !LCP_ALPHABET_PROTEIN

/**
 * @brief Emits every length-3 window of a reverse-complement segment.
 *
 * `seg_begin` is the highest-address character of the segment; the segment
 * extends downwards from there.
 *
 * @param cores      Destination array.
 * @param core_index Current write index.
 * @param origin     End of the whole sequence (index origin).
 * @param offset     Global index shift.
 * @param seg_begin  Highest-address character of the segment.
 * @param seg_len    Length of the segment, in characters.
 * @param start      Unused in this mode.
 * @param end        Unused in this mode.
 * @return Updated write index.
 */
static inline int emit_core2(struct core *cores, int core_index, const char *origin, uint64_t offset, const char *seg_begin, lcp_pos seg_len, lcp_pos start, lcp_pos end) {

    (void)start;
    (void)end;

    for (lcp_pos i = 0; i + 3 <= seg_len; i++) {
        const char *p = seg_begin - i;
        init_core2(&(cores[core_index]), p, 3, (uint64_t)(origin - p - 1) + offset, (uint64_t)(origin - p + 2) + offset);
        core_index++;
    }

    return core_index;
}

#endif /* !LCP_ALPHABET_PROTEIN */

/**
 * @brief Emits every length-3 window of a segment of cores.
 *
 * @param cores      Destination array.
 * @param core_index Current write index.
 * @param seg_begin  First core of the segment.
 * @param seg_len    Length of the segment, in cores.
 * @return Updated write index.
 */
static inline int emit_core3(struct core *cores, int core_index, struct core *seg_begin, lcp_pos seg_len) {

    for (lcp_pos i = 0; i + 3 <= seg_len; i++) {
        init_core3(&(cores[core_index]), seg_begin + i, 3);
        core_index++;
    }

    return core_index;
}

#endif
