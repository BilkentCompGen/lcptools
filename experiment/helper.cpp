/**
 * @file    helper.cpp
 * @brief   Statistical Analysis of Genomic Sequence Alignments
 *
 * This program performs statistical analysis on genomic sequence alignment data.
 * It includes functionality for calculating the mean and standard deviation of distances
 * and lengths in the alignments. The analysis is done at multiple levels, defined by
 * LCP_LEVEL, and caters to a range of distance lengths, specified by DISTANCE_LENGTH.
 * The results are printed to a file and include detailed statistics at each level of
 * analysis.
 *
 */

#include <cmath>
#include <map>
#include <vector>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>

#define DISTANCE_LENGTH 10000
#define LCP_LEVEL 8
#define MAXIMUM_FREQ_THRESHOLD 256
#define MAX_DISTANCE_THRESHOLD 30

typedef unsigned long stats_type;

/**
 * @brief Formats an integer with thousands separators for better readability.
 *
 * This function converts an integer to a string representation that includes thousands
 * separators, based on the current locale settings. The resulting string makes large
 * numbers more readable by grouping digits into thousands using separators (e.g., commas
 * in English locales).
 *
 * @param value The integer to be formatted.
 * @return A std::string containing the formatted integer with thousands separators.
 *
 * Example:
 *     format_int(1234567) returns "1,234,567" (in an English locale).
 */
std::string format_int(int value) {
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << value;
    return ss.str();
};

/**
 * @brief Formats a double value with two decimal places and thousands separators.
 *
 * This function converts a double to a string representation that includes thousands
 * separators and ensures that the number is formatted to two decimal places. The
 * thousands separators are applied based on the current locale settings, making the
 * number easier to read, especially for large values.
 *
 * @param value The double value to be formatted.
 * @return A std::string containing the formatted double with thousands separators
 *         and two decimal places.
 *
 * Example:
 *     format_double(1234567.89123) returns "1,234,567.89" (in an English locale).
 */
std::string format_double(double value, size_t precision = 2) {
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
};

/**
 * @brief Calculates the mean of integers.
 *
 * This function calculates the mean of integers in a given range. It supports
 * processing of standard integers within a predefined range (DISTANCE_LENGTH)
 * and additional larger integers provided as a vector.
 * 
 * @param numbers Array of integers within the predefined range.
 * @param numbersXL Vector of integers outside the predefined range.
 * @return The mean of all provided integers.
 */
double mean(int (&numbers)[DISTANCE_LENGTH], std::vector<int> numbersXL = {}) {
    double sum = 0;
    double count = 0;
    for (size_t i = 0; i < DISTANCE_LENGTH; i++) {
        sum += (i * numbers[i]);
        count += numbers[i];
    }
    for (size_t i = 0; i < numbersXL.size(); i++) {
        sum += numbersXL[i];
    }
    count += numbersXL.size();
    return sum / count;
};

/**
 * @brief Calculates the standard deviation of integers.
 *
 * This function calculates the standard deviation of integers in a given range.
 * It leverages the mean function for its calculation and handles both standard
 * and larger integers.
 *
 * @param numbers Array of integers within the predefined range.
 * @param numbersXL Vector of integers outside the predefined range.
 * @return The standard deviation of all provided integers.
 */
double stdev(int (&numbers)[DISTANCE_LENGTH], std::vector<int> numbersXL = {}) {
    double mean_value = mean(numbers, numbersXL);
    double count = 0;
    for (size_t i = 0; i < DISTANCE_LENGTH; i++) {
        count += numbers[i];
    }
    count += numbersXL.size();
    double variance = 0;
    for (size_t i = 0; i < DISTANCE_LENGTH; i++) {
        variance += ((mean_value - i) * (mean_value - i) * numbers[i]);
    }
    for (size_t i = 0; i < numbersXL.size(); i++) {
        variance += ((mean_value - numbersXL[i]) * (mean_value - numbersXL[i]));
    }
    return sqrt(variance / count);
};

/**
 * @brief Calculates the mean of distances.
 *
 * This function computes the weighted mean of an array of distances. Each index
 * in the array corresponds to a distance, and the value at that index represents
 * how many times that distance occurs.
 *
 * @param distances Pointer to an array of integers where each element represents
 *                  the count of occurrences of the corresponding index distance.
 * @param n         The size of the distances array.
 * @return The weighted mean of the distances.
 */
double mean(int* distances, int n) {
    double numOfDistances = 0;
    double totalDistance = 0;
    for (int i = 0; i < n; i++)
    {
        totalDistance += i * distances[i];
        numOfDistances += distances[i];
    }
    return totalDistance / numOfDistances;
};

/**
 * @brief Calculates the standard deviation of distances.
 *
 * This function computes the standard deviation of distances based on their
 * frequency in the array. The mean value is provided as a parameter for use in the
 * calculation.
 *
 * @param distances Pointer to an array of integers representing the counts of each
 *                  distance.
 * @param n         The size of the distances array.
 * @param mean      The precomputed mean of the distances.
 * @return The standard deviation of the distances.
 */
double stdev(int* distances, int n, double mean) {
    double sum = 0;
    double N = 0;
    for (int i = 0; i < n; i++)
    {
        N += distances[i];
        sum += pow(((i) - mean), 2) * distances[i];
    }
    return sqrt(sum / N);
};

void reverse_complement(std::string &str) {
    std::reverse(str.begin(), str.end());
    for (size_t i = 0; i < str.length(); i++) {
        if (str[i] == 'A') {
            str[i] = 'T';
        } else if (str[i] == 'T') {
            str[i] = 'A';
        } else if (str[i] == 'G') {
            str[i] = 'C';
        } else if (str[i] == 'C') {
            str[i] = 'G';
        }
    }
};

template<typename kmer_type, typename index_type>
void alignment_pairwise_stats(const std::map<kmer_type, std::vector<index_type>> &mapGTRead, 
                              const std::map<kmer_type, std::vector<index_type>> &mapSimRead, 
                              stats_type &true_positive, stats_type &false_positive, stats_type &false_negative, stats_type &total) {
    // Process ground truth read to find matches and mismatches
    for (const std::pair<const kmer_type, std::vector<index_type>> &gt_pair : mapGTRead) {
        
        kmer_type id = gt_pair.first;
        const std::vector<index_type> &gt_indices = gt_pair.second;

        bool isInSim = mapSimRead.find(id) != mapSimRead.end();
        int currMatch = 0;

        if (isInSim) {
            const std::vector<index_type> &sim_indices = mapSimRead.at(id);
            size_t index1 = 0; // gt index
            size_t index2 = 0; // sim index
            while(index1 < gt_indices.size() && index2 < sim_indices.size()) {
                if ((sim_indices[index2] >= gt_indices[index1] && sim_indices[index2] - gt_indices[index1] < MAX_DISTANCE_THRESHOLD) ||
                    (sim_indices[index2] < gt_indices[index1] && gt_indices[index1] - sim_indices[index2] < MAX_DISTANCE_THRESHOLD)) {
                    index1++;
                    index2++;
                    currMatch++;
                } else if (sim_indices[index2] > gt_indices[index1]) {
                    index1++;
                } else {
                    index2++;
                }
            }
            true_positive += currMatch;
            false_positive += sim_indices.size() - currMatch;
            total += sim_indices.size() - currMatch;
        } 
        false_negative += gt_indices.size() - currMatch;
        total += gt_indices.size();
    }
    // Consider sketches that are in simulated read but not in original sequence
    for (const std::pair<const kmer_type, std::vector<index_type>> &sim_pair : mapSimRead) {
        kmer_type id = sim_pair.first;
        if (mapGTRead.find(id) == mapGTRead.end()) {
            false_positive += sim_pair.second.size(); // false positive
            total += sim_pair.second.size();
        }
    }
};

template<typename kmer_type, typename index_type>
void alignment_global_stats(const std::map<kmer_type, std::vector<index_type>> &mapReference, 
                            const std::map<kmer_type, std::vector<index_type>> &mapGTRead, 
                            const std::map<kmer_type, std::vector<index_type>> &mapSimRead, 
                            stats_type &true_positive, stats_type &false_positive, stats_type &false_negative) {
    
    for (const std::pair<const kmer_type, std::vector<index_type>> &sim_pair : mapSimRead) {
        kmer_type id = sim_pair.first;
        
        bool isInChr = mapReference.find(id) != mapReference.end();
        if (!isInChr) {
            continue;
        }

        const std::vector<index_type> &ref_indices = mapReference.at(id); // we are sure that it exists somewhere
        const std::vector<index_type> &sim_indices = sim_pair.second;

        if (MAXIMUM_FREQ_THRESHOLD <= ref_indices.size()) {
            continue;
        }

        bool isInGTRead = mapGTRead.find(id) != mapGTRead.end();
        if (!isInGTRead) {
            false_positive += ref_indices.size() * sim_indices.size(); // as each sketch will hit to incorrect spot
            continue;
        }

        const std::vector<index_type> &gt_indices = mapGTRead.at(id); // we are sure that sketch exists inside the read now
        
        // we process each sketch to find the true positive and false positive counts
        for (const index_type &sim_index : sim_indices) {   
            bool is_tp = false;

            for (const index_type &gt_index : gt_indices) {
                if (gt_index+MAX_DISTANCE_THRESHOLD < sim_index) { // not came to the correct place yet
                    continue;
                } 
                if (sim_index+MAX_DISTANCE_THRESHOLD < gt_index) { // exceeds the read's original locations
                    break;
                }
                is_tp = true;
                break;
            }

            if (!is_tp) {
                false_positive += ref_indices.size(); // as there no match and we have ref_count number of incorrect hit
            } else {
                true_positive++; // there is single correct match
                false_positive += ref_indices.size()-1; // the rest is incorrect match
            }
        }
    }
    // Consider ids that are in read but not in sequence -> FN
    for (const std::pair<const kmer_type, std::vector<index_type>> &gt_pair : mapGTRead) {
        kmer_type id = gt_pair.first;
        const std::vector<index_type> &ref_indices = mapReference.at(id);

        // check if the id was filtered out
        if (MAXIMUM_FREQ_THRESHOLD <= ref_indices.size()) {
            continue;
        }

        bool isInSimRead = mapSimRead.find(id) != mapSimRead.end();
        const std::vector<index_type> &gt_indices = gt_pair.second;

        if (!isInSimRead) {
            false_negative += gt_indices.size();
            continue;
        }

        const std::vector<index_type> &sim_indices = mapSimRead.at(id); // we know that the id exists in both

        // now we need to find how many of them are true positive, so we can exclude them in false negative
        for (const index_type &gt_index : gt_indices) {
            for (const index_type &sim_index : sim_indices) {
                if (sim_index+MAX_DISTANCE_THRESHOLD < gt_index) {
                    continue;
                } 
                else if (gt_index+MAX_DISTANCE_THRESHOLD < sim_index) { // exceeds the read's original locations, if not foun until now, then FN
                    false_negative++;
                    break;
                }
                // indeces matches, then it is false_positive and should have been counted earlier as TP
                break;
            }
        }
    }
};