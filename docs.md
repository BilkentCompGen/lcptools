# `lps` Data Structure Documentation

The `lps` data structure is designed for efficient parsing and representation of strings with advanced features such as reverse complement transformations, offset-based indexing, and binary file serialization/deserialization. 

This document provides detailed instructions on initializing, using, and managing the `lps` data structure.

---

## Overview

The `lps` structure is composed of two main components:

1. **`lps` structure**: Represents the main object with attributes like `level`, `size`, and an array of `core` objects.
2. **`core` structure**: Represents the individual elements stored in the `lps`, with attributes for bit representation, labels, and positional information.

### Struct Definitions
```c
struct lps {
    int level;
    int size;
    struct core *cores;
};

struct core {
   lcp_pos start;       // Start index in the string   (omitted when POS=0)
   lcp_pos end;         // End index in the string     (omitted when POS=0)
   uint64_t bit_rep;    // Bit representation of the core
   ubit_size bit_size;  // Number of bits held in bit_rep
   lcp_label label;     // Unique label for the core   (omitted when LABEL=0)
};
```

`start`/`end` and `label` are present only when `POS` and `LABEL` are non-zero; see the build configuration section of the README. 
Code that touches them should be guarded with `#if LCP_POS_BITS != 0` / `#if LCP_LABEL_BITS != 0` so it keeps compiling across configurations.

---

Besides the original implementation, this version is highly optimized for bioinformatics tools.

LCP level-1 cores are encoded in a single 64-bit word, most significant bit first:

```txt
<sentinel><middle-char-count><left-char><middle-char><right-char>
     1 bit        rest           SB bits    SB bits    SB bits
```

`SB` is `LCP_SYMBOL_BITS`, set by the `ALPHABET` build variable: **2** for DNA, leaving 57 bits for the count, and **8** for protein, leaving 39. 
The sentinel at bit 63 marks a core as level 1.

Upper-level cores are encoded as concatenated forms of compressed lower-level cores.
The functionalities are implemented to handle both types internally, but a symbol's encoding must fit in `SB` bits - a wider one would overflow into the neighbouring character's field.

## Functions

### Initialization

#### `init_lps`
Initializes an `lps` object using a given string.

**Parameters**:
- `struct lps *lps_ptr`: Pointer to the `lps` object to initialize.
- `const char *str`: Input string to be parsed.
- `int len`: Length of the input string.

**Usage**:
```c
struct lps my_lps;
init_lps(&my_lps, "ACGTACGT", 8);
```

---

#### `init_lps_offset`
Initializes an `lps` object with offset-based indexing. This function adds the `offset` value to each cores' indices.

**Parameters**:
- `struct lps *lps_ptr`: Pointer to the `lps` object to initialize.
- `const char *str`: Input string to be parsed.
- `int len`: Length of the input string.
- `uint64_t offset`: Offset value for indexing.

**Usage**:
```c
struct lps my_lps;
init_lps_offset(&my_lps, "ACGTACGT", 8, 100);
```

---

#### `init_lps2`
Initializes an `lps` object and includes a reverse complement transformation.

Available only when `ALPHABET=dna`. 
A protein build does not declare it, since proteins have no reverse complement.

**Parameters**:
- `struct lps *lps_ptr`: Pointer to the `lps` object to initialize.
- `const char *str`: Input string to be parsed.
- `int len`: Length of the input string.

**Usage**:
```c
struct lps my_lps;
init_lps2(&my_lps, "ACGTACGT", 8);
```

---

#### `init_lps3`
Initializes an `lps` object by reading from a binary file.

**Parameters**:
- `struct lps *lps_ptr`: Pointer to the `lps` object to initialize.
- `FILE *in`: File pointer to the binary file containing serialized `lps` data.

**Usage**:
```c
struct lps my_lps;
FILE *input_file = fopen("lps_data.bin", "rb");
init_lps3(&my_lps, input_file);
fclose(input_file);
```

---

#### `init_lps4`
Initializes an `lps` object using divide and conquer approach.

**Parameters**:
- `struct lps *lps_ptr`: Pointer to the `lps` object to initialize.
- `const char *str`: Input string to be parsed.
- `int len`: Length of the input string.
- `int lcp_level`: Level to deepen each chunk to before merging.
- `int chunk_size`: Size of the chunks to be processed.

Available only when `POS` is non-zero, since merging chunks relies on core positions.

**Usage**:
```c
struct lps my_lps;
init_lps4(&my_lps, "ACGTACGT", 8, 2, 100000);
```

---

### Memory Management

#### `free_lps`
Frees dynamically allocated memory associated with an `lps` object.

**Parameters**:
- `struct lps *lps_ptr`: Pointer to the `lps` object to free.

**Usage**:
```c
free_lps(&my_lps);
```

---

### File Operations

#### `write_lps`
Serializes and writes an `lps` object to a binary file.

**Parameters**:
- `struct lps *lps_ptr`: Pointer to the `lps` object to write.
- `FILE *out`: File pointer to the binary file for writing.

**Usage**:
```c
FILE *output_file = fopen("lps_data.bin", "wb");
write_lps(&my_lps, output_file);
fclose(output_file);
```

---

# Alphabet Encoding

This section provides functions to manage the encoding of standard DNA bases (A, C, G, T) and their complements, used in the Locally Consistent Parsing (LCP) tool. 
Please note that any custom alphabet encoding can be provided to the program.

## Functions

### ``LCP_SUMMARY``
   - **Description**: Displays the summary of the alphabet encoding, including coefficients and dictionary bit size. This function helps you understand the encoding setup and verify that the parameters are correctly configured.
   - **Usage**: Simply call `LCP_SUMMARY()` to print the encoding details to the console.
   
```c
void LCP_SUMMARY();
```

---

### ``LCP_INIT``
   - **Description**: Initializes the encoding coefficients for the alphabet this build was configured with: the standard DNA bases (A, C, G, T) and their complements at `ALPHABET=dna`, or the 25 amino-acid codes at `ALPHABET=protein`. The function sets default values for the coefficients and dictionary bit size.
   - **Usage**: Call `LCP_INIT()` to initialize the encoding with the default values.
   
```c
void LCP_INIT();
```

---

### ``LCP_INIT2``
   - **Description**: Same as `LCP_INIT`, for whichever alphabet the build was configured with. Sets default values for coefficients and dictionary bit size. The verbosity of the output can be controlled by the `verbose` parameter:
     - If `verbose` is 0, no encoding summary will be printed.
     - If `verbose` is 1, the encoding summary is printed after initialization.
   - **Usage**: Call `LCP_INIT2(0)` or `LCP_INIT2(1)` based on whether you want the encoding summary printed.
   
```c
void LCP_INIT2(int verbose);
```

---

### ``LCP_INIT_FILE``
   - **Description**: Initializes the encoding coefficients by reading them from a file. The column count depends on the alphabet: three at `ALPHABET=dna` (character, encoding, complement encoding) and two at `ALPHABET=protein`, which has no complement. Characters absent from the file become invalid. After initializing the encoding coefficients, the function prints the encoding summary if `verbose` is set to 1.
   - **Parameters**:
     - `filename`: The path to the file containing the character encodings, in the format shown below.
     - `verbose`: If true (1), the encoding summary will be printed after initialization. If false (0), no summary will be printed.
   - **Returns**: 0 upon successful initialization, -1 if the file cannot be opened.
   - **Exits**: calls `exit(EXIT_FAILURE)` if any encoding needs more than `LCP_SYMBOL_BITS` bits, since it would not fit a level-1 core's symbol field.
   - **Usage**: Call `LCP_INIT_FILE("path/to/encoding_file.txt", 1)` to initialize the encoding from a file and print the summary.
   
```c
int LCP_INIT_FILE(const char *filename, int verbose);
```

## File Format for LCP_INIT_FILE

At `ALPHABET=dna` the file has three columns:

```
A 0 1
C 1 0
G 2 3
T 3 2
```

Where:
- The first column represents the DNA base (A, C, G, T).
- The second column represents the encoding value.
- The third column represents the complement encoding value.

At `ALPHABET=protein` there is no complement, so the file has two columns:

```
M 0
K 1
T 2
A 3
```

---

## Example Workflow

```c
#include <stdio.h>
#include "lps.h"

int main() {
    LCP_INIT();

    struct lps my_lps;

    // Initialize from a string
    init_lps(&my_lps, "ACGTACGT", 8);

    // Serialize to a file
    FILE *output_file = fopen("lps_data.bin", "wb");
    write_lps(&my_lps, output_file);
    fclose(output_file);

    // Free memory
    free_lps(&my_lps);

    // Deserialize from a file
    FILE *input_file = fopen("lps_data.bin", "rb");
    init_lps3(&my_lps, input_file);
    fclose(input_file);

    // Free memory again
    free_lps(&my_lps);

    return 0;
}
```

---

## Notes
- The encoding coefficients define how each DNA base and its complement are represented internally. These values are essential for performing efficient parsing and comparison of DNA sequences.
- Ensure that the encoding file is properly formatted to avoid errors during initialization.
- Always ensure memory is freed after use to avoid leaks.
- Use valid and open binary files for serialization/deserialization functions.
- Offset and complement functions provide advanced features for genomic data analysis.

---
