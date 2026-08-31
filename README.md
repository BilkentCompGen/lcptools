# `LCP` (Locally Consistent Parsing) Algorithm Implementation <br>
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/BilkentCompGen/lcptools)
![GitHub last commit](https://img.shields.io/github/last-commit/BilkentCompGen/lcptools)
![GitHub](https://img.shields.io/github/license/BilkentCompGen/lcptools)

This repository contains an implementation of the Locally Consistent Parsing (LCP) algorithm, applied to strings using a specific binary alphabet encoding. The implementation is in C and is designed for efficient computation of LCP on large datasets.

For additional details about the data structure, please refer to [this document](https://github.com/BilkentCompGen/lcptools/blob/main/docs.md).

If you use LCP (or *lcptools*) in your work, please cite:
- LCPan: efficient variation graph construction using locally consistent parsing. Akmuhammet Ashyralyyev, Zülal Bingöl, Begüm Filiz Öz, Salem Malikic, Uzi Vishkin, S. Cenk Sahinalp, Can Alkan. [Genome Biol (2026)](https://doi.org/10.1186/s13059-026-04088-w).

## Features

- **Efficient LCP Computation:** Implemented in C for efficient and scalable computation on large datasets.
- **High Accuracy:** Achieves highly accurate comparisons by leveraging the unique LCP approach.
- **Designed for Genomics:** Specifically caters to the needs of genomic researchers and bioinformaticians.

## Installation

You can install lcptools either system-wide (requires sudo privileges) or in a user-specific directory (no sudo required).

### System-wide Installation

To install lcptools system-wide, you need sudo privileges. This will install the library in `/usr/local`.

1. **Install the repository:**
    ```sh
    git clone https://github.com/BilkentCompGen/lcptools.git
    cd lcptools
    git checkout biolcp

    # Install the library
    sudo make install
    ```

2. **Uninstall the library (if needed):**
    ```sh
    sudo make uninstall
    ```

### User-specific Installation

To install lcptools in your home directory (or another custom directory), you don't need sudo privileges.

1. **Install the repository:**
    ```sh
    git clone https://github.com/BilkentCompGen/lcptools.git
    cd lcptools
    git checkout biolcp

    # Install the library to a custom directory (e.g., `~/.local`):**
    make install PREFIX=$(HOME)/.local
    ```

2. **Uninstall the library from the custom directory (if needed):**
    ```sh
    make uninstall PREFIX=$(HOME)/.local
    ```
#### Build Configuration

This project uses build-time configuration to generate type definitions for core data structures. Type widths are determined at compile time via the Makefile, allowing you to optimize struct sizes for your specific use case.

##### Configurable Parameters

| Variable   | Values           | Default | Effect                                                                           |
| ---------- | ---------------- | ------- | -------------------------------------------------------------------------------- |
| `ALPHABET` | `dna`, `protein` | `dna`   | Symbol table and bits per symbol. See [Character Encoding](#character-encoding). |
| `CORE`     | `var`, `fixed`   | `var`   | How RINT and SSEQ segments become cores. See [Core length](#core-length).        |
| `LABEL`    | `0`, `32`, `64`  | `32`    | Width of `lcp_label` in bits. `0` removes the field from `struct core`.          |
| `POS`      | `0`, `32`, `64`  | `32`    | Width of `lcp_pos` in bits. `0` removes `start`/`end` from `struct core`.        |
| `DCT`      | integer          | `1`     | Default DCT iteration count.                                                     |
| `PREFIX`   | path             | `/usr/local` | Install location.                                                           |
| `SUFFIX`   | string           | empty   | Appended to the library name, for installing several variants side by side.      |

Two values are **derived** from the above and are not settable directly, because setting them inconsistently would corrupt the output:

- `LCP_SYMBOL_BITS` - 2 for DNA, 8 for protein.
- `LCP_CONSTANT_FACTOR` - the allocation divisor, 1.5 for `CORE=var` and 1 for `CORE=fixed`. Fixed-length emission produces up to one core per position, so a larger divisor would under-allocate.

##### Building with Custom Configuration

Pass configuration variables to `make`:

```sh
make install LABEL=64 POS=64 DCT=1
make install ALPHABET=protein
make install CORE=fixed
```

`make install` echoes the variant it built, for example `dna-var-l32-p32-d1`, so you can confirm what an installed library actually is.

`make test` runs the test suite against the installed library under the same configuration; pass the same variables you installed with:

```sh
make install ALPHABET=protein PREFIX=$(HOME)/.local
make test    ALPHABET=protein PREFIX=$(HOME)/.local
```

##### Installing several variants

Every configuration installs the same `liblcptools` name and the same `config.h`, so installing a second configuration over the first leaves one library paired with the other's header - a silently misread `struct core` rather than a link error. 
Use `SUFFIX` to keep them apart:

```sh
make install ALPHABET=protein SUFFIX=-protein   # liblcptools-protein.a
```

**Note:** the build generates `config.h` from `config.h.in` by substituting your values, and the headers installed alongside the library must be the ones it was built with. 
If you change parameters, re-run `make install`; mismatched binaries and headers cause undefined behavior.

## Usage

To compile your program with your program, you need to specify the include and library paths based on your installation method.

### Compile with System-wide installed library

If you want to link static library, please use as follows:

```sh
g++ your_program.cpp -static -llcptools -o your_program
```

If you want to link dynamic library, please use as follows:

```sh
g++ your_program.cpp -llcptools -o your_program
```

### Compile with User-specific installed library

If you want to link static library, please use as follows:

```sh
g++ your_program.cpp -static -I$(HOME)/.local/include -L$(HOME)/.local/lib -llcptools -o your_program
```

If you want to link dynamic library, please use as follows:

```sh
g++ your_program.cpp -I$(HOME)/.local/include -L$(HOME)/.local/lib -llcptools -Wl,-rpath,$(HOME)/.local/lib -o your_program
```

**Note**: Make sure that paths are correct.

## Character Encoding

The symbol table is chosen at build time with `ALPHABET`. This default encoding is used unless a custom encoding is provided.

### `ALPHABET=dna` (default)

Four symbols at 2 bits each, plus a reverse-complement table:

| Character | Encoding | Reverse complement |
| --------- | -------- | ------------------ |
| A, a      | 00       | 11                 |
| C, c      | 01       | 10                 |
| G, g      | 10       | 01                 |
| T, t      | 11       | 00                 |

### `ALPHABET=protein`

The 20 standard amino acids plus the ambiguity and non-standard codes B (Asx),
Z (Glx), X (any), U (selenocysteine) and O (pyrrolysine) - 25 symbols at 8 bits
each, numbered alphabetically over the whole set, upper and lower case alike:

| Character | Encoding | | Character | Encoding | | Character | Encoding |
| --------- | -------- |-| --------- | -------- |-| --------- | -------- |
| A         | 0        | | K         | 9        | | T         | 18       |
| B         | 1        | | L         | 10       | | U         | 19       |
| C         | 2        | | M         | 11       | | V         | 20       |
| D         | 3        | | N         | 12       | | W         | 21       |
| E         | 4        | | O         | 13       | | X         | 22       |
| F         | 5        | | P         | 14       | | Y         | 23       |
| G         | 6        | | Q         | 15       | | Z         | 24       |
| H         | 7        | | R         | 16       | |           |          |
| I         | 8        | | S         | 17       | |           |          |

`J` and every other character are invalid, exactly as non-ACGT characters are in a DNA build.

**The ordering is part of the algorithm, not a display detail.** LMIN and LMAX are decided by comparing these values, so a different ordering produces a different parse. 
Alphabetical is the default precisely because it implies nothing; if the parse should follow a biochemical property instead, change `PROTEIN_SYMBOLS` in `encoding.c` or supply a custom encoding file.

Proteins have no reverse complement, so a protein build does not provide `init_lps2`, `init_core2` or `rc_alphabet` at all. 
Referring to them is a compile error rather than a silently meaningless result.

### Level 1 core layout

A level 1 core packs its three significant symbols and a run length into one 64-bit word, which is why the symbol width matters beyond the encoding table:

| Bits            | Contents                    |
| --------------- | --------------------------- |
| 63              | level 1 sentinel            |
| 62 ... 3*SB     | run length (`distance - 2`) |
| 3*SB-1 ... 2*SB | first symbol                |
| 2*SB-1 ... SB   | second-to-last symbol       |
| SB-1 ... 0      | last symbol                 |

`SB` is `LCP_SYMBOL_BITS`: 2 for DNA, giving a 57-bit length, and 8 for protein, giving 39 bits. 
Upper-level cores do not use this layout; they hold compressed bit strings.

### Core length

`CORE` selects how a RINT or SSEQ segment becomes cores. LMIN and LMAX cores span exactly three symbols in both modes.

- `CORE=var` (default) - one core spanning the whole segment, so cores vary in length.
- `CORE=fixed` - every length-3 sliding window of the segment, so every core spans exactly three symbols.

### Initialization

To initialize the encodings, use the following function call at the beginning of your program.

```cpp
LCP_INIT();
```

In the above code, defaults the verbose to `0`.

A integer parameter `verbose` can be provided, which, when set to `1`, prints a summary of the encoding:

```cpp
LCP_INIT2(verbose);
```

To display the encoding summary separately, use:

```cpp
LCP_SUMMARY();
```

### Custom encodings

`LCP_INIT_FILE(path, verbose)` replaces the built-in table with one read from a file. 
The column count depends on the alphabet, since only DNA has a complement:

```
# ALPHABET=dna       <symbol> <encoding> <reverse complement>
A 0 3
C 1 2

# ALPHABET=protein   <symbol> <encoding>
M 0
K 1
```

Symbols absent from the file become invalid. An encoding wider than `LCP_SYMBOL_BITS` is rejected, because it would overflow into the neighbouring symbol's field in a level 1 core.

## Usage Example

Below is an example demonstrating the usage of the LCP algorithm implementation:

```c
#include "lps.h"

int main() {

    // Initialize alphabet coefficients
    LCP_INIT();

    // Example string
    const char *str = "GGGACCTGGTGACCCCAGCCCACGACAGCCAAGCGCCAGCTGAGCTCAGGTGTGAGGAGATCACAGTCCT";

    // Create LCP string object
    struct lps lcp_str;
    init_lps(&lcp_str, str, strlen(str));

    // Deepen the LCP analysis
    int isSuccess = lps_deepen(&lcp_str, 2);

    // Output LCP string
    print_lps(&lcp_str);

    // Clean up to prevent memory leaks
    free_lps(&lcp_str);

    return 0;
}
```

## LCP Algorithm Description

The LCP algorithm operates as follows:

### Constructor:

Processes the input string and identifies cores that adhere to specific rules:

1. (LMIN) The subsequent characters should not be the same, and the middle character is local minima.

	Ex: $w = xyz$ where $x \neq y$ and $y \neq z$, and $x \gt y$ and $y \lt z$

3. (LMAX) The subsequent characters should not be the same, and the middle character local maxima, and its neighbors are not local minima. 

	Ex: $w = sxyzt$ where $s \neq x$ and $x \neq y$ and $y \neq z$ and $z \neq t$, and $s \leq x$ and $x \lt y$ and $y \gt z$ and $z \geq t$.

4. (RINT) The characters, except the front and back, are the same.

	Ex: $w=xy^iz$ where $i > 1$ and $|w| \gt 3$, , and $x\neq y$, $y\neq z$.

5. (SSEQ) The subsequent characters are either strictly increasing or decreasing with respect to the lexicographic order, and only the first and last characters are part of either a LMIN, LMAX, or a RINT. 
    
    Ex: $w = xyza_1 . . . a_nklm$, where $n \geq 1$ and $xyz$ and $klm$ are identified as cores, and $z \lt a_1 \lt \dots \lt a_n \lt k$ or $z \gt a_1 \gt \dots \gt a_n \gt k$.

### Deterministic Coin Tossing:

The dct function in the LCP algorithm is crucial for processing binary sequences. It starts by pinpointing the initial point of difference between two binary strings, beginning from the right-end. The function then assesses the difference based on the position and value of the divergent bit. This detail is transformed into a new binary sequence, which establishes the foundation of a newly generated 'core'. This core is a clear representation of the differences between the original sequences, integral to the algorithm's deepening process. Essentially, the dct function effectively consolidates and encapsulates the information, ensuring efficient further analysis within the LCP framework.

Ex: 11101**0**00 vs 00010**1**00 -> **10**0 as the position is **2 (10)** and the bit is **0**. Position index start from 0.

### Deepen Function:

The deepen function in the LCP algorithm primarily focuses on the compression of 'cores' alongside their left neighbors. The purpose of this repeated compression (dct) is to manage the length of the cores, preventing them from becoming large. After a compression, the LCP algorithm is re-applied. This re-application aims to identify new cores within the compressed data. In this context, each compressed core is treated as a discrete value, represented in binary form. This representation facilitates efficient processing and analysis within the algorithm.

This function iteratively compresses and processes cores to find new cores in compliance with the rules stated above.

## Default Variables

A plain `make install` builds `ALPHABET=dna CORE=var LABEL=32 POS=32 DCT=1`: the DNA alphabet, one core per segment, a 32-bit `lcp_label`, 32-bit `lcp_pos` begin and end offsets, and one DCT iteration per deepening. 
See [Build Configuration](#build-configuration) to change any of them.

## Platform Notes

Linux and macOS are both supported. 
The build probes `uname -s` and adjusts: macOS produces `liblcptools.dylib` with its absolute install path recorded, so programs linked against it resolve the library after `make install` without setting `DYLD_LIBRARY_PATH`.
 Linux produces `liblcptools.so` with a soname, and relies on the `-rpath` shown in the [Usage](#usage) examples.
