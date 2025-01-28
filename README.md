# `LCP` (Locally Consistent Parsing) Algorithm Implementation <br>
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/BilkentCompGen/lcptools)
![GitHub last commit](https://img.shields.io/github/last-commit/BilkentCompGen/lcptools)
![GitHub](https://img.shields.io/github/license/BilkentCompGen/lcptools)

This repository contains an implementation of the Locally Consistent Parsing (LCP) algorithm, applied to strings using a specific binary alphabet encoding. The implementation is in C and is designed for efficient computation of LCP on large datasets.

For additional details about the data structure, please refer to [this document](https://github.com/BilkentCompGen/lcptools/blob/main/docs.md).

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

    # Install the library to a custom directory (e.g., `~/.local`):**
    make install PREFIX=$(HOME)/.local
    ```

2. **Uninstall the library from the custom directory (if needed):**
    ```sh
    make uninstall PREFIX=$(HOME)/.local
    ```

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
g++ your_program.cpp -static -I$(HOME)/.local/include/lcptools -L$(HOME)/.local/lib -llcptools -o your_program
```

If you want to link dynamic library, please use as follows:

```sh
g++ your_program.cpp -I$(HOME)/.local/include/lcptools -L$(HOME)/.local/lib -llcptools -Wl,-rpath,$(HOME)/.local/lib -o your_program
```

**Note**: Make sure that paths are correct.

## Character Encoding

The binary encoding of the alphabet is defined as follows. This default encoding is used unless a custom encoding is provided:

| Character | Binary Encoding |
| --------- | --------------- |
| A, a      | 00              |
| T, t      | 11              |
| G, g      | 10              |
| C, c      | 01              |

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
    
    Ex: $w = xyza_0 . . . a_nklm$, where $n \geq 1$ and $xyz$ and $klm$ are identified as cores, and $z \lt a_0 \lt \dots \lt a_n \lt k$ or $z \gt a_0 \gt \dots \gt a_n \gt k$.

[<img src="https://github.com/BilkentCompGen/lcptools/blob/main/img/lcp-visual.png?raw=true" width="50%"/>](https://github.com/BilkentCompGen/lcptools/blob/main/img/lcp-visual.png?raw=true)

**Processing a string using LCP.** Here, blue underlines the core that satisfies the Local Minimum core, green represents a Local Maximum core,  red corresponds to a Repetitive Interior core, and yellow denotes a Stranded Sequence core.

### Deterministic Coin Tossing:

The dct function in the LCP algorithm is crucial for processing binary sequences. It starts by pinpointing the initial point of difference between two binary strings, beginning from the right-end. The function then assesses the difference based on the position and value of the divergent bit. This detail is transformed into a new binary sequence, which establishes the foundation of a newly generated 'core'. This core is a clear representation of the differences between the original sequences, integral to the algorithm's deepening process. Essentially, the dct function effectively consolidates and encapsulates the information, ensuring efficient further analysis within the LCP framework.

[<img src="https://github.com/BilkentCompGen/lcptools/blob/main/img/dct-visual.png?raw=true" width="50%"/>](https://github.com/BilkentCompGen/lcptools/blob/main/img/dct-visual.png?raw=true)

**DCT for reducing bitstreams to a new alphabet.** Each block above corresponds to the bitstream of a core. The DCT compares each core's bitstream to its left neighbor to form a shorter alphabet. E.g., the least significant bit of the core bitstream *11101011*, which differs from its left neighbor *011011* at the fourth index (counting from the right and starting at index 0). The value of this bit is *0*. Therefore, DCT replaces the core bitstream *11101011* with the concatenation of the bits of *4*, and the value of the differing bit, resulting in *1000*. The figure shows the reduced alphabet inferred by the DCT for each core bitstream, on which LCP can be applied.

### Overview of LCP

[<img src="https://github.com/BilkentCompGen/lcptools/blob/main/img/lcp-overview-visual.png?raw=true" width="75%"/>](https://github.com/BilkentCompGen/lcptools/blob/main/img/lcp-overview-visual.png?raw=true)

**Overview of the LCP with DCT and Labeling Paradigm.**

### Deepen Function:

The deepen function in the LCP algorithm primarily focuses on the compression of 'cores' alongside their left neighbors. The purpose of this repeated compression (dct) is to manage the length of the cores, preventing them from becoming large. After a compression, the LCP algorithm is re-applied. This re-application aims to identify new cores within the compressed data. In this context, each compressed core is treated as a discrete value, represented in binary form. This representation facilitates efficient processing and analysis within the algorithm.

This function iteratively compresses and processes cores to find new cores in compliance with the rules stated above.

## Experimental Results of LCP in CHM13v2.0

| **LCP Level**         | **1**         | **2**         | **3**         | **4**         | **5**         | **6**         | **7**         | **8**         |
|-----------------------|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|
| **Total #**           | 1,385,807,259 | 595,674,177   | 253,384,737   | 108,531,321   | 46,426,965    | 19,885,308    | 8,514,640     | 3,646,315     |
| **Unique**            | 1,395         | 517,688       | 159,701,016   | 93,601,379    | 42,859,941    | 18,795,370    | 8,191,525     | 3,563,854     |
| **Exec. Time (sec)**  | 120.12        | 104.07        | 48.80         | 20.23         | 8.98          | 3.64          | 1.53          | 0.64          |
| |
| **Avg Dist.**         | 2.25          | 5.23          | 12.30         | 28.72         | 67.14         | 156.76        | 366.10        | 854.87        |
| **StdDev Dist.**      | 1.04          | 1.85          | 4.21          | 8.74          | 19.95         | 45.24         | 104.24        | 241.08        |
| **Avg Len.**          | 3.67          | 10.47         | 26.58         | 63.84         | 151.23        | 354.89        | 830.88        | 1,941.85      |
| **StdDev Len.**       | 1.01          | 1.90          | 4.32          | 8.83          | 19.49         | 43.02         | 97.57         | 222.96        |
| |
| **Decr. in Core #**   | 0.44          | 0.43          | 0.43          | 0.43          | 0.43          | 0.43          | 0.43          | 0.43          |
| **Incr. in Avg Len.** | 3.67          | 2.85          | 2.54          | 2.40          | 2.37          | 2.35          | 2.34          | 2.34          |
| **Incr. in Avg Dist.**| 2.25          | 2.33          | 2.35          | 2.33          | 2.34          | 2.33          | 2.34          | 2.34          |
| |
| **Total Size (GB)**   | 56.79         | 24.41         | 10.38         | 4.45          | 1.90          | 0.81          | 0.35          | 0.15          |

### Notes
- **Dataset:** The [CHM13v2.0](https://www.science.org/doi/10.1126/science.abj6987) genome was used to generate LCP statistics with a single iteration of DCT.
- **Decr. in Core #:** Proportional reduction in the number of cores compared to the previous level.
- **Incr. in Avg Len.:** Proportional increase in the average core length compared to the previous level.
- **Incr. in Avg Dist.:** Proportional increase in the average distance between the start positions of consecutive cores compared to the previous level.
- As the LCP level increases:
  - The decrease in core count tends to **0.43**.
  - The average length and distance increase tends to **1/0.43 ≈ 2.33**.

## Default Variables

The default iteration count for compression in each deepening is set to 1.

