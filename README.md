# klue

**klue** is a C++ program for pseudoassembly with k-mers. The program uses the [Bifrost](https://github.com/pmelsted/Bifrost) colored de Bruijn graph for indexing. For details on the method and its applications, see the preprint describing klue: [Pseudoassembly of k-mers](https://doi.org/10.1101/2025.05.11.653354) by Delaney Sullivan, Mayuko Boffelli, and Lior Pachter.

## Prerequisites

To build `klue`, you need the following:

- A C++17-compatible compiler (e.g., AppleClang 15+, GCC 7+, or newer)
- [Bifrost](https://github.com/pmelsted/Bifrost) installed with:
  - Headers in `/usr/local/include/bifrost`
  - Library in `/usr/local/lib/libbifrost.a` or `.dylib`
- [Zlib](https://zlib.net) (preinstalled on macOS and most Linux distributions)

## Build and Install

Clone the repository and build using CMake:

```bash
git clone https://github.com/pachterlab/klue.git
cd klue
mkdir build
cd build
cmake ..
make -j
sudo make install
```

## Usage

After installation, run `klue` from the command line:

```bash
klue 0.29.0

Usage: klue <CMD> [arguments] ..

Where <CMD> can be one of:

    distinguish   Extracts distinguishing contigs from FASTA/FASTQ files 
    refine        Refine contigs based on certain criteria
    version       Prints version information
    cite          Prints citation information
```
