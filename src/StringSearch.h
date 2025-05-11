#ifndef KLUE_STRINGSEARCH_H
#define KLUE_STRINGSEARCH_H

#include "common.h"
#include "Kmer.hpp"
#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <map>
#include <set>

struct ContigInfo { // 16 bytes
  int16_t color; // Color ID of contig = 2 bytes
  char s[MAX_KMER_SIZE*2]; // Holds more sequence info 62 bits = 8 bytes
  uint8_t l; // Holds sequence length = 1 byte
  uint8_t rule; // How to process the sequence when found [padding for now]; = 1 byte
  bool fwd; // Strandedness; = 1 byte
  bool curr_color; // Whether it's found in current file iteration; = 1 byte
  uint16_t colors_found; // How many colors (sequence files) it's found in; = 2 bytes
};

class AhoCorasick {
private:
  int dictionary_index;

  void search(const char* corpus, size_t len, std::vector<ContigInfo*>& info_vec);
  void processWord(const char* corpus, size_t len, const Kmer& km, int position, std::vector<ContigInfo*>& info_vec);

public:
  AhoCorasick();
  ~AhoCorasick();
  AhoCorasick( const AhoCorasick& ) = delete;
  AhoCorasick& operator=( const AhoCorasick& ) = delete;
  void searchInCorpus(const char* corpus, size_t len, std::vector<ContigInfo*>& info_vec);
  void add(const std::string& contig, uint16_t color);
  u_map_<Kmer, ContigInfo, KmerHash> infomap;
};

#endif // KLUE_STRINGSEARCH_H

