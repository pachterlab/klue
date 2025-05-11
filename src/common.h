#ifndef KLUE_COMMON_H
#define KLUE_COMMON_H

#define KLUE_VERSION "0.29.0"

// NOTE: MAKE SURE THIS FILE GETS INCLUDED FIRST IN ALL OTHER FILES AND BEFORE ANY EXTERNAL LIBRARIES

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include "kseq.h"

#if defined(_MSVC_LANG)
#define KLUE_CPP_VERSION _MSVC_LANG
#else
#define KLUE_CPP_VERSION __cplusplus
#endif
#if KLUE_CPP_VERSION < 201703L
#include "robin_hood.h"
#define u_map_ robin_hood::unordered_flat_map
#define u_set_ robin_hood::unordered_set
#else
#include "unordered_dense.h"
#define u_map_ ankerl::unordered_dense::map
#define u_set_ ankerl::unordered_dense::set
#endif


#ifdef _WIN64
typedef unsigned int uint;
#endif

struct ProgramOptions {
  bool verbose;
  int threads;
  std::string distinguish_output_fasta;
  std::string bubble_left_output_fasta;
  std::string bubble_right_output_fasta;
  std::vector<std::string> bubble_variation_output_fasta;
  std::string input_fasta_contig;
  std::string map_file;
  std::string input_set_operations;
  std::string tmp_dir;
  int k;
  int g;
  int min_found_colors;
  bool stream_out;
  bool distinguish_all_but_one_color;
  bool distinguish_union;
  bool distinguish_combinations;
  bool extend;
  bool bubble;
  bool flanking_bubble;
  int distinguish_range_begin;
  int distinguish_range_end;
  std::vector<int> kmer_multiplicity;
  std::vector<int> inner;
  std::vector<std::string> transfasta;

ProgramOptions() :
  verbose(false),
  threads(1),
  k(31),
  g(0),
  stream_out(false),
  distinguish_all_but_one_color(false),
  distinguish_union(false),
  distinguish_combinations(false),
  extend(false),
  bubble(false),
  flanking_bubble(false),
  distinguish_range_begin(0),
  distinguish_range_end(0),
  min_found_colors(-1),
  tmp_dir("tmp")
  {}
};

std::string pretty_num(size_t num);
std::string pretty_num(int64_t num);
std::string pretty_num(int num);


#ifdef KLUE_USE_ZLIB_NG
#include "zlib-ng/zlib.h"
#else
#include <zlib.h>
#endif

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif
  
extern std::string revcomp(const std::string& s);


#endif // KLUE_COMMON_H
