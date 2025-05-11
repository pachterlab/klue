#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
#include <thread>
#include <time.h>
#include <algorithm>
#include <limits>

#include <cstdio>

#include "common.h"
#include "ProcessReads.h"
#include "KmerIndex.h"
#include <CompactedDBG.hpp>

//#define ERROR_STR "\033[1mError:\033[0m"
#define ERROR_STR "Error:"

using namespace std;

int my_mkdir(const char *path, mode_t mode) {
  #ifdef _WIN64
  return mkdir(path);
  #else
  return mkdir(path,mode);
  #endif
}

bool checkFileExists(std::string fn) {
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0;
}


std::string revcomp(const std::string& s) {
  std::string r(s);
  std::transform(s.rbegin(), s.rend(), r.begin(), [](char c) {
    switch(c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'a': return 'T';
    case 'c': return 'G';
    case 'g': return 'C';
    case 't': return 'A';
    default: return 'N';
    }
    return 'N';
  });
  return r;
}


void ParseOptionsDistinguish(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int pipe_flag = 0;
  int distinguish_all_flag = 0;
  int distinguish_all_but_one_flag = 0;
  int distinguish_combinations_flag = 0;
  int extension_flag = 0;
  int bubble_flag = 0;
  int flanking_bubble_flag = 0;
  const char *opt_string = "o:L:R:V:k:m:t:r:M:g:ps:T:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"all", no_argument, &distinguish_all_flag, 1},
    {"all-but-one", no_argument, &distinguish_all_but_one_flag, 1},
    {"combinations", no_argument, &distinguish_combinations_flag, 1},
    {"extend", no_argument, &extension_flag, 1},
    {"bubble", no_argument, &bubble_flag, 1},
    {"flanking_bubble", no_argument, &flanking_bubble_flag, 1},
    // short args
    {"output", required_argument, 0, 'o'},
    {"left", required_argument, 0, 'L'},
    {"right", required_argument, 0, 'R'},
    {"variation", required_argument, 0, 'V'},
    {"pipe", no_argument, &pipe_flag, 'p'},
    {"kmer-size", required_argument, 0, 'k'},
    {"kmer-multiplicity", required_argument, 0, 'M'},
    {"min-size", required_argument, 0, 'm'},
    {"map-file", required_argument, 0, 'g'},
    {"threads", required_argument, 0, 't'},
    {"tmp", required_argument, 0, 'T'},
    {"range", required_argument, 0, 'r'},
    {"set", required_argument, 0, 's'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 'o': {
      opt.distinguish_output_fasta = optarg;
      break;
    }
    case 'L': {
      opt.bubble_left_output_fasta = optarg;
	  break;
	}
    case 'R': {
        opt.bubble_right_output_fasta = optarg;
        break;
    }
    case 'V': {
        std::string str;
        std::stringstream(optarg) >> str;
        std::stringstream ss(str);
        while (ss.good()) {
            std::string s;
            getline(ss, s, ',');
            opt.bubble_variation_output_fasta.push_back(s);
        }
        break;
    }
    case 'p': {
      pipe_flag = 1;
      break;
    }
    case 'k': {
      stringstream(optarg) >> opt.k;
      break;
    }
    case 'g': {
      opt.map_file = optarg;
      break;
    }
    case 's': {
      opt.input_set_operations = optarg;
      break;
    }
    case 'm': {
      stringstream(optarg) >> opt.g;
      break;
    }
    case 'M': {
      std::string str;
      stringstream(optarg) >> str;
      std::stringstream ss(str);
      while(ss.good()) {
        std::string s;
        getline(ss, s, ',');
        opt.kmer_multiplicity.push_back(std::atoi(s.c_str()));
      }
      break;
    }
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'T': {
      stringstream(optarg) >> opt.tmp_dir;
      break;
    }
    case 'r': {
      std::string range_input_str;
      stringstream(optarg) >> range_input_str;
      stringstream ss(range_input_str);
      std::string range_val;
      int i = 0;
      char delimeter = '-';
      if (range_input_str.find(',') < range_input_str.length()) {
        delimeter = ','; // If string contains commas, use commas as delimeter
      }
      while (std::getline(ss, range_val, delimeter)) { 
        if (i == 0) opt.distinguish_range_begin = std::atoi(range_val.c_str());
        if (i == 1) opt.distinguish_range_end = std::atoi(range_val.c_str());
        i++;
      }
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
  if (pipe_flag) {
    opt.stream_out = true;
  }
  
  if (distinguish_all_flag) {
    opt.distinguish_union = true;
  }
  if (distinguish_all_but_one_flag) {
    opt.distinguish_all_but_one_color = true;
  }
  if (distinguish_combinations_flag) {
    opt.distinguish_combinations = true;
  }
  if (extension_flag) {
    opt.extend = true;
  }
  if (bubble_flag) {
	opt.bubble = true;
	opt.stream_out = true;
  }
  if (flanking_bubble_flag) {
	opt.bubble = true;
	opt.flanking_bubble = true;
	opt.stream_out = true;
  }
  for (int i = optind; i < argc; i++) {
    opt.transfasta.push_back(argv[i]);
  }
  
  if (opt.kmer_multiplicity.size() == 0) {
    opt.kmer_multiplicity.resize(opt.transfasta.size(), 1);
  }
}


void ParseOptionsRefineUnitigs(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int pipe_flag = 0;
  int distinguish_all_flag = 0;
  int distinguish_all_but_one_flag = 0;
  const char *opt_string = "o:L:R:V:k:m:t:r:pT:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"all", no_argument, &distinguish_all_flag, 1},
    {"all-but-one", no_argument, &distinguish_all_but_one_flag, 1},
    // short args
    {"output", required_argument, 0, 'o'},
    {"left", required_argument, 0, 'L'},
    {"right", required_argument, 0, 'R'},
    {"variation", required_argument, 0, 'V'},
    {"pipe", no_argument, &pipe_flag, 'p'},
    {"kmer-size", required_argument, 0, 'k'},
    {"min-size", required_argument, 0, 'm'},
    {"threads", required_argument, 0, 't'},
    {"tmp", required_argument, 0, 'T'},
    {"range", required_argument, 0, 'r'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);
    
    if (c == -1) {
      break;
    }
    
    switch (c) {
    case 0:
      break;
    case 'o': {
      opt.distinguish_output_fasta = optarg;
      break;
    }
    case 'L': {
        opt.bubble_left_output_fasta = optarg;
        break;
    }
    case 'R': {
        opt.bubble_right_output_fasta = optarg;
        break;
    }
    case 'V': {
        std::string str;
        std::stringstream(optarg) >> str;
        std::stringstream ss(str);
        while (ss.good()) {
            std::string s;
            getline(ss, s, ',');
            opt.bubble_variation_output_fasta.push_back(s);
        }
        break;
    }
    case 'p': {
      pipe_flag = 1;
      break;
    }
    case 'k': {
      stringstream(optarg) >> opt.k;
      break;
    }
    case 'm': {
      stringstream(optarg) >> opt.g;
      break;
    }
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'T': {
      stringstream(optarg) >> opt.tmp_dir;
      break;
    }
    case 'r': {
      std::string range_input_str;
      stringstream(optarg) >> range_input_str;
      stringstream ss(range_input_str);
      std::string range_val;
      int i = 0;
      char delimeter = '-';
      if (range_input_str.find(',') < range_input_str.length()) {
        delimeter = ','; // If string contains commas, use commas as delimeter
      }
      while (std::getline(ss, range_val, delimeter)) { 
        if (i == 0) opt.distinguish_range_begin = std::atoi(range_val.c_str());
        if (i == 1) opt.distinguish_range_end = std::atoi(range_val.c_str());
        i++;
      }
      break;
    }
    default: break;
    }
  }
  
  if (verbose_flag) {
    opt.verbose = true;
  }
  if (pipe_flag) {
    opt.stream_out = true;
  }
  
  if (distinguish_all_flag) {
    opt.distinguish_union = true;
  }
  if (distinguish_all_but_one_flag) {
    opt.distinguish_all_but_one_color = true;
  }
  
  for (int i = optind; i < argc; i++) {
    opt.transfasta.push_back(argv[i]);
  }
}


void ParseOptionsRefine(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int pipe_flag = 0;
  opt.k = 29; // New default for k-mer size
  const char *opt_string = "o:L:R:V:i:I:t:r:m:k:pT:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    // short args
    {"output", required_argument, 0, 'o'},
    {"pipe", no_argument, &pipe_flag, 'p'},
    {"left", required_argument, 0, 'L'},
    {"right", required_argument, 0, 'R'},
    {"variation", required_argument, 0, 'V'},
    {"input", required_argument, 0, 'i'},
    {"inner", required_argument, 0, 'I'},
    {"threads", required_argument, 0, 't'},
    {"tmp", required_argument, 0, 'T'},
    {"range", required_argument, 0, 'r'},
    {"min", required_argument, 0, 'm'},
    {"kmer-size", required_argument, 0, 'k'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);
    
    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 'o': {
      opt.distinguish_output_fasta = optarg;
      break;
    }
    case 'L': {
        opt.bubble_left_output_fasta = optarg;
        break;
    }
    case 'R': {
        opt.bubble_right_output_fasta = optarg;
        break;
    }
    case 'V': {
        std::string str;
        std::stringstream(optarg) >> str;
        std::stringstream ss(str);
        while (ss.good()) {
		  std::string s;
		  getline(ss, s, ',');
		  opt.bubble_variation_output_fasta.push_back(s);
		}
        break;
    }
    case 'p': {
      pipe_flag = 1;
      break;
    }
    case 'm': {
      stringstream(optarg) >> opt.min_found_colors;
      break;
    }
    case 'k': {
      stringstream(optarg) >> opt.k;
      break;
    }
    case 'i': {
      opt.input_fasta_contig = optarg;
      break;
    }
    case 'I': {
      std::string str;
      stringstream(optarg) >> str;
      std::stringstream ss(str);
      while(ss.good()) {
        std::string s;
        getline(ss, s, ',');
        opt.inner.push_back(std::atoi(s.c_str()));
      }
      break;
    }
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'T': {
      stringstream(optarg) >> opt.tmp_dir;
      break;
    }
    case 'r': {
      std::string range_input_str;
      stringstream(optarg) >> range_input_str;
      stringstream ss(range_input_str);
      std::string range_val;
      int i = 0;
      char delimeter = '-';
      if (range_input_str.find(',') < range_input_str.length()) {
        delimeter = ','; // If string contains commas, use commas as delimeter
      }
      while (std::getline(ss, range_val, delimeter)) { 
        if (i == 0) opt.distinguish_range_begin = std::atoi(range_val.c_str());
        if (i == 1) opt.distinguish_range_end = std::atoi(range_val.c_str());
        i++;
      }
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
  if (pipe_flag) {
    opt.stream_out = true;
  }
  
  for (int i = optind; i < argc; i++) {
    opt.transfasta.push_back(argv[i]);
  }
}

bool CheckOptionsDistinguish(ProgramOptions& opt) {

  bool ret = true;
  
  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
      opt.threads = n;
    }
  }

  if (opt.k <= 1 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error: invalid k-mer length " << opt.k << ", minimum is 3 and maximum is " << (MAX_KMER_SIZE -1) << endl;
    ret = false;
  }

  if (opt.k % 2 == 0) {
    cerr << "Error: k needs to be an odd number" << endl;
    ret = false;
  }

  if (opt.transfasta.empty()) {
    cerr << "Error: no FASTA files specified" << endl;
    ret = false;
  } else {

    for (auto& fasta : opt.transfasta) {
      // we want to generate the index, check k, index and transfasta
      struct stat stFileInfo;
      auto intStat = stat(fasta.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: FASTA file not found " << fasta << endl;
        ret = false;
      }
    }
  }

  // TODO file validation for --bubble, also make sure that -V option and #colors agree
  if (opt.distinguish_output_fasta.empty() && !opt.stream_out) {
    cerr << "Error: need to specify output FASTA file name or use --pipe" << endl;
    ret = false;
  } else if (!opt.distinguish_output_fasta.empty() && opt.stream_out) {
    cerr << "Error: cannot supply both output FASTA file name and --pipe" << endl;
    ret = false;
  }
  if (opt.stream_out && opt.verbose) {
    cerr << "Error: cannot supply both --verbose and --pipe" << endl;
    ret = false;
  }
  
  if (opt.distinguish_union && opt.distinguish_all_but_one_color) {
    cerr << "Error: Cannot use multiple distinguish options" << endl;
    ret = false;
  }
  
  if (opt.kmer_multiplicity.size() != opt.transfasta.size()) {
    cerr << "Error: k-mer multiplicities must match the number of input files" << endl;
    ret = false;
  } else {
    for (auto m : opt.kmer_multiplicity) {
      if (m != 1 && m != 2) {
        cerr << "Error: k-mer multiplicity must be either 1 or 2" << endl;
        ret = false;
        break;
      }
    }
  }

  if (opt.g != 0) {
    if (opt.g <= 2 || opt.g > opt.k - 2) {
      cerr << "Error: invalid minimizer size " << opt.g << ", minimum is 3 and maximum is k - 2" << endl;
      ret = false;
    }
  }

  if (opt.distinguish_range_end == 0) {
    opt.distinguish_range_end = opt.distinguish_range_begin;
  }
  if (opt.distinguish_range_end < opt.distinguish_range_begin) {
    cerr << "Error: invalid range supplied" << endl;
    ret = false;
  }
  if (opt.bubble && opt.bubble_variation_output_fasta.size() != opt.transfasta.size()) {
      cerr << "Error: Number of var files invalid" << endl;
      ret = false;
  }
  return ret;
}



bool CheckOptionsRefineUnitigs(ProgramOptions& opt) {
  
  bool ret = true;
  
  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
      opt.threads = n;
    }
  }

  if (opt.k <= 1 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error: invalid k-mer length " << opt.k << ", minimum is 3 and maximum is " << (MAX_KMER_SIZE -1) << endl;
    ret = false;
  }

  if (opt.k % 2 == 0) {
    cerr << "Error: k needs to be an odd number" << endl;
    ret = false;
  }

  if (opt.transfasta.empty()) {
    cerr << "Error: no FASTA files specified" << endl;
    ret = false;
  } else {
    for (auto& fasta : opt.transfasta) {
      // we want to generate the index, check k, index and transfasta
      struct stat stFileInfo;
      auto intStat = stat(fasta.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: FASTA file not found " << fasta << endl;
        ret = false;
      }
    }
  }
  
  if (opt.distinguish_output_fasta.empty() && !opt.stream_out) {
    cerr << "Error: need to specify output FASTA file name or use --pipe" << endl;
    ret = false;
  } else if (!opt.distinguish_output_fasta.empty() && opt.stream_out) {
    cerr << "Error: cannot supply both output FASTA file name and --pipe" << endl;
    ret = false;
  }
  if (opt.stream_out && opt.verbose) {
    cerr << "Error: cannot supply both --verbose and --pipe" << endl;
    ret = false;
  }
  
  if (opt.distinguish_union && opt.distinguish_all_but_one_color) {
    cerr << "Error: Cannot use multiple distinguish options" << endl;
    ret = false;
  }
  
  if (opt.g != 0) {
    if (opt.g <= 2 || opt.g > opt.k - 2) {
      cerr << "Error: invalid minimizer size " << opt.g << ", minimum is 3 and maximum is k - 2" << endl;
      ret = false;
    }
  }
  
  if (opt.distinguish_range_end == 0) {
    opt.distinguish_range_end = opt.distinguish_range_begin;
  }
  if (opt.distinguish_range_end < opt.distinguish_range_begin) {
    cerr << "Error: invalid range supplied" << endl;
    ret = false;
  }
  
  return ret;
}

bool CheckOptionsRefine(ProgramOptions& opt) {

  bool ret = true;
  
  if (opt.k <= 1 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error: invalid k-mer length " << opt.k << ", minimum is 3 and maximum is " << (MAX_KMER_SIZE -1) << endl;
    ret = false;
  }
  
  if (opt.k % 2 == 0) {
    cerr << "Error: k needs to be an odd number" << endl;
    ret = false;
  }

  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
      opt.threads = n;
    }
  }
  
  if (opt.transfasta.empty()) {
    cerr << "Error: no FASTA files specified" << endl;
    ret = false;
  } else {
    for (auto& fasta : opt.transfasta) {
      struct stat stFileInfo;
      auto intStat = stat(fasta.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: FASTA file not found " << fasta << endl;
        ret = false;
      }
    }
  }
  
  if (opt.input_fasta_contig.size() == 0) {
    cerr << "Error: no input contig FASTA file specified" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.input_fasta_contig.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << "Error: FASTA file not found " << opt.input_fasta_contig << endl;
      ret = false;
    }
  }

  if (opt.distinguish_output_fasta.empty() && !opt.stream_out) {
    cerr << "Error: need to specify output FASTA file name or use --pipe" << endl;
    ret = false;
  } else if (!opt.distinguish_output_fasta.empty() && opt.stream_out) {
    cerr << "Error: cannot supply both output FASTA file name and --pipe" << endl;
    ret = false;
  }
  if (opt.stream_out && opt.verbose) {
    cerr << "Error: cannot supply both --verbose and --pipe" << endl;
    ret = false;
  }

  if (opt.distinguish_range_end == 0) {
    opt.distinguish_range_end = opt.distinguish_range_begin;
  }
  if (opt.distinguish_range_end < opt.distinguish_range_begin) {
    cerr << "Error: invalid range supplied" << endl;
    ret = false;
  }

  return ret;
}


void PrintCite() {
  cout << "When using this program in your research, please cite" << endl << endl
       << "  [todo]" << endl
       << endl;
}

void PrintVersion() {
  cout << "klue, version " << KLUE_VERSION << endl;
}

void usage() {
  cout << "klue " << KLUE_VERSION << endl << endl
       << "Usage: klue <CMD> [arguments] .." << endl << endl
       << "Where <CMD> can be one of:" << endl << endl
       << "    distinguish   Extracts distinguishing contigs from FASTA/FASTQ files "<< endl
       << "    refine        Refine contigs based on certain criteria"<< endl
       << "    version       Prints version information" << endl
       << "    cite          Prints citation information" << endl << endl
       << "Running klue <CMD> without arguments prints usage information for <CMD>"<< endl << endl;
}

void usageDistinguish() {
  cout << "klue " << KLUE_VERSION << endl
       << "Extracts distinguishing contigs from FASTA/FASTQ files" << endl << endl
       << "Usage: klue distinguish [arguments] FASTA-files" << endl << endl
       << "Required argument (choose one of the following):" << endl
       << "-p, --pipe                 Direct output to standard output" << endl
       << "-o, --output=STRING        Filename for the output FASTA" << endl << endl
       << "Optional argument:" << endl
       << "-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: " << (MAX_KMER_SIZE-1) << ")" << endl
       << "-M, --kmer-multiplicity=INT Number of times a k-mer must be encountered (default: 1)" << endl
       << "                            Can specify multiple numbers (one for each input file)" << endl
       << "-r, --range=INT-INT         Set the range of of length of output sequences (format: begin-end)" << endl
       << "                            Can specify multiple numbers (one for each input file)" << endl
       << "    --all                   Set the mode to be extracting sequences found in any input" << endl
       << "    --all-but-one           Set the mode to be extracting all sequences except those found across all inputs" << endl
       << "    --combinations          Set the mode to be listing every combination of sets" << endl
       << "    --extend                Extend the contigs" << endl
       << "    --bubble                Extract bubble sequences" << endl
       << "-g, --map-file=STRING       Filename for mapping output FASTA headers to file names" << endl
       << "-t, --threads=INT           Number of threads to use (default: 1)" << endl
       << "-m, --min-size=INT          Length of minimizers (default: automatically chosen)" << endl
       << "-T, --tmp=STRING            Directory for temporary files" << endl
       << "Bubble-specific options:" << endl
       << "-L, --left=FILE             File to output the left-bubble-end sequences in" << endl
       << "-R, --right=FILE            File to output the right-bubble-end sequences sequences in" << endl
       << "-V, --var=STRING            Files (separated by commas) to output the variations (i.e. sequences within bubbles)" << endl
       << "                            (number of var files should equal number of FASTA-files, i.e. number of colors)" << endl
       << endl;

}

void usageRefineUnitigs() {
  cout << "klue " << KLUE_VERSION << endl
       << "Reconstruct contigs in FASTA file, outputting extended contig sequences" << endl << endl
       << "Usage: klue reconstruct [arguments] FASTA-files" << endl << endl
       << "Required argument (choose one of the following):" << endl
       << "-p, --pipe                 Direct output to standard output" << endl
       << "-o, --output=STRING        Filename for the output FASTA" << endl << endl
       << "Optional argument:" << endl
       << "-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: " << (MAX_KMER_SIZE-1) << ")" << endl
       << "-r, --range=INT-INT         Set the range of of length of output sequences (format: begin-end)" << endl
       << "                            Can specify multiple numbers (one for each input file)" << endl
       << "    --all                   Set the mode to be extracting sequences found in any input" << endl
       << "    --all-but-one           Set the mode to be extracting all sequences except those found across all inputs" << endl
       << "-t, --threads=INT           Number of threads to use (default: 1)" << endl
       << "-m, --min-size=INT          Length of minimizers (default: automatically chosen)" << endl
       << "-T, --tmp=STRING            Directory for temporary files" << endl
       << endl;
  
}

void usageRefine() {
  cout << "klue " << KLUE_VERSION << endl
       << "Refines contigs in FASTA file, outputting the ones that meet certain criteria" << endl << endl
       << "Usage: klue refine [arguments] FASTA-files" << endl << endl
       << "Required argument (choose one of the following):" << endl
       << "-p, --pipe                 Direct output to standard output" << endl
       << "-i, --input=STRING         Filename for the input contig FASTA" << endl
       << "-o, --output=STRING        Filename for the output FASTA" << endl << endl
       << "Optional argument:" << endl
       << "-k, --kmer-size=INT         k-mer (odd) length (default: 29, max value: " << (MAX_KMER_SIZE-1) << ")" << endl
       << "-I, --inner=INT             The length of the inner region between two flanking regions of a contig" << endl
       << "-r, --range=INT-INT         Set the range of of length of output sequences (format: begin-end)" << endl
       << "-t, --threads=INT           Number of threads to use (default: 1)" << endl
       << "-m, --min=INT               Minimum number of files contig must be found in (default: number of FASTA-files)" << endl
       << "-T, --tmp=STRING            Directory for temporary files" << endl
       << endl;
  
}

std::string argv_to_string(int argc, char *argv[]) {
  std::string res;
  for (int i = 0; i < argc; ++i) {
    res += argv[i];
    if (i + 1 < argc) {
      res += " ";
    }
  }

  return res;
}

std::string get_local_time() {
  time_t rawtime;
  struct tm * timeinfo;

  time( &rawtime );
  timeinfo = localtime( &rawtime );
  std::string ret(asctime(timeinfo));

  // chomp off the newline
  return ret.substr(0, ret.size() - 1);
}

int main(int argc, char *argv[]) {
  std::cout.sync_with_stdio(false);
  setvbuf(stdout, NULL, _IOFBF, 1048576);

  if (argc < 2) {
    usage();
    exit(1);
  } else {
    auto start_time(get_local_time());
    ProgramOptions opt;
    string cmd(argv[1]);
    if (cmd == "version") {
      PrintVersion();
    } else if (cmd == "cite") {
      PrintCite();
    } else if (cmd == "distinguish") {
      if (argc==2) {
        usageDistinguish();
        return 0;
      }
      ParseOptionsDistinguish(argc-1,argv+1,opt);
      if (!CheckOptionsDistinguish(opt)) {
        usageDistinguish();
        exit(1);
      } else {
        Kmer::set_k(opt.k);
        KmerIndex index(opt);
        index.BuildDistinguishingGraph(opt, opt.transfasta);
        if (!opt.map_file.empty()) { // Write out mapping file (a type of t2g file)
          std::ofstream out_map_f;
          out_map_f.open(opt.map_file);
          for (int i = 0; i < opt.transfasta.size(); i++) {
            out_map_f << std::to_string(i) << "\t" << opt.transfasta[i] << "\n";
          }
          out_map_f.close();
        }
      }
    } else if (cmd == "reconstruct") {
      if (argc==2) {
        usageRefineUnitigs();
        return 0;
      }
      ParseOptionsRefineUnitigs(argc-1,argv+1,opt);
      if (!CheckOptionsRefineUnitigs(opt)) {
        usageRefineUnitigs();
        exit(1);
      } else {
        Kmer::set_k(opt.k);
        KmerIndex index(opt);
        index.BuildReconstructionGraph(opt);
      }
    } else if (cmd == "refine") {
      if (argc==2) {
        usageRefine();
        return 0;
      }
      ParseOptionsRefine(argc-1,argv+1,opt);
      if (!CheckOptionsRefine(opt)) {
        usageRefine();
        exit(1);
      } else {
        if (!opt.stream_out) opt.verbose = true;
        Kmer::set_k(opt.k);
        MasterProcessor MP(opt);
        int numreads = ProcessReads(MP, opt);
        FILE* out;
        auto f = opt.distinguish_output_fasta;
        if (opt.stream_out) {
          out = stdout;
        } else {
          out = fopen(f.c_str(), "wb");
        }
        MP.writeContigs(out, opt.min_found_colors < 0 ? opt.transfasta.size() : opt.min_found_colors);
        fflush(out);
        fflush(stdout);
        fclose(out);
      }
    } else {
      cerr << "Error: invalid command " << cmd << endl;
      usage();
      exit(1);
    }
    rmdir(opt.tmp_dir.c_str()); // Remove temp directory if non-empty
  }

  fflush(stdout);

  return 0;
}
