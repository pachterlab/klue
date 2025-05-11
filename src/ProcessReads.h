#ifndef KLUE_PROCESSREADS_H
#define KLUE_PROCESSREADS_H

#include "common.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include "StringSearch.h"

class MasterProcessor;

int64_t ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt);

class SequenceReader {
public:
  
  SequenceReader(const ProgramOptions& opt) :
  readbatch_id(-1) {};
  SequenceReader() : state(false), readbatch_id(-1) {};
  virtual ~SequenceReader() {}
  
  virtual bool empty() = 0;
  virtual void reset();
  virtual void reserveNfiles(int n) = 0;
  virtual bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                              std::vector<std::pair<const char*, int>>& names,
                              std::vector<std::pair<const char*, int>>& quals,
                              std::vector<uint32_t>& flags,
                              int &readbatch_id,
                              bool full=false,
                              bool comments=false) = 0;
  
  
public:
  bool state; // is the file open
  int readbatch_id = -1;
};

class FastqSequenceReader : public SequenceReader {
public:
  
  FastqSequenceReader(const ProgramOptions& opt, const std::vector<std::string>& _files) : 
  SequenceReader(opt), current_file(0), files(_files) {
    SequenceReader::state = false;
    nfiles = _files.size();
    reserveNfiles(nfiles);
  }
  FastqSequenceReader() : SequenceReader(), 
  current_file(0) {};
  FastqSequenceReader(FastqSequenceReader &&o);
  ~FastqSequenceReader();
  
  bool empty();  
  void reset();
  void reserveNfiles(int n);
  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<uint32_t>& flags,
                      int &readbatch_id,
                      bool full=false,
                      bool comments=false);
  
public:
  uint32_t numreads = 0;
  std::vector<gzFile> fp;
  std::vector<int> l;
  std::vector<int> nl;
  std::vector<std::string> files;
  int current_file;
  std::vector<kseq_t*> seq;
  int interleave_nfiles;
  int nfiles = 1;
};

class MasterProcessor {
public:
  MasterProcessor (const ProgramOptions& opt)
    : opt(opt), numreads(0), numcontigs(0), bufsize(1ULL<<23), curr_readbatch_id(0) { 

    readSeqs = false;
    std::vector<std::string> in_files;
    in_files.push_back(opt.input_fasta_contig);
    inSR = new FastqSequenceReader(opt, in_files);
    verbose = opt.verbose;
    nfiles = opt.transfasta.size();
    rangefilteredcount = 0;
    numchars = 0;
  }
  
  ~MasterProcessor() {
    if (SR != nullptr) delete SR;
    if (inSR != nullptr) delete inSR;
  }
  
  std::mutex reader_lock;
  std::mutex writer_lock;
  std::condition_variable cv;

  bool verbose;
  bool readSeqs; // If we should start reading FASTA sequences, not contigs
  
  FastqSequenceReader *SR; // Reading the FASTA files
  FastqSequenceReader *inSR; // Reading the input contigs FASTA

  const ProgramOptions& opt;
  int64_t numreads, numcontigs;
  size_t bufsize;
  int nfiles;
  int curr_readbatch_id;
  size_t numchars, numchars_;
  std::atomic<int> rangefilteredcount;
  std::vector<ContigInfo*> info_vec;
  
  AhoCorasick ac;
  
  void processReads();
  void processContigs();
  void update(int n,
              size_t numchars,
              std::vector<ContigInfo*>& _info_vec,
              std::vector<std::pair<const char*, int>>& seqs,
              std::vector<std::pair<const char*, int>>& names,
              std::vector<std::pair<const char*, int>>& quals,
              std::vector<uint32_t>& flags,
              int readbatch_id);
  void writeContigs(FILE* out, int min_colors_found);
};

class ReadProcessor {
public:
  ReadProcessor(const ProgramOptions& opt, MasterProcessor& mp, int file_no = -1);
  ReadProcessor(ReadProcessor && o);
  ~ReadProcessor();
  char *buffer;
  
  size_t bufsize;
  MasterProcessor& mp;
  int64_t numreads;
  
  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<uint32_t> flags;
  std::vector<ContigInfo*> info_vec;
  
  bool full;
  bool comments;
  int file_no;
  size_t numchars, numchars_;
  
  /*std::vector<std::vector<int>> newIDs;
   std::vector<std::vector<int>> IDs;*/
  
  void operator()();
  void processBuffer();
  void processBufferContigs();
  void clear();
};

std::string pretty_num(size_t num);

#endif // KLUE_PROCESSREADS_H
