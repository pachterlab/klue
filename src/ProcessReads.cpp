#include <fstream>
#include <limits>
#include <iomanip>
#include "common.h"
#include "ProcessReads.h"
#include "kseq.h"
#include <unordered_set>

std::string pretty_num(size_t num) {
  auto s = std::to_string(num);
  auto ret = std::string("");
  
  if (s.size() <= 3) {
    return s;
  }
  
  int remainder = s.size() % 3;
  if (remainder == 0) {
    remainder = 3;
  }
  
  size_t start_pos = 0;
  while (start_pos + remainder < s.size() - 1) {
    ret += s.substr(start_pos, remainder) + ",";
    start_pos += remainder;
    remainder = 3;
  }
  
  ret += s.substr(start_pos, 3);
  
  return ret;
}


//methods


int64_t ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt) {
  std::ios_base::sync_with_stdio(false);
  size_t numreads = 0;
  
  if (MP.verbose) {
    std::cerr << "* processing the contigs ..."; std::cerr.flush();
  }
  MP.processContigs();
  numreads = MP.numreads;
  if (MP.verbose) {
    std::cerr << std::endl << "done " << std::endl;
  }
  if (MP.verbose && MP.rangefilteredcount > 0) {
    std::cerr << "* " << MP.rangefilteredcount << " contigs filtered out due to length" << std::endl;
  }
  if (MP.verbose) {
    std::cerr << "* processed " << pretty_num(numreads) << " contigs";
    std::cerr << std::endl;
  }
  
  MP.numcontigs = numreads;
  
  if (MP.verbose) {
    std::cerr << "* processing the sequences ..."; std::cerr.flush();
  }
  MP.processReads();
  numreads = MP.numreads;
  if (MP.verbose) {
    std::cerr << std::endl << "done " << std::endl;
  }
  if (MP.verbose) {
    std::cerr << "* processed " << pretty_num(numreads) << " sequences";
    std::cerr << std::endl;
  }

  return numreads;
}


/** -- read processors -- **/

void MasterProcessor::processContigs() {
  
  numreads = 0;
  curr_readbatch_id = 0;

  // start worker threads
  
  std::vector<std::thread> workers;
  
  // Set number of threads to 1 because race condition when storing contigs
  
  for (int i = 0; i < /*opt.threads*/1; i++) {
    workers.emplace_back(std::thread(ReadProcessor(opt,*this)));
  }
    
  // let the workers do their thing
  for (int i = 0; i < /*opt.threads*/1; i++) {
    workers[i].join(); //wait for them to finish
  }
  delete inSR;
  inSR = nullptr;
}


void MasterProcessor::processReads() {
  
  numreads = 0;
  curr_readbatch_id = 0;
  readSeqs = true;
  
  // start worker threads
  
  for (int j = 0; j < nfiles; j++) {
    std::vector<std::thread> workers;
    std::vector<std::string> _files;
    if (verbose) {
      std::cerr << std::endl;
      std::cerr << "* Processing file #" << j << ": " << opt.transfasta[j] << std::endl;
    }
    _files.push_back(opt.transfasta[j]);
    SR = new FastqSequenceReader(opt, _files); // New sequence reader for next file!
    for (int i = 0; i < opt.threads; i++) {
      workers.emplace_back(std::thread(ReadProcessor(opt,*this, j)));
    }
    // let the workers do their thing
    for (int i = 0; i < opt.threads; i++) {
      workers[i].join(); //wait for them to finish
    }
    delete SR; // Delete current sequence reader
    SR = nullptr;
    // Iterate through info_vec and update its counts and color book-keeping
    for (auto &info_ptr : info_vec) {
      if ((*info_ptr).curr_color) {
        (*info_ptr).colors_found++;
        (*info_ptr).curr_color = false; // Reset it to false (accounts for the fact that the same sequence might be found multiple times in corpus)
      }
    }
    info_vec.clear(); // Clear it
  }
}

void MasterProcessor::update(int n, 
                             size_t b,
                             std::vector<ContigInfo*>& _info_vec,
                             std::vector<std::pair<const char*, int>>& seqs,
                             std::vector<std::pair<const char*, int>>& names,
                             std::vector<std::pair<const char*, int>>& quals,
                             std::vector<uint32_t>& flags,
                             int readbatch_id) {
  // acquire the writer lock
  std::unique_lock<std::mutex> lock(this->writer_lock);

  //while (readbatch_id != curr_readbatch_id) {
    //cv.wait(lock, [this, readbatch_id]{ return readbatch_id == curr_readbatch_id; });
  //}

  numreads += n;
  numchars += b;

  info_vec.insert(info_vec.end(), _info_vec.begin(), _info_vec.end());
  _info_vec.clear();
  
  //curr_readbatch_id++;
  lock.unlock(); // releases the lock
  //cv.notify_all(); // Alert all other threads to check their readbatch_id's!
}

void MasterProcessor::writeContigs(FILE* out, int min_colors_found) {
  size_t i = 0;
  for (const auto& entry : ac.infomap) {
    const auto& info = entry.second;
    if (info.colors_found >= min_colors_found) {
      i++;
      std::string header = ">" + std::to_string(info.color);
      fwrite(header.c_str(), sizeof(char), header.length(), out);
      fwrite("\n", sizeof(char), 1, out);
      fwrite(info.fwd ? entry.first.toString().c_str() : entry.first.twin().toString().c_str(), sizeof(char), entry.first.k, out);
      fwrite(info.s, sizeof(char), info.l, out);
      fwrite("\n", sizeof(char), 1, out);
    }
    // Debug:
    // std::cout << "DEBUG: " << info.s << " colors_found=" << info.colors_found.size() << " " << std::endl;
  }
  if (verbose) {
    std::cerr << "* " << pretty_num(i) << " contigs retained and written to output" << std::endl;
  }
}

ReadProcessor::ReadProcessor(const ProgramOptions& opt, MasterProcessor& mp, int file_no) : 
  mp(mp), numreads(0), file_no(file_no) {
  // initialize buffer
  bufsize = mp.bufsize;
  buffer = new char[bufsize];
  seqs.reserve(bufsize/50);
  numchars = 0;
  numchars_ = 0;
  clear();
}

ReadProcessor::ReadProcessor(ReadProcessor && o) :
  bufsize(o.bufsize),
  mp(o.mp),
  numreads(o.numreads),
  seqs(std::move(o.seqs)),
  names(std::move(o.names)),
  quals(std::move(o.quals)),
  flags(std::move(o.flags)),
  full(o.full),
  comments(o.comments),
  file_no(o.file_no),
  numchars(o.numchars), 
  numchars_(o.numchars_){
  buffer = o.buffer;
  o.buffer = nullptr;
  o.bufsize = 0;
}

ReadProcessor::~ReadProcessor() {
  if (buffer != nullptr) {
    delete[] buffer;
    buffer = nullptr;
  }
}

void ReadProcessor::operator()() {
  while (true) {
    int readbatch_id;
    std::vector<std::string> umis;
    FastqSequenceReader *SR;
    bool full = false;
    bool comments = false;
    if (mp.readSeqs) {
      SR = mp.SR;
    } else {
      full = true; // We need to get the name
      SR = mp.inSR;
    }
    // grab the reader lock
    {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (SR->empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        SR->fetchSequences(buffer, bufsize, seqs, names, quals, flags, readbatch_id, full, comments);
        if (seqs.size() == 0) { // resize buffer if overflown
          delete[] buffer;
          bufsize *= 2;
          buffer = new char[bufsize];
        }
      }
      // release the reader lock
    }
    
    // process our sequences
    if (mp.readSeqs) {
      processBuffer();
    } else {
      processBufferContigs();
    }

    // update the results, MP acquires the lock
    int nfiles = SR->nfiles;
    mp.update(seqs.size() / nfiles, numchars_, info_vec, seqs, names, quals, flags, readbatch_id);
    clear();
  }
}

void ReadProcessor::processBuffer() {
  // actually process the sequence
  
  int incf, jmax, nfiles;
  nfiles = mp.nfiles;
  incf = 0;
  jmax = 1;
  
  std::vector<const char*> s(jmax, nullptr);
  std::vector<int> l(jmax,0);

  for (int i = 0; i + incf < seqs.size(); i++) {
    for (int j = 0; j < jmax; j++) {
      s[j] = seqs[i+j].first;
      l[j] = seqs[i+j].second;
      numchars += l[j];
      numchars_ += l[j];
      // DEBUG:
      //std::cout << std::to_string(i) << ": " << std::string(s[j], l[j]) << std::endl;
      mp.ac.searchInCorpus(s[j], l[j], info_vec);
    }
    i += incf;
    numreads++;
    if (((numreads > 0 && numreads % 1000000 == 0) || (numchars > 0 && numchars % 100000000 == 0)) && mp.verbose) {
      numreads = 0; // reset counter
      numchars = 0; // reset counter
      std::cerr << '\r' << (mp.numreads/1000000) << "M reads, " << (mp.numchars/100000000) << "00M bases processed";
      std::cerr << "                             ";
      std::cerr.flush();
    }
  }
}

void ReadProcessor::processBufferContigs() {
  // actually process the sequence
  
  int incf, jmax, nfiles;
  nfiles = 1; //nfiles = mp.nfiles;
  incf = nfiles-1;
  jmax = nfiles;
  uint32_t rb = std::max(mp.opt.distinguish_range_begin,0); // range begin filter
  uint32_t re = mp.opt.distinguish_range_end == 0 ? rb : std::max(mp.opt.distinguish_range_end,0); // range end filter
  if (rb == 0 && re == 0) re = std::numeric_limits<uint32_t>::max();

  for (int i = 0; i + incf < seqs.size(); i++) {
    for (int j = 0; j < jmax; j++) {
      // Debug:
      // std::cout << std::string(names[i+j].first, names[i+j].second) << " " << std::string(seqs[i+j].first, seqs[i+j].second) << std::endl;
      uint32_t len = seqs[i+j].second;
      if (len >= rb && len <= re) {
        mp.ac.add(seqs[i+j].first, std::atoi(std::string(names[i+j].first, names[i+j].second).c_str())); // race condition, so this function should NOT be multithreaded
      } else {
        mp.rangefilteredcount++;
      }
    }
    i += incf;
    numreads++;
    if (numreads > 0 && numreads % 1000000 == 0 && mp.verbose) {
      numreads = 0; // reset counter
      std::cerr << '\r' << (mp.numreads/1000000) << "M contigs processed";
      std::cerr << "         ";
      std::cerr.flush();
    }
  }
}

void ReadProcessor::clear() {
  memset(buffer,0,bufsize);
  numchars_ = 0;
}

/** -- sequence readers -- **/

void SequenceReader::reset() {
  state = false;
  readbatch_id = -1;
}

FastqSequenceReader::~FastqSequenceReader() {
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
  }
  
  for (auto &s : seq) {
    kseq_destroy(s);
  }
}


bool FastqSequenceReader::empty() {
  return (!state && current_file >= files.size());
}

void FastqSequenceReader::reset() {
  SequenceReader::reset();
  
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
    f = nullptr;
  }
  
  for (auto &ll : l) {
    ll = 0;
  }
  for (auto &nll : nl) {
    nll = 0;
  }
  
  current_file = 0;
  for (auto &s : seq) {
    kseq_destroy(s);
    s = nullptr;
  }
}

void FastqSequenceReader::reserveNfiles(int n) {
  fp.resize(nfiles);
  l.resize(nfiles, 0);
  nl.resize(nfiles, 0);
  seq.resize(nfiles, nullptr);
}

// returns true if there is more left to read from the files
bool FastqSequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
                                         std::vector<std::pair<const char *, int> > &names,
                                         std::vector<std::pair<const char *, int> > &quals,
                                         std::vector<uint32_t>& flags,
                                         int& read_id,
                                         bool full,
                                         bool comments) {
  
  std::string line;
  readbatch_id += 1; // increase the batch id
  read_id = readbatch_id; // copy now because we are inside a lock
  seqs.clear();
  if (full) {
    names.clear();
    quals.clear();
  }
  flags.clear();
  
  int bufpos = 0;
  int count = 0; // for interleaving
  int pad = nfiles;
  while (true) {
    if (!state) { // should we open a file
      if (current_file >= files.size()) {
        // nothing left
        return false;
      } else {
        // close the current files
        for (auto &f : fp) {
          if (f) {
            gzclose(f);
          }
        }
        
        // open the next one
        for (int i = 0; i < nfiles; i++) {
          fp[i] = files[0] == "-" && nfiles == 1 && files.size() == 1 ? gzdopen(fileno(stdin), "r") : gzopen(files[current_file+i].c_str(), "r");
          seq[i] = kseq_init(fp[i]);
          l[i] = kseq_read(seq[i]);
          
        }
        current_file+=nfiles;
        state = true; 
      }
    }
    // the file is open and we have read into seq1 and seq2
    int bufadd = nfiles;
    int num_finished = 0;
    for (int i = 0; i < nfiles; i++) {
      if (l[i] < 0) num_finished++;
      bufadd += l[i]; // includes seq
    }
    if (num_finished != nfiles) {      
      // fits into the buffer
      if (full) {
        for (int i = 0; i < nfiles; i++) {
          if (l[i] < 0) continue;
          nl[i] = seq[i]->name.l + (comments ? seq[i]->comment.l+1 : 0);
          bufadd += nl[i]; // includes name
          //bufadd += l[i]; // include qual
        }
        bufadd += 2*pad;
      }
      
      if (bufpos+bufadd< limit) {
        if (interleave_nfiles != 0) { // Hack to allow interleaving
          // assert(nfiles == 1);
          if (bufpos+bufadd >= limit-262144 && count % interleave_nfiles == 0) {
            return true;
          }
          count++;
        }
        
        for (int i = 0; i < nfiles; i++) {
          if (l[i] < 0) continue;
          char *pi = buf + bufpos;
          memcpy(pi, seq[i]->seq.s, l[i]+1);
          bufpos += l[i]+1;
          seqs.emplace_back(pi,l[i]);
          
          if (full && !comments) {
            //pi = buf + bufpos;
            //memcpy(pi, seq[i]->qual.s,l[i]+1);
            //bufpos += l[i]+1;
            //quals.emplace_back(pi,l[i]);
            pi = buf + bufpos;
            memcpy(pi, seq[i]->name.s, nl[i]+1);
            bufpos += nl[i]+1;
            names.emplace_back(pi, nl[i]);
          } else if (full && comments) {
            //pi = buf + bufpos;
            //memcpy(pi, seq[i]->qual.s,l[i]+1);
            //bufpos += l[i]+1;
            //quals.emplace_back(pi,l[i]);
            pi = buf + bufpos;
            memcpy(pi, seq[i]->name.s, (nl[i]-(seq[i]->comment.l+1)));
            names.emplace_back(pi, nl[i]);
            bufpos += (nl[i]-(seq[i]->comment.l+1));
            pi = buf + bufpos;
            const char* blank_space = " ";
            memcpy(pi, blank_space, 1);
            bufpos += 1;
            pi = buf + bufpos;
            memcpy(pi, seq[i]->comment.s, seq[i]->comment.l+1);
            bufpos += seq[i]->comment.l+1;
          }
        }
        
        numreads++;
        flags.push_back((current_file-nfiles) / nfiles); // flags.push_back(numreads-1);
      } else {
        if (interleave_nfiles != 0) {
          std::cerr << "Error: There was an error processing interleaved FASTQ input. Exiting..." << std::endl;
          exit(1);
        }
        return true; // read it next time
      }
      
      // read for the next one
      for (int i = 0; i < nfiles; i++) {
        l[i] = kseq_read(seq[i]);
      }        
    } else {
      state = false; // haven't opened file yet
    }
  }
}

// move constructor

FastqSequenceReader::FastqSequenceReader(FastqSequenceReader&& o) :
  nfiles(o.nfiles),
  numreads(o.numreads),
  fp(std::move(o.fp)),
  l(std::move(o.l)),
  nl(std::move(o.nl)),
  files(std::move(o.files)),
  current_file(o.current_file),
  seq(std::move(o.seq)) {
  
  o.fp.resize(nfiles);
  o.l.resize(nfiles, 0);
  o.nl.resize(nfiles, 0);
  o.seq.resize(nfiles, nullptr);
  o.state = false;
}

