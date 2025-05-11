#include "StringSearch.h"
#include "KmerIterator.hpp"
#include <stdexcept>

  // TODO: check threading error for contig processing
  // TODO: reduce children to A,T,C,G (4 TrieNode* pointers)
  // TODO: see if we can eliminate std::vector<std::string> matchedWords
  // TODO: see if we can store DNA as 2-bits
  // TODO: char counter (1M chars processed instead of 1M reads processed)
  // TODO: naive dictionary lookups (if we only have no more than 3 variations)
  // TODO: better AC implementation
  // TODO: prevent buffer from going out of control (or do system other than jmax)
  // TODO: implement other features (other flanks)
  // TODO: Check contig lengths (allow longer contigs and shorter contigs, and output them, but just return false from the add() function to notify)
  // TODO: Deal with color_found logic
  // TODO: Deal with corpus to handle N's that are stored in the Kmer
  // TODO: strandedness
  // TODO: bubble
  // Alternative: Build a colored de bruijn graph of 29-mers over contigs (w/ the next 29-mer stored in node and 3-inbetween chars) and ensure the sliding window works [and remove N's]
  // Microsatellites: k-mer must cross [f][r][r][f] entirely to be unique, put [r][r] in dictionary, and look up and extend to find [f].
  // // GATATATATATATATATATATATATATATC; others can be shorter (e.g. GATATATC) but others can't be longer (violates uniqueness rule)
  // // // how might the flanks be unique?
  // // // 31-bp of AT's uniqueness = long unique in one and only short ones elsewhere; need to find the long one, extract flanking 31-mer, find flanking in others
  // // // Contig of GGGGATATATATATATATATATATCCCC = how unique? 

void AhoCorasick::processWord(const char* corpus, size_t len, const Kmer& km, int position, std::vector<ContigInfo*>& info_vec) {
  auto it = infomap.find(km.rep());
  if (it == infomap.end()) return; // Nothing more to do
  // TODO: forward working, reverse not working
  const char* corpus_ = corpus+position; // move to position
  auto& info = it->second;
  /*auto l = info.l;
  auto contig_size = Kmer::k + l;
  bool lex = (km == km.rep());
  bool found_fwd = lex;
  auto& contig_substr = info.s;
  int middle_start = Kmer::k; // 29
  int middle_end = contig_size-Kmer::k; // We consider the middle segment within [start,end), e.g. [29,32)
  int middle_segment_length = middle_end-middle_start;
  int flank, flank_left, flank_right;
  bool odd_flank = false;
  //std::cout << "A: " << middle_start << " " << middle_segment_length << " " << middle_end << " : " << contig_size << " " << std::to_string(l) << std::endl;
  if (middle_segment_length < 0) return; // No way to have flanks
  if (middle_segment_length % 2 == 0) {
    flank = 0;
    int middle_two = (middle_start+middle_end)/2; // e.g. position 30 in [29,31)
    int middle_one = middle_two-1; // e.g. position 29 in [29,31)
    flank_left = middle_one-flank;
    flank_right = middle_two+flank;
    //std::cout << "B: " << middle_start << " " << middle_one << " " << middle_two << " " << middle_end << " " << flank_left << " " << flank_right << std::endl;
    
    // [29,33) = 29,30,31,32
  } else {
    odd_flank = true;
    flank = 1;
    int middle = (middle_start+middle_end-1)/2; // e.g. position 30 in [29,32)
    flank_left = middle-flank; // e.g. position 29
    flank_right = middle+flank; // e.g. position 31
    std::cout << "C: " << middle_start << " " << middle << " " << middle_end << " " << flank_left << " " << flank_right << std::endl;
  }
  
  // Debug:
  std::cout << "Word: " << km.toString() << ", Canonical Word: " << km.rep().toString() << ", Position: " << position << ", Color: " << info.color << ", Rest of word: " << std::string(contig_substr, l) << ", Strand: " << std::to_string(info.fwd) << ", Corpus: " << std::string(corpus, len) << std::endl;

  bool success = false;
  if ((found_fwd && info.fwd) || (!found_fwd && !info.fwd)) {
    
    // Pretend AAAAA is canonical and TTTTT is non-canonical
    // Reference: TTTTT GGGG (TTTTT = !info.fwd; how it exists in reference; exists as AAAAA in reference but stored in hash table as TTTTT)
    // Read:      TTTTT GGGG (TTTTT = !found_fwd; exists as TTTTT in read but, again, is non-canonical)
    // 
    // Read:      AAAAA GGGG (AAAAA = !found_fwd; how it exists in read; also non-canonical)
    // Therefore: Do not to reverse complement TTTTT
    // 
    
    // Check if the surroundings actually match between contig and sequencing read
    //std::cout << "1: " << std::to_string(middle_start) << " " << std::string(contig_substr,l) << " " << std::to_string(flank_left) << " " << std::to_string(middle_start) << std::endl; // 2: 29 CATAC?????GTTGTG 28672 29
    if (strncmp(corpus_+middle_start, contig_substr, flank_left - middle_start + 1) != 0) return; // left flank
    //std::cout << "2: " << flank_right << " " << flank_left << " " << middle_start << " " << middle_end << " " << flank << " :" <<  (flank_left-middle_start+flank) << std::endl;
    // 2: 17 16 9 25 0 :7
    // [29]M---L----R---M
    //     ^cs
    // middle_end=100-31 = 69 = 69 - flank_right = 69 - 17 = 52     vs.     0+31 = 31 = middle_start  ; 19+19
    // 34: 34-17= 17
    if (strncmp(corpus_+flank_right, contig_substr-middle_start+flank_right, middle_end-flank_right) != 0) return; // right flank
    //std::cout << "3" << std::endl;
    if (strncmp(corpus_+middle_end, contig_substr+middle_segment_length, l-middle_segment_length) != 0) return; // contig substring right of right flank
    //std::cout << "4" << std::endl;
    // Now check "inside" the flank
    bool non_ATCG = false;
    const char* c = corpus_;
    for (size_t i_c = 0; i_c < contig_size; c++, i_c++) {
      if (*c != 'A' && *c != 'T' && *c != 'C' && *c != 'G' && *c != 'a' && *c != 't' && *c != 'c' && *c != 'g') {
        non_ATCG = true; // No N's allowed in the region in the corpus segment overlapped by the contig
      }
    }
    
    c = corpus_+flank_left+1;
    for (size_t i_c = 0; i_c < flank; c++, i_c++)  std::cout << (*c); // Debug

    if (!non_ATCG) {
      success = true;
      std::cout << "\n^Success FWD match" << std::endl; // Debug
    }
  } else {
    // Pretend AAAAA is canonical and TTTTT is non-canonical
    // Reference: TTTTT GGGG (TTTTT = !info.fwd; how it exists in reference; exists as AAAAA in reference but stored in hash table as TTTTT)
    // Read: CCCC AAAAA      (AAAAA = found_fwd; exists as AAAAA in read canonically)
    // 
    // Therefore: Need to reverse complement the GGGG (how it's stored in the reference) to be CCCC (how it appears in the read)
    std::string revcomp_s = revcomp(std::string(contig_substr, l)); // Need to store the string as a separate variable so the const char* doesn't point to something out of bounds
    const char* contig_substr_revcomp = revcomp_s.c_str();
    std::cout << "1" << " " << contig_substr_revcomp << " " << std::to_string(l) << " " << revcomp_s.c_str() << std::endl;
    if (strncmp(corpus_+middle_start, contig_substr_revcomp, flank_left - middle_start + 1) == 0) return; // left flank
    std::cout << "2: " << std::to_string(middle_start) << " " << contig_substr_revcomp << " " << flank_right << " " << flank_left << " " << middle_end << " " << middle_segment_length << " " << flank << " " << std::to_string(l) << std::endl;
    // 2: 9 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN 21 19 32 23 1 32
    // contig: AAA[ATCCCC] (k=3)
    // sequen: CCCCCATATCCGCACCCCCTGCGGATATGGGGTATAATGCGAAACGCATTATACCCC[GGGGAT]TTTAAAAAAAATCCCC
    // 2: middle_start=3 GGGGAT flank_right=5 flank_left=3 middle_end=6 middle_segment_length=3 flank=1 l=6
    // Might be better to just reverse the sequence and treat it as above
    std::cout << (corpus_+flank_right)[0]; // corpus_ = where the TTT starts
    std::cout << (corpus_+flank_right)[1];
    std::cout << (corpus_+flank_right)[2];
    std::cout << (corpus_+flank_right)[3];
    std::cout << (corpus_+flank_right)[4] << std::endl; // AAAAA
    std::cout << (contig_substr_revcomp-middle_start+flank_right)[0];
    std::cout << (contig_substr_revcomp-middle_start+flank_right)[1];
    std::cout << (contig_substr_revcomp-middle_start+flank_right)[2] << std::endl; // GGA
    if (strncmp(corpus_+flank_right, contig_substr_revcomp-middle_start+flank_right, middle_end-flank_right) != 0) return; // right flank
    std::cout << "3" << std::endl; // Nope, doesn't work, only 1 2 printed
    
    bool non_ATCG = false;
    const char* c = corpus_;
    for (size_t i_c = 0; i_c < contig_size; c++, i_c++) {
      if (*c != 'A' && *c != 'T' && *c != 'C' && *c != 'G' && *c != 'a' && *c != 't' && *c != 'c' && *c != 'g') {
        non_ATCG = true; // No N's allowed in the region in the corpus segment overlapped by the contig
      }
    }
    // TODO: 
    if (!non_ATCG) {
      success = true;
      std::cout << "\n^Success FWD match" << std::endl; // Debug
    }
  }*/

  //if (success) {
    // Don't do below because calling vector::size() is very slow
    /*int sz = info_vec.size();
    if (sz >= 4096) {
      info_vec.reserve(sz*1.2); // grow slowly in capacity
    }*/
    info.curr_color = true;
    info_vec.push_back(&info);
  //}
}

void AhoCorasick::search(const char* corpus, size_t len, std::vector<ContigInfo*>& info_vec) {
  KmerIterator kit(corpus), kit_end;
  for (int i = 0;  kit != kit_end; ++i,++kit) {
    processWord(corpus, len, kit->first, i, info_vec);
  }
}

AhoCorasick::AhoCorasick() {
}

void AhoCorasick::add(const std::string& contig, uint16_t color) {
  if (contig.length() > (Kmer::k + MAX_KMER_SIZE*2)) {
    throw std::invalid_argument("Contig too long, must be less than " + std::to_string(Kmer::k + MAX_KMER_SIZE*2) + ": " + contig);
    return;
  } else if (contig.length() < Kmer::k) {
    throw std::invalid_argument("Contig short long, must be at least " + std::to_string(Kmer::k) + ": " + contig);
    return;
  }
  Kmer km(contig.c_str()); // 29-mer
  int flank = Kmer::k;
  std::string word = contig.substr(0, flank);
  ContigInfo info;
  info.color = color;
  info.l = contig.length()-flank;
  memcpy(info.s, contig.c_str() + flank, info.l);
  info.rule = 0;
  bool lex = (km == km.rep()); // We want to index the canonical (i.e. lexicographically smaller) word
  info.fwd = lex; // If we're storing the sequence as-is (i.e. the non-reverse-complemented sequence is the canonical one)
  info.colors_found = 0;
  info.curr_color = false;
  infomap[km.rep()] = info; // Index into hash map
}

AhoCorasick::~AhoCorasick() {
}

void AhoCorasick::searchInCorpus(const char* corpus, size_t len, std::vector<ContigInfo*>& info_vec) {
  search(corpus, len, info_vec);
}
