#include <algorithm>
#include <random>
#include <sstream>
#include <ctype.h>
#include <unordered_set>
#include <functional>
#include "common.h"
#include "KmerIndex.h"
#include <iostream>
#include <unordered_map>
#include <string>
#include "ColoredCDBG.hpp"
#include <ExpressionParser.h>
#include <iomanip>
#include <queue>
#include <algorithm>
#include <cctype>

// other helper functions
// pre: u is sorted
bool isUnique(const std::vector<int>& u) {
    for (int j = 1; j < u.size(); j++) {
        if (u[j - 1] == u[j]) {
            return false;
        }
    }
    return true;
}

std::vector<int> unique(const std::vector<int>& u) {
    std::vector<int> v;
    v.reserve(u.size());
    v.push_back(u[0]);
    for (int j = 1; j < u.size(); j++) {
        if (u[j - 1] != u[j]) {
            v.push_back(u[j]);
        }
    }
    return v;
}

const char Dna(int i) {
    static const char* dna = "ACGT";
    return dna[i & 0x03];
}

int hamming(const char* a, const char* b) {
    int h = 0;
    while (*a != 0 && *b != 0) {
        if (*a != *b) {
            h++;
        }
        a++;
        b++;
    }
    return h;
}

int my_mkdir_kmer_index(const char* path, mode_t mode) {
#ifdef _WIN64
    return mkdir(path);
#else
    return mkdir(path, mode);
#endif
}

std::string generate_tmp_file(std::string seed, std::string tmp_dir) {
    struct stat stFileInfo;
    auto intStat = stat(tmp_dir.c_str(), &stFileInfo);
    if (intStat == 0) {
        // file/dir exits
        if (!S_ISDIR(stFileInfo.st_mode)) {
            cerr << "Error: file " << tmp_dir << " exists and is not a directory" << endl;
            exit(1);
        }
    }
    else {
        // create directory
        if (my_mkdir_kmer_index(tmp_dir.c_str(), 0777) == -1) {
            cerr << "Error: could not create directory " << tmp_dir << endl;
            exit(1);
        }
    }
    std::string base = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    std::string tmp_file = "klue.";
    srand((unsigned int)std::hash<std::string>{}(seed));
    int pos;
    while (tmp_file.length() < 32) {
        pos = ((rand() % (base.size() - 1)));
        tmp_file += base.substr(pos, 1);
    }
    return tmp_dir + "/" + tmp_file;
}

struct Extend {
    std::string result;
    int color;
};

Extend extendUnitig(const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current,
    std::unordered_set<std::string>& visited,
    const int& color,
    std::string current_str,
    const std::unordered_set<int>& superset_colors,
    const int& k) {

    Extend traversal_result;
    bool hasValidSuccessor = false;

    std::string colored_unitig = current.getUnitigHead().toString() + std::to_string(color);
    if (visited.find(colored_unitig) != visited.end() || current.isEmpty) {
        traversal_result.color = color;
        return traversal_result;
    }
    visited.insert(colored_unitig);

    auto successors = current.getSuccessors();
    if (successors.begin() == successors.end()) {
        traversal_result.color = color;
        return traversal_result;
    }

    for (const auto& next : successors) {
        if (visited.find(next.getUnitigHead().toString() + std::to_string(color)) == visited.end()) {
            UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
            UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
            std::unordered_set<int> colors;
            for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
            if (colors.find(color) != colors.end()) {
                hasValidSuccessor = true;
                std::string str;
                if (next.strand) {
                    for (int i = next.dist; i < next.len; ++i) {
                        str += next.getUnitigKmer(i).toString().substr(k - 1);
                    }
                }
                else {
                    for (int i = next.dist; i < next.len; ++i) {
                        str = next.getUnitigKmer(i).twin().toString().substr(k - 1) + str;
                    }
                }
                Extend extendResult = extendUnitig(next, visited, color, str, superset_colors, k);
                current_str += extendResult.result;
            }
        }
    }
    traversal_result.result = current_str;
    traversal_result.color = color;
    return traversal_result;
}

std::string generate_revcomp(const std::string& s) {
    std::string r(s.size(), ' ');
    std::transform(s.rbegin(), s.rend(), r.begin(), [](char c) {
        switch (c) {
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
        });
    return r;
}

// Check for (k-1) overlap between two sequences
bool overlap(const std::string& seq1, const std::string& seq2, int k) {
    if (k <= 1 || k - 1 > seq1.size() || k - 1 > seq2.size()) {
        return false;
    }
    return seq1.substr(seq1.size() - (k - 1)) == seq2.substr(0, (k - 1));
}

// Checks all permutations of pairs from the input sequences, accounting for both the sequences and their reverse complements
bool permute(const std::vector<std::string>& sequences, int k) {
    for (const auto& left : sequences) {
        for (const auto& right : sequences) {
            if (overlap(left, right, k) || overlap(generate_revcomp(left), right, k) ||
                overlap(left, generate_revcomp(right), k) || overlap(generate_revcomp(left), generate_revcomp(right), k)) {
                return true;
            }
        }
    }
    return false;
}


struct Bubble {
    std::string bubble_left;
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::string>>> variations; // maybe try another data structure?
    std::string bubble_right;
};

/*
 *
 */
void findPredecessorBaseNodes(ColoredCDBG<void>& ccdbg,
    Kmer curr_node,
    std::unordered_set<Kmer, KmerHash>& base_nodes,
    std::unordered_set<Kmer, KmerHash> visited_nodes,
    std::unordered_set<Kmer, KmerHash> visited_nodes_second_time,
    const std::unordered_set<int>& color_profile,
    const std::unordered_set<int>& ideal_color_profile) {
    const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current = ccdbg.find(curr_node); // Unitig associated with current node
    //std::cout << ":strand=" << std::to_string(current.strand) << " " << curr_node.toString() << std::endl;
    UnitigColors::const_iterator it_colors = current.getData()->getUnitigColors(current)->begin(current);
    UnitigColors::const_iterator it_colors_end = current.getData()->getUnitigColors(current)->end();
    std::unordered_set<int> color_profile_;
    for (; it_colors != it_colors_end; ++it_colors) { color_profile_.insert(it_colors.getColorID()); }

    if (color_profile_ == ideal_color_profile) {
        base_nodes.insert(curr_node.rep()); // Found a base node (anchoring the bubble on the left, aka predecessor, side)
        return;
    }
    if (color_profile_ != color_profile) { // Color switching; therefore we won't continue...
        return;
    }
    //std::cout << ":E " << curr_node.rep().toString() << std::endl;
    if (visited_nodes.find(curr_node.rep()) != visited_nodes.end()) { // Already found current node
        if (visited_nodes_second_time.find(curr_node.rep()) != visited_nodes_second_time.end()) { // Also already found it a second time!
            //std::cout << ":EEE " << std::to_string(current.strand) << " " << curr_node.toString() << std::endl;
            return; // If this is our third time visiting current node, we don't want to continue...
        }
        else {
            //std::cout << ":EE " << std::to_string(current.strand) << " " << curr_node.toString() << std::endl;
            visited_nodes_second_time.insert(curr_node.rep()); // Second time visiting current node, which is ok (we can loop once in the graph)
        }
    }
    else {
        visited_nodes.insert(curr_node.rep()); // First time visiting current node
    }
    auto predecessors_ = current.getPredecessors();
    if (predecessors_.begin() == predecessors_.end()) { // No more predecessors
        //std::cout << ": NOMORE" << std::endl;
        return;
    }
    for (const auto& prev : predecessors_) {
        findPredecessorBaseNodes(ccdbg, prev.getUnitigHead().rep(), base_nodes, visited_nodes, visited_nodes_second_time, color_profile, ideal_color_profile);
    }
}

void traverseBubble(ColoredCDBG<void>& ccdbg,
    Kmer curr_node,
    std::vector<std::vector<std::vector<std::string>>>& path, // outer index = color ID; inside vector = paths associated with color where each path is a vector of strings; so we have 1) color index, 2) path index, and 3) string index
    std::vector<std::string> tmp_path,
    int color,
    bool start,
    int cycle,
    int prev_strand,
    int rand_id,
    int k,
    std::unordered_set<Kmer, KmerHash> visited_nodes,
    std::unordered_set<Kmer, KmerHash> visited_nodes_second_time,
    const std::unordered_set<int>& ideal_color_profile,
    bool skip_successor_flag = false) {

    UnitigMap<DataAccessor<void>, DataStorage<void>, false> current = ccdbg.find(curr_node.rep()); // Unitig associated with current node
    if (current.isEmpty) { return; } // Empty mapping
    if (cycle > 3) { return; } // Set cycle limit
    current.strand = prev_strand;
    UnitigColors::const_iterator it_uc = current.getData()->getUnitigColors(current)->begin(current);
    UnitigColors::const_iterator uc_end = current.getData()->getUnitigColors(current)->end();

    std::string contig = "";
    int curr_pos = -1;
    std::unordered_map<Kmer, std::unordered_set<int>, KmerHash> kmers_color_map; // Key = k-mer; Value = Set of all colors associated with that k-mer
    bool skip_successor = skip_successor_flag;

    if (visited_nodes.find(curr_node.rep()) != visited_nodes.end()) { // Already found current node
        if (visited_nodes_second_time.find(curr_node.rep()) != visited_nodes_second_time.end()) { // Also already found it a second time!
            skip_successor_flag = true;
        }
        else {
            visited_nodes_second_time.insert(curr_node.rep()); // Second time visiting current node, which is ok (we can loop once in the graph)
        }
    }
    else {
        visited_nodes.insert(curr_node.rep()); // First time visiting current node
    }

    for (; it_uc != uc_end; ++it_uc) { // TODO: this is what we need to navigate well
        Kmer km = current.getUnitigKmer(it_uc.getKmerPosition());
        kmers_color_map[km.rep()].insert(it_uc.getColorID());
    }
    if (cycle == 0) {
        std::unordered_set<int> c; // Intersection
        for (auto e : ideal_color_profile) {
            if (kmers_color_map[curr_node.rep()].count(e)) { c.insert(e); }
        }
        if (c.size() != ideal_color_profile.size()) {
            return; // We only can start from unitigs that contain the same number of colors as the bubble ends
        }
    }

    std::unordered_set<int> color_profile;
    color_profile.insert(color); // This is the color of the path that we're currently navigating through

    bool end_of_bubble = false;
    bool initial_start_status = start;
    int i = 0;
    std::unordered_set<int> current_color;
    auto unitig_len_in_kmers = current.size - ccdbg.getK() + 1;
    
    while (true) {
        if (i == unitig_len_in_kmers) {
            break; // We've reached the end of the unitig so we're done
        }
        Kmer kmer = Kmer(current.referenceUnitigToString().c_str() + i);
        i++;
        if (kmers_color_map.find(kmer.rep()) != kmers_color_map.end()) {
            current_color = kmers_color_map[kmer.rep()]; // Update color because we're onto the next color of the unitig!
        }
        assert(!current_color.empty());
        std::string km = kmer.toString();
        if (color == -1) {
            if (!start) { // Make sure we have continuous segment of ideal_color_profile color profile before doing anything further
                if (current_color == ideal_color_profile) {
                    start = true;
                }
            }
            else if (current_color != ideal_color_profile) { // We're done with the beginning of the bubble
                // the next line only allows superbubble exploration (i.e. NO nested bubbles)
                //if (current_color.size() != 1)  return; // We only want one color for our path!                
                color = *(current_color.begin());
                color_profile.clear();
                color_profile.insert(color); // Now, we have our "official" path with our "official" color
                for (int i = 1; i < current_color.size(); ++i)  {
                    color_profile.insert(*(std::next(current_color.begin(), i)));
                }                
                tmp_path.push_back(contig); // Push back the contig we have so far (i.e. the beginning of the bubble)
                contig = ""; // reset contig
            }
        }
        else if (current_color == ideal_color_profile || end_of_bubble) {
            // Yay, we've encountered the end of the bubble!
            if (!end_of_bubble) {
                if (visited_nodes_second_time.find(kmer.rep()) != visited_nodes_second_time.end()) return; // Can't have a loop be the end of a bubble
                tmp_path.push_back(contig); // Push back the contig we have so far
                contig = ""; // reset contig
                end_of_bubble = true;
            }
            else if (current_color != ideal_color_profile) { // There's a color change within the unitig so we want to stop extending
                break; // Yup, stop right here
            }   
        }        
        else if (current_color != color_profile) {
            if (!contig.empty()) return; // Oops, ran into a color switch (aka a "dead end"), just return without updating anything
            color_profile.clear();
            for (const auto& c : current_color) {
                color_profile.insert(c);
                if (color_profile.size() == 1) color = c;
            }            
        }
        std::string _km = prev_strand ? kmer.twin().toString() : kmer.toString();
        if (contig == "") contig = _km;
        else {            
            // if last (k-1) of _km and first (k-1) of contig are the same, we want to take the first character of _km and append it to contig
            if (_km.substr(1) == contig.substr(0, k - 1)) {
                contig = _km[0] + contig;
            }
            // if last (k-1) of _km and last (k-1) of contig are the same, we want to take the last character of _km and append it to contig
            else {
                contig += _km[k - 1]; // k = _km.length()               
            }
        }
    }
    if (contig != "" && start) tmp_path.push_back(contig);
    if (end_of_bubble) {
        path[color].push_back(std::move(tmp_path)); // Note: Make sure we allocate in advance the path vector with the number of colors in advance
        return;
    }
    auto successors = current.getSuccessors();
    auto predecessors = current.getPredecessors();
    ++cycle;

    for (const auto& successor : successors) {
        bool move_next = true;
        Kmer next_kmer = curr_node;
        if (skip_successor) break;
        if (true) traverseBubble(ccdbg, successor.getUnitigHead().rep(), path, tmp_path, color, start, cycle, successor.strand, rand_id, k, visited_nodes, visited_nodes_second_time, ideal_color_profile, skip_successor_flag);
    }
    for (const auto& predecessor : predecessors) {
        break; // Don't do predecessors? 
        bool move_next = true;
        Kmer next_kmer = curr_node;
        if (move_next) traverseBubble(ccdbg, predecessor.getUnitigHead(), path, tmp_path, color, start, cycle, predecessor.strand, rand_id, k, visited_nodes, visited_nodes_second_time, ideal_color_profile);
    }
}

Bubble exploreBubble(ColoredCDBG<void>& ccdbg,
    const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current,
    std::unordered_set<std::string>& visited,
    std::unordered_set<Kmer, KmerHash> all_nodes,
    const std::unordered_set<int>& superset_colors,
    int color,
    int k,
    std::ostringstream& left_stream,
    std::vector<std::string>& var_stream,
    std::ostringstream& right_stream,
    std::mutex& visited_mutex) {

    if (current.isEmpty) { return {}; } // Empty node

    auto successors_ = current.getSuccessors();
    auto predecessors_ = current.getPredecessors();

    if (successors_.begin() == successors_.end() && predecessors_.begin() == predecessors_.end()) { return {}; } // Isolated node [no bubble structure]

    // Navigate to the very left anchor of bubble
    std::unordered_set<Kmer, KmerHash> visited_nodes;
    std::unordered_set<Kmer, KmerHash> visited_nodes_second_time;
    std::unordered_set<int> ideal_color_profile;
    for (int i = 0; i < ccdbg.getColorNames().size(); i++) {
        ideal_color_profile.insert(i);
    }
    Kmer curr_node = current.getUnitigHead();

    // Do the bubble traversion!
    std::vector<std::vector<std::vector<std::string>>> path;
    path.resize(ccdbg.getColorNames().size()); // Resize it to the number of colors
    std::vector<std::string> tmp_path;   

    traverseBubble(ccdbg, curr_node, path, tmp_path, -1, false, 0, -1, rand(), k, visited_nodes, visited_nodes_second_time, ideal_color_profile);

    for (int color = 0; color < path.size(); color++) {
        for (int path_i = 0; path_i < path[color].size(); path_i++) {
            int i = 0;
            while (i < path[color][path_i].size() - 2) { // stop before the last element
                auto x = path[color][path_i][i];
                // first/last element are source/sink
                // want to stitch together the middle elements
                std::string x_first = x.substr(0, k - 1); // first (k-1) characters
                std::string x_last;
                if (k - 1 > 0 && x.length() >= k - 1) {
                    x_last = x.substr(x.length() - (k - 1)); // last (k-1) characters
                }
                auto& next = path[color][path_i][i + 1];
                std::string next_first = next.substr(0, k - 1); // first (k-1) characters
                std::string next_last;
                std::string substr1;
                std::string substr2;
                if (k - 1 > 0 && next.length() >= k - 1) {
                    next_last = next.substr(next.length() - (k - 1)); // last (k-1) characters
                    substr1 = next.substr(0, next.length() - (k - 1));
                    substr2 = next.substr(k - 1);
                }
                if (x_first == next_last) {
                    std::string stitched_path = substr1 + x;
                    path[color][path_i][i + 1] = stitched_path;
                }
                else if (x_last == next_first) {
                    std::string stitched_path = x + substr2;
                    path[color][path_i][i + 1] = stitched_path;
                }
                i++;
            }
        }
    }

    // print for debugging
    /*
    for (int color = 0; color < path.size(); color++) {
         if (path[color].size() == 0) { continue; }
         std::cout << ":COLOR: " << color << std::endl;
         for (int path_i = 0; path_i < path[color].size(); path_i++) {
             std::cout << ":::";
             for (auto x : path[color][path_i]) {
                 std::cout << x << " ";
             }
             std::cout << std::endl;
         }
    }
    */

    std::string header;
    for (int color = 0; color < path.size(); color++) {
        header += std::to_string(color);
        if (color != path.size() - 1) { header += "_"; }
    }

    bool outputted_left_right = false;
    std::vector<std::string> visited_local;
    for (int color = 0; color < path.size(); color++) {
        bool outputted_color = false;
        if (path[color].empty()) { continue; }
        if (var_stream.size() == 1) { header = color; }
        for (int path_i = 0; path_i < path[color].size(); path_i++) {
            bool outputted_once = false;
            if (!outputted_once) {
                std::string left = path[color][path_i][0];
                std::string right = path[color][path_i].back();
                bool already_visited = false;
                {
                    std::lock_guard<std::mutex> lock(visited_mutex);
                    if (std::find_if(path[color][path_i].begin(), path[color][path_i].end(), [&](const std::string& node) {
                        return visited.count(node) || visited.count(revcomp(node));
                        }) != path[color][path_i].end()) {
                        already_visited = true;
                    }
		}
                if (!already_visited) {
                    {
                        std::lock_guard<std::mutex> lock(visited_mutex);
                        visited.insert(left);
                        visited.insert(right);
		    }
                    visited_local.push_back(left);
                    visited_local.push_back(right);
                }
                if (visited_local.empty()) continue;
                // Check if the bubble is valid
                std::vector<std::string> sequences = { path[color][path_i][0], path[color][path_i][path[color][path_i].size() - 3], path[color][path_i].back() };
                if (!permute(sequences, k)) continue; // is there (k-1) overlap between left/variation and right/variation sequences? 
                if (!outputted_left_right) left_stream << ">" << header << "\n" << path[color][path_i][0] << "\n"; // First element
                if (!outputted_color) var_stream[color] += ">" + std::to_string(color) + "\n" + path[color][path_i][path[color][path_i].size() - 3] + "\n"; // Stitched element
                outputted_color = true;
                if (!outputted_left_right) right_stream << ">" << header << "\n" << path[color][path_i].back() << "\n"; // Last element
                outputted_left_right = true; // We only want to output left/right once
                outputted_once = true;
            }
        }
    }
    return {};
}

// Begin set operations
// Split user-inputted set operation commands by " "
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Recursively compute set operations (union and intersection are nodes) that correspond to specific sets (leaves)
std::set<int> computeSetOperation(const Node* root, const std::map<int, std::set<int>>& k_map) {
    if (!root) { return {}; }
    if (root->value != 'U' && root->value != 'I' && root->value != '\\' && root->value != 'N' && root->value != 'X') {
        auto it = k_map.find(root->value - 'A'); // 'A' maps to 0, 'B' to 1, etc.
        return (it != k_map.end()) ? it->second : std::set<int>{};
    }
    // Construct the universal set from all sets in k_map
    std::set<int> universalSet;
    for (const auto& pair : k_map) {
        universalSet.insert(pair.second.begin(), pair.second.end());
    }
    auto leftSet = computeSetOperation(root->left, k_map);
    auto rightSet = computeSetOperation(root->right, k_map);
    std::set<int> result;
    if (root->value == 'U') { // (Logical AND) union of sets A,B. User-inputted as "AUB"
        std::set_union(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == 'I') { // (Logical OR) intersection of sets A,B. User-inputted as "AIB"
        std::set_intersection(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == '\\') { // Difference of sets A,B. User-inputted as "A\B"
        std::set_difference(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == 'N') { // Not both of sets A,B. User-inputted as "ANB"
        std::set<int> intersection;
        std::set_intersection(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(intersection, intersection.end()));
        std::set_difference(universalSet.begin(), universalSet.end(), intersection.begin(), intersection.end(), std::inserter(result, result.end()));
    }
    // DEBUG: add logical NOR operator -- the dual of NAND
    else if (root->value == 'X') { // Exclusive OR of sets A,B. User-inputted as "AXB"
        std::set<int> tempUnion, tempIntersection, xorResult;
        std::set_union(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(tempUnion, tempUnion.end()));
        std::set_intersection(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(tempIntersection, tempIntersection.end()));
        std::set_difference(tempUnion.begin(), tempUnion.end(), tempIntersection.begin(), tempIntersection.end(), std::inserter(xorResult, xorResult.end()));
        return xorResult;
    }
    return result;
}

std::string stringOpt(size_t n, std::string option) {
    std::string result;
    // Skip letters 'I', 'N', 'X', 'U' because they are used in set operations
    auto skipLetters = [&](char& ch) {
        if (ch >= 'I') { ch += 1; }
        if (ch >= 'N') { ch += 1; }
        if (ch >= 'X') { ch += 1; }
        if (ch >= 'U') { ch += 1; }
        };
    if (option == "default") {
        for (size_t i = 0; i < n; ++i) {
            char currentChar = 'A' + i;
            skipLetters(currentChar);
            result += currentChar;
            result += "\\(";
            for (size_t j = 0; j < n; ++j) {
                char innerChar = 'A' + j;
                skipLetters(innerChar);
                if (i != j) {
                    result += innerChar;
                    if (j < n - 1) {
                        if (!(i == n - 1 && j == n - 2)) {
                            result += 'U';
                        }
                    }
                }
            }
            result += ')';
            if (i != n - 1) { result += ' '; }
        }
    }
    else if (option == "all-but-one") {
        for (size_t i = 0; i < n; ++i) {
            char currentChar = 'A' + i;
            skipLetters(currentChar);
            result += currentChar;
            result += "\\(";
            for (size_t j = 0; j < n; ++j) {
                char innerChar = 'A' + j;
                skipLetters(innerChar);
                result += innerChar;
                if (j < n - 1) { result += 'I'; }
            }
            result += ')';
            if (i != n - 1) { result += ' '; }
        }
    }
    return result;
}
// End set operations

void KmerIndex::BuildReconstructionGraph(const ProgramOptions& opt) {
    // Read in FASTAs and output each color into a new separate temp file
    std::vector<std::string> transfasta = opt.transfasta;
    std::vector<std::string> tmp_files;
    std::vector<std::ofstream*> ofs; // Store pointers to circumvent certain compiler bugs where ofstream is non-movable
    gzFile fp = 0;
    kseq_t* seq;
    int l = 0;
    num_trans = 0;
    size_t range_discard = 0;
    uint32_t rb = std::max(opt.distinguish_range_begin, 0); // range begin filter
    uint32_t re = opt.distinguish_range_end == 0 ? rb : std::max(opt.distinguish_range_end, 0); // range end filter
    if (rb == 0 && re == 0) re = std::numeric_limits<uint32_t>::max();
    for (auto& fasta : transfasta) {
        std::cerr << "[build] Preparing FASTA file: " << fasta << std::endl;
        fp = transfasta.size() == 1 && transfasta[0] == "-" ? gzdopen(fileno(stdin), "r") : gzopen(fasta.c_str(), "r");
        seq = kseq_init(fp);
        while (true) {
            l = kseq_read(seq);
            if (l <= 0) {
                break;
            }
            std::string name = seq->name.s;
            std::string str = seq->seq.s;
            int color;
            try {
                color = std::stoi(name);
            }
            catch (std::exception const& e) {
                std::cerr << "Error: Non-numerical name found in " << fasta << std::endl;
                exit(1);
            }
            if (color < 0 || color > 4096) {
                std::cerr << "Error: Invalid number name, " << std::to_string(color) << ", found in " << fasta << std::endl;
                exit(1);
            }
            if (color >= tmp_files.size()) {
                tmp_files.resize(color + 1);
                ofs.resize(color + 1);
                for (int i = 0; i < tmp_files.size(); i++) {
                    if (tmp_files[i].empty()) {
                        tmp_files[i] = generate_tmp_file(opt.distinguish_output_fasta + fasta + std::to_string(i), opt.tmp_dir);
                        ofs[i] = new std::ofstream(tmp_files[i]); // Store pointers to circumvent certain compiler bugs where ofstream is non-movable
                    }
                }
            }
            if (str.length() < k) {
                continue;
            }
            if (str.length() >= rb && str.length() <= re) {
                *(ofs[color]) << ">" << std::to_string(color) << "\n" << str << "\n";
                num_trans++;
            }
            else {
                range_discard++;
            }
        }
        gzclose(fp);
        fp = 0;
    }
    for (auto& of : ofs) (*of).close(); // Close files now that we've outputted everything
    for (auto& of : ofs) delete of; // Free pointer memory

    if (range_discard > 0) {
        std::cerr << "[build] Number of input sequences filtered out due to length: " << range_discard << std::endl;
    }
    BuildDistinguishingGraph(opt, tmp_files, true); // TODO: modify to handle temporary files
}

void KmerIndex::BuildDistinguishingGraph(const ProgramOptions& opt, const std::vector<std::string>& transfasta, bool reconstruct) {
    k = opt.k;
    std::cerr << "[build] k-mer length: " << k << std::endl;
    size_t ncolors = 0;
    std::string out_file = opt.distinguish_output_fasta;
    std::vector<std::string> tmp_files;
    std::mutex visited_mutex;
    if (reconstruct) {
        tmp_files = transfasta;
    }
    else {
        for (auto& fasta : transfasta) {
            std::cerr << "[build] loading fasta file " << fasta
                << std::endl;
            tmp_files.push_back(generate_tmp_file(opt.distinguish_output_fasta + fasta, opt.tmp_dir));
        }
        std::vector<std::ofstream*> ofs; // Store pointers to circumvent certain compiler bugs where ofstream is non-movable
        for (auto tmp_file : tmp_files) ofs.push_back(new std::ofstream(tmp_file));
        num_trans = 0;

        // read fasta file using kseq
        gzFile fp = 0;
        kseq_t* seq;
        int l = 0;
        std::mt19937 gen(42);
        int countNonNucl = 0;
        int countTrim = 0;
        int countUNuc = 0;

        int i = 0;
        for (auto& fasta : transfasta) {
            if (opt.kmer_multiplicity[i] != 1) {
                i++;
                continue; // Don't do any processing or anything, use input file as-is
            }
            fp = transfasta.size() == 1 && transfasta[0] == "-" ? gzdopen(fileno(stdin), "r") : gzopen(fasta.c_str(), "r");
            seq = kseq_init(fp);
            while (true) {
                l = kseq_read(seq);
                if (l <= 0) {
                    break;
                }
                int trimNonNuclStart = 0;
                int trimNonNuclEnd = 0;
                int runningValidNuclLength = 0;
                bool finishTrimStart = false;
                std::string str = seq->seq.s;
                auto n = str.size();
                for (auto i = 0; i < n; i++) {
                    char c = str[i];
                    c = ::toupper(c);
                    if (c == 'U') {
                        str[i] = 'T';
                        countUNuc++;
                    }
                    if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'U') {
                        countNonNucl++;
                        if (runningValidNuclLength >= k && finishTrimStart) trimNonNuclEnd = i - 1; // We trim the sequence end from the last valid k-mer onward
                        runningValidNuclLength = 0;
                    }
                    else { // Valid nucleotide
                        runningValidNuclLength++;
                        if (!finishTrimStart) {
                            if (runningValidNuclLength >= k) {
                                finishTrimStart = true;
                                trimNonNuclStart = i + 1 - k; // We trim the sequence beginning until we encounter k valid nucleotides (first valid k-mer)
                            }
                        }
                        if (runningValidNuclLength >= k && finishTrimStart) trimNonNuclEnd = i; // We trim the sequence end from the last valid k-mer onward
                    }
                }
                std::transform(str.begin(), str.end(), str.begin(), ::toupper);
                str = (trimNonNuclEnd == 0 ? str.substr(trimNonNuclStart) : str.substr(trimNonNuclStart, trimNonNuclEnd + 1 - trimNonNuclStart));
                countTrim += n - str.length();
                if (str.length() >= k) {
                    *(ofs[i]) << ">" << num_trans++ << "\n" << str << "\n";
                }
                //target_lens_.push_back(seq->seq.l);
                //std::string name(seq->name.s);
                //size_t p = name.find(' ');
                //if (p != std::string::npos) {
                //  name = name.substr(0,p);
                //}
            }

            gzclose(fp);
            fp = 0;
            i++;
        }

        for (auto& of : ofs) (*of).close(); // Close files now that we've outputted everything
        for (auto& of : ofs) delete of; // Free pointer memory
        ofs.clear();

        if (countNonNucl > 0) {
            std::cerr << "[build] warning: counted " << countNonNucl << " non-ACGUT characters in the input sequence" << std::endl;
        }

        if (countTrim > 0) {
            std::cerr << "[build] warning: trimmed " << countTrim << " characters from ends of input sequences" << std::endl;
        }

        if (countUNuc > 0) {
            std::cerr << "[build] warning: replaced " << countUNuc << " U characters with Ts" << std::endl;
        }
    }

    CCDBG_Build_opt c_opt;
    c_opt.k = k;
    c_opt.nb_threads = opt.threads;
    c_opt.build = true;
    c_opt.clipTips = false;
    c_opt.deleteIsolated = false;
    c_opt.verbose = opt.verbose;
    for (int i = 0; i < tmp_files.size(); i++) {
        if (reconstruct || opt.kmer_multiplicity[i] == 1) {
            c_opt.filename_ref_in.push_back(tmp_files[i]);
        }
        else {
            c_opt.filename_seq_in.push_back(transfasta[i]);
        }
    }

    if (opt.g > 0) { // If minimizer length supplied, override the default
        c_opt.g = opt.g;
    }
    else { // Define minimizer length defaults
        int g = k - 8;
        if (k <= 13) {
            g = k - 2;
        }
        else if (k <= 17) {
            g = k - 4;
        }
        else if (k <= 19) {
            g = k - 6;
        }
        c_opt.g = g;
    }

    std::cerr << "[build] Building colored graph" << std::endl;
    ColoredCDBG<void> ccdbg = ColoredCDBG<void>(k, c_opt.g);
    ccdbg.buildGraph(c_opt);
    ccdbg.buildColors(c_opt);
    auto color_names = ccdbg.getColorNames();
    std::vector<int> color_map;
    color_map.resize(tmp_files.size());
    if (!reconstruct) {
        for (int i = 0; i < color_names.size(); i++) {
            auto color_name = color_names[i];
            for (int j = 0; j < tmp_files.size(); j++) {
                std::string fname = tmp_files[j];
                std::string real_fname = transfasta[j];
                if (opt.kmer_multiplicity[j] != 1) {
                    fname = real_fname;
                }
                if (color_name == fname) {
                    std::cerr << "        " << real_fname << ": " << j << " (multiplicity: " << opt.kmer_multiplicity[j] << ")" << std::endl;
                    color_map[i] = j;
                }
            }
        }
    }
    else {
        for (int i = 0; i < color_names.size(); i++) {
            auto color_name = color_names[i];
            for (int j = 0; j < tmp_files.size(); j++) {
                if (color_name == tmp_files[j]) {
                    color_map[i] = j;
                }
            }
        }
    }

    std::cerr << "[build] Extracting k-mers from graph" << std::endl;
    std::streambuf* buf = nullptr;
    std::ofstream of;
    if (!opt.stream_out) {
        of.open(out_file); // Write color contigs into another file
        buf = of.rdbuf();
    }
    else {
        buf = std::cout.rdbuf();
    }
    std::ostream o(buf);
    size_t max_threads_read = opt.threads;
    std::vector<std::vector<std::pair<const UnitigColors*, const UnitigMap<DataAccessor<void>, DataStorage<void>, false> > > > unitigs_v(max_threads_read);
    size_t n = 0;
    const size_t thresh_size = 50000; // Max number of unitigs across all threads
    std::mutex mutex_unitigs; // Lock for multithreading writing output FASTA file
    std::vector<std::thread> workers; // Worker threads
    uint32_t rb = std::max(opt.distinguish_range_begin, 0); // range begin filter
    uint32_t re = opt.distinguish_range_end == 0 ? rb : std::max(opt.distinguish_range_end, 0); // range end filter
    if (rb == 0 && re == 0) re = std::numeric_limits<uint32_t>::max();
    int range_discard = 0;
    int num_written = 0;
    std::string input_str = opt.input_set_operations;
    bool perform_set_operations = false;
    auto expressions = split(input_str, ' ');
    std::map<std::string, int> expr_to_int;
    if (!input_str.empty()) {
        perform_set_operations = true;
        ExpressionParser parser(input_str);  // create parser instance
        auto tokens = parser.tokenize(input_str);
        for (int i = 0; i < expressions.size(); ++i) { expr_to_int[expressions[i]] = i; }
        int maxWidth = 0;;
        for (const auto& pair : expr_to_int) { maxWidth = std::max(maxWidth, static_cast<int>(pair.first.length())); }
        for (const auto& pair : expr_to_int) { std::cout << std::setw(maxWidth + 8) << pair.first << ": " << pair.second << "\n"; }
    }
    // for bubble
    std::ofstream const_left_file;
    std::ofstream const_right_file;
    std::vector<std::ofstream*> var_files;
    if (opt.bubble) {
        for (auto f : opt.bubble_variation_output_fasta) {
            var_files.push_back(new std::ofstream(f));
            if (!(*(var_files[var_files.size() - 1])).is_open()) {
                std::cerr << "Warning: Error opening variation output files." << std::endl;
                return;
            }
        }
        // Leave out -L and -R for --flanking
        if (!opt.flanking_bubble) {
            const_left_file.open(opt.bubble_left_output_fasta);
            const_right_file.open(opt.bubble_right_output_fasta);
            if (!const_left_file.is_open() || !const_right_file.is_open()) {
                std::cerr << "Warning: Error opening anchor output files." << std::endl;
                return;
            }
        }
    }
    // TODO: Reconstruct below
    std::unordered_set<std::string> bubble_visited_anchor;
    std::unordered_set<Kmer, KmerHash> all_nodes;
    for (const auto& unitig : ccdbg) { // Iterate through all the unitigs in the de bruijn graph
        //std::cout << ":--" << unitig.referenceUnitigToString() << std::endl;
        const UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
        const UnitigMap<DataAccessor<void>, DataStorage<void>, false> unitig_ = unitig;
        unitigs_v[n % unitigs_v.size()].push_back(std::make_pair(uc, unitig_)); // unitigs_v = vector of vectors of unitigs (b/c each thread contains a vector of unitigs)
        n++;
        if (unitigs_v[unitigs_v.size() - 1].size() >= thresh_size || n >= ccdbg.size()) { // If we're ready to start processing the unitigs using a series of workers
            for (size_t u_i = 0; u_i < unitigs_v.size(); u_i++) { // u_i = specific batch of unitigs that a worker will act on
                workers.emplace_back(
                    [&, u_i] { // [&] = capture all variables by reference; u_i = specific batch of unitigs that a worker will act on
                        std::ostringstream oss;
                        std::ostringstream const_left_stream, const_right_stream; // for bubble    
                        std::vector<std::string> variation_stream; // Since stringstreams are not copyable, just use an ordinary string here
                        variation_stream.resize(var_files.size()); // Resize it to be the number of colors
                        int _num_written = 0;
                        int _range_discard = 0;
                        for (auto unitig_x : unitigs_v[u_i]) { // Go through unitigs in the batch labeled u_i
                            auto uc = unitig_x.first;
                            auto& unitig = (unitig_x.second);
                            UnitigColors::const_iterator it_uc = uc->begin(unitig);
                            UnitigColors::const_iterator it_uc_end = uc->end();
                            std::map<int, std::set<int>> k_map; // key = color; value = list of positions (i.e. k-mers) along the current unitig (note: a k-mer is a position along a unitig)
                            std::unordered_set<int> superset_colors;
                            std::map<std::string, std::set<int>> colorsetmap; // Key = Canonical k-mer in string form; Value = Set of colors; Use for the "combinations" workflow
                            for (; it_uc != it_uc_end; ++it_uc) {
                                superset_colors.insert(it_uc.getColorID());
                                int color = color_map[it_uc.getColorID()];
                                k_map[color].insert(it_uc.getKmerPosition());
                                if (opt.distinguish_combinations) colorsetmap[unitig.getUnitigKmer(it_uc.getKmerPosition()).rep().toString()].insert(color);
                                // DEBUG:
                                // std::cout << color << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).rep().toString() << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).toString() << " " << it_uc.getKmerPosition() << " " << unitig.strand << std::endl;
                            }
                            if (opt.distinguish_combinations) {
                                for (auto elem = colorsetmap.begin(); elem != colorsetmap.end(); elem++) {
                                    oss << ">";
                                    std::string color_key = "";
                                    for (auto color : elem->second) {
                                        color_key += std::to_string(color) + "_"; // Set of colors concatenated by underscores
                                    }
                                    if (color_key.size() == 0) color_key = "NULL "; // Should never happen...
                                    color_key.pop_back(); // Remove final underscore
                                    oss << color_key;
                                    oss << "\n";
                                    oss << elem->first << "\n";
                                }
                            }
                            std::string to_append = "";
                            // begin --extend 
                            if (opt.extend) {
                                int color = *superset_colors.begin();
                                std::unordered_set<std::string> visited;
                                Extend traversal = extendUnitig(unitig, visited, color, "", superset_colors, k);
                                to_append = traversal.result;
                            }
                            // end --extend
                            // begin --bubble
                            if (opt.bubble) {
                                // default --bubble
                                if (!opt.flanking_bubble) {
                                    int color = *superset_colors.begin();
                                    //std::cout << ":-" << superset_colors.size() << "-" << unitig.referenceUnitigToString() << std::endl;
                                    Bubble result = exploreBubble(ccdbg, unitig, bubble_visited_anchor, all_nodes, superset_colors, color, k, const_left_stream, variation_stream, const_right_stream, visited_mutex);
                                }
                                // end default --bubble
                                // begin --flanking
                                else {
                                    if (superset_colors.size() > 1) continue; // We only check whether single-colored unitigs are flanked by bi-colored unitigs
                                    int current_color = *superset_colors.begin(); // this is the color of the variation

                                    auto successors = unitig.getSuccessors();
                                    auto predecessors = unitig.getPredecessors();
                                    if (successors.begin() == successors.end() && predecessors.begin() == predecessors.end()) continue; // If the unitig has no successors or predecessors, skip it
                                    
                                    bool all_successors_bi_colored = true;   // Check if all successors are bi-colored
                                    int valid_length = 40;
                                    bool all_successors_valid_length = true; // Check if all successors have length >= 40 (want >= k at least)
                                    bool check_next_successor_color = true;  // Check if the next successor does not have the same color as the variation (current color)
                                    for (const auto& next : successors) {
                                        UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
                                        UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
                                        std::unordered_set<int> colors;
                                        for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
                                        if (colors.size() <= 1) {
                                            all_successors_bi_colored = false;
                                            break;
                                        }
                                        if (next.size < valid_length) {
                                            all_successors_valid_length = false;
                                            //break;                                        
                                        }
                                        // for each next, we also want to check that its successor does not have the same color as variation (current color)
                                        for (const auto& next_next: next.getSuccessors()) {
											UnitigColors::const_iterator it_next_next = next_next.getData()->getUnitigColors(next_next)->begin(next_next);
											UnitigColors::const_iterator it_next_next_end = next_next.getData()->getUnitigColors(next_next)->end();
											std::unordered_set<int> next_colors;
											for (; it_next_next != it_next_next_end; ++it_next_next) { next_colors.insert(it_next_next.getColorID()); }
                                            // check if the previously defined color is not in the set
                                            if (next_colors.find(current_color) == next_colors.end()) { // color = color of variation
                                                check_next_successor_color = false;
                                                break;
                                            }
										}
                                    }
                                    bool all_predecessors_bi_colored = true;
                                    bool all_predecessors_valid_length = true;
                                    bool check_prev_predecessor_color = true;
                                    for (const auto& prev : predecessors) {
                                        UnitigColors::const_iterator it_prev = prev.getData()->getUnitigColors(prev)->begin(prev);
                                        UnitigColors::const_iterator it_prev_end = prev.getData()->getUnitigColors(prev)->end();
                                        std::unordered_set<int> colors;
                                        for (; it_prev != it_prev_end; ++it_prev) { colors.insert(it_prev.getColorID()); }
                                        if (colors.size() <= 1) {
                                            all_predecessors_bi_colored = false;
                                            break;
                                        }
                                        if (prev.size < valid_length) {
                                        	all_predecessors_valid_length = false;
											//break; 
                                        }
                                        // for each prev, we also want to check that its predecessor does not have the same color as variation (current color)
                                        for (const auto& prev_prev : prev.getPredecessors()) {
                                            UnitigColors::const_iterator it_prev_prev = prev_prev.getData()->getUnitigColors(prev_prev)->begin(prev_prev);
                                            UnitigColors::const_iterator it_prev_prev_end = prev_prev.getData()->getUnitigColors(prev_prev)->end();
                                            std::unordered_set<int> prev_colors;
                                            for (; it_prev_prev != it_prev_prev_end; ++it_prev_prev) { prev_colors.insert(it_prev_prev.getColorID()); }
                                            // check if the previously defined color is not in the set
                                            if (prev_colors.find(current_color) == prev_colors.end()) { // color = color of variation
                                                check_prev_predecessor_color = false;
                                                break;
                                            }
                                        }
                                    }
                                    // Output the current unitig if all successors and predecessors are bi-colored
                                    if (all_successors_bi_colored && all_predecessors_bi_colored &&
                                        //all_successors_valid_length && all_predecessors_valid_length &&
                                        check_next_successor_color && check_prev_predecessor_color) {
                                        variation_stream[current_color] += ">" + std::to_string(current_color) + "\n" + unitig.referenceUnitigToString() + "\n";
                                    }
                                }
                                // end --flanking
                            }
                            // end --bubble
                            std::set<int> positions_to_remove; // Positions (i.e. k-mers) along the current unitig that will be cut out
                            std::map<std::vector<int>, int> result_map; // Key = colors; Value = Position (i.e. k-mer)
                            std::stringstream ss; // For --combinations outputting aggregated colors
                            int int_to_print = 1;
                            if (perform_set_operations) {
                                for (const auto& expr : expressions) {
                                    ExpressionParser expr_parser(expr);
                                    Node* root = expr_parser.parse();
                                    std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                    for (const auto& k_elem : k_map) {
                                        int curr_pos = -1;
                                        auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                        std::string colored_contig = "";
                                        for (const auto& pos : k_elem.second) {
                                            if (set_operation_result.count(pos)) {
                                                std::string km = unitig.getUnitigKmer(pos).toString();
                                                if (curr_pos == -1) { colored_contig = km; }
                                                else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                else {
                                                    if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                    else _range_discard++;
                                                    colored_contig = km;
                                                }
                                                curr_pos = pos;
                                            }
                                        }
                                        if (colored_contig != "") {
                                            if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                            else _range_discard++;
                                        }
                                    }
                                }
                                continue; // colored contigs already extracted, continue to next iteration
                            }
                            { // Special operations (e.g. --all, --all-but-one, --combinations, or reconstruct)
                                if (!opt.distinguish_all_but_one_color && !opt.distinguish_combinations) { // Workflow: Find k-mers unique/exclusive to each color, or reconstruct, or --all
                                    std::string default_str = stringOpt(tmp_files.size(), "default"); // for n=8, generate "A\(AIBI...IH) B\(AIBI...IH) ... H\(AIBI...IH)"
                                    auto all_expr = split(default_str, ' ');
                                    for (const auto& expr : all_expr) {
                                        ExpressionParser expr_parser(expr);
                                        Node* root;
					std::set<int> set_operation_result;
					if (opt.distinguish_union) reconstruct = true; // --all
					if (!reconstruct) {
					    root = expr_parser.parse();
					    set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
					}
                                        // write out what remains among the contigs
                                        for (const auto& k_elem : k_map) {
                                            int curr_pos = -1;
                                            auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                            std::string colored_contig = "";
                                            for (const auto& pos : k_elem.second) {
                                                if (reconstruct || set_operation_result.count(pos)) {
                                                    std::string km = unitig.getUnitigKmer(pos).toString();
                                                    if (curr_pos == -1) { colored_contig = km; }
                                                    else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                    else {
                                                        if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                        else _range_discard++;
                                                        colored_contig = km;
                                                    }
                                                    curr_pos = pos;
                                                }
                                            }
                                            if (colored_contig != "") {
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                else _range_discard++;
                                            }
                                        }
                                    }
                                    continue; // colored contigs already extracted, continue to next iteration
                                }
                                else if (!opt.distinguish_union && !opt.distinguish_combinations) { // workflow: opt.distinguish_all_but_one_color (e.g. if we have 8 colors, for each color, output all k-mers except those that are 8-colored)
                                    std::string all_but_one_str = stringOpt(tmp_files.size(), "all-but-one"); // generate "A\(AIBI...IH) B\(AIBI...IH) ... H\(AIBI...IH)"
                                    auto all_expr = split(all_but_one_str, ' ');
                                    for (const auto& expr : all_expr) {
                                        ExpressionParser expr_parser(expr);
                                        Node* root = expr_parser.parse();
                                        std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                        // write out what remains among the contigs
                                        for (const auto& k_elem : k_map) {
                                            int curr_pos = -1;
                                            auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                            std::string colored_contig = "";
                                            for (const auto& pos : k_elem.second) {
                                                if (set_operation_result.count(pos)) {
                                                    std::string km = unitig.getUnitigKmer(pos).toString();
                                                    if (curr_pos == -1) { colored_contig = km; }
                                                    else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                    else {
                                                        if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                        else _range_discard++;
                                                        colored_contig = km;
                                                    }
                                                    curr_pos = pos;
                                                }
                                            }
                                            if (colored_contig != "") {
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                else _range_discard++;
                                            }
                                        }
                                    }
                                    continue; // colored contigs already extracted, continue to next iteration
                                }
                                else if (!opt.distinguish_all_but_one_color) { // workflow: opt.distinguish_combinations (output every combination [i.e. each area in the Venn diagram], with colors separated by underscores)
                                    std::set<int> exclusive_positions;
                                    for (const auto& k_elem : k_map) { // Iterate over each color and get exclusive positions for that color
                                        std::set<int> current_positions = k_elem.second;
                                        if (current_positions.size() == 1) { // Check if the k-mer is exclusive to one color
                                            // This k-mer is unique to one color, so process it accordingly
                                            exclusive_positions.insert(current_positions.begin(), current_positions.end());
                                            std::vector<int> single_color_key = { k_elem.first };
                                            int integer_value = *(current_positions.begin());
                                            result_map[single_color_key] = integer_value; // Add this single-color k-mer to the result map
                                        }
                                        else {
                                            bool is_exclusive = true;
                                            for (const auto& other_elem : k_map) { // Check if the positions of this color appear in other colors
                                                if (k_elem.first != other_elem.first) {
                                                    std::set<int> intersection;
                                                    std::set_intersection(current_positions.begin(), current_positions.end(),
                                                        other_elem.second.begin(), other_elem.second.end(),
                                                        std::inserter(intersection, intersection.begin()));
                                                    if (!intersection.empty()) {
                                                        is_exclusive = false;
                                                        break;
                                                    }
                                                }
                                            }
                                            if (is_exclusive) { exclusive_positions.insert(current_positions.begin(), current_positions.end()); }
                                        }
                                    }
                                    positions_to_remove.insert(exclusive_positions.begin(), exclusive_positions.end()); // Add the exclusive positions to the positions_to_remove
                                    // updated color key method below
                                   /*
                                    if (k_map.size() >= 2) {
                                        std::vector<int> consolidated_key;
                                        for (auto it = k_map.begin(); it != k_map.end(); ++it) { // Iterate over k_map and build the consolidated key
                                            consolidated_key.push_back(it->first);
                                        }
                                        int integer_value = *(k_map.begin()->second.begin());
                                        result_map[consolidated_key] = integer_value; // {0 1 2} : 0 (color : position)
                                    }
                                    */
                                    std::vector<int> consolidated_key;
                                    for (auto it = k_map.begin(); it != k_map.end(); ++it) {
                                        consolidated_key.push_back(it->first);
                                    }
                                    int integer_value = *(k_map.begin()->second.begin());
                                    result_map[consolidated_key] = integer_value; // {0} or {0 1 2} : 0 (color : position)

                                    for (const auto& k_elem : result_map) {
                                        int curr_pos = -1;
                                        std::string colored_contig;
                                        auto color = k_elem.first;
                                        auto pos = k_elem.second;
                                        for (int i = 0; i < color.size(); ++i) {
                                            ss << color[i];
                                            if (i < color.size() - 1) {
                                                ss << "_";
                                            }
                                        }
                                    }
                                }
                            } // end if !opt.distinguish_union
                            // Now, write out what remains among the contigs
                            for (const auto& k_elem : k_map) {
                                if (opt.distinguish_combinations) break; // Don't need this loop; we're just doing combinations which we have already stored in output stream
                                int curr_pos = -1;
                                std::string colored_contig = "";
                                auto color = k_elem.first;
                                std::string color_key = std::to_string(color);
                                //std::string contig_metadata = " :" + unitig.dist + "," + unitig.len + "," + unitig.size + "," + unitig.strand;
                                if (opt.distinguish_combinations) { color_key = ss.str(); int_to_print++; }
                                if (int_to_print == k_map.size()) {
                                    for (const auto& pos : k_elem.second) {
                                        if (!positions_to_remove.count(pos)) {
                                            std::string km = unitig.getUnitigKmer(pos).toString();
                                            if (curr_pos == -1) { // How to correspond color?
                                                colored_contig = km;
                                            }
                                            else if (pos == curr_pos + 1) {
                                                colored_contig += km[km.length() - 1];
                                            }
                                            else {
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re) {
                                                    oss << ">" << color_key << "\n" << colored_contig + to_append << "\n"; _num_written++;
                                                }
                                                else _range_discard++;
                                                colored_contig = km;
                                            }
                                            curr_pos = pos;
                                        }
                                    }
                                    if (colored_contig != "") {
                                        if (colored_contig.length() >= rb && colored_contig.length() <= re) {
                                            // if opt.extend, then we need to append the traversal.result to the colored_contig
                                            oss << ">" << color_key << "\n" << colored_contig + to_append << "\n"; _num_written++;
                                        }
                                        else _range_discard++;
                                    }
                                }
                            }
                        }
                        if (var_files.size() == 0) {
                            std::unique_lock<std::mutex> lock(mutex_unitigs);
                            o << oss.str();
                        }
                        // write bubble results to respective files
                        if (var_files.size() > 0) {
                            std::unique_lock<std::mutex> lock(mutex_unitigs);
                            const_left_file << const_left_stream.str();
                            for (int i = 0; i < var_files.size(); i++) {
                                *(var_files[i]) << variation_stream[i]; // Output variation within bubbles
                            }
                            const_right_file << const_right_stream.str();
                        }
                        num_written += _num_written;
                        range_discard += _range_discard;
                    }
                );
            }
            for (auto& t : workers) t.join();
            workers.clear();
            for (size_t u_i = 0; u_i < unitigs_v.size(); u_i++) {
                unitigs_v[u_i].clear();
            }
        }
    }
    o.flush();
    if (opt.bubble) {
        const_left_file.flush();
        const_right_file.flush(); // for bubble
        const_left_file.close();
        const_right_file.close();
        // Close bubble files
        for (auto& of : var_files) {
            (*of).flush();
            (*of).close();
        } // Close files now that we've outputted everything
        for (auto& of : var_files) delete of; // Free pointer memory
    }
    ccdbg.clear(); // Free memory associated with the colored compact dBG
    ncolors = tmp_files.size(); // Record the number of "colors"
    for (auto tmp_file : tmp_files) std::remove(tmp_file.c_str()); // Remove temp files needed to make colored graph
    if (!opt.stream_out && !opt.bubble) {
        of.close();
    }
    tmp_files.clear();
    if (opt.verbose) {
        if (range_discard > 0) std::cerr << "[build] Number of output sequences filtered out due to length: " << range_discard << std::endl;
        std::cerr << "[build] Number of output sequences written: " << num_written << std::endl;
    }
}
