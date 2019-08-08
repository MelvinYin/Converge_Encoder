#include <iostream>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <math.h>
#include <assert.h>


std::vector<std::string> read_file(std::string const &fileName) {
  std::vector<std::string> vecOfStrs;
  // Open the File
  std::ifstream ifs;
  ifs.open(fileName);
  std::string line;
  while (std::getline(ifs, line))
  {
    // Line contains string of length > 0 then save it in vector
    if (!line.empty()) {
      vecOfStrs.push_back(line);
    }
  }
  ifs.close();
  return vecOfStrs;
}


void load_fasta_sequences(const std::string& filename,
  std::vector<std::string>& headers, std::vector<std::vector<int>>& sequences) {
  std::set<char> kAlphabtets_set = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
                                    'K',
                                    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                                    'W', 'Y'};
  
  std::array<char, 20> kAlphabtets = {'A', 'C', 'D', 'E', 'F', 'G',
                                      'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                                      'R', 'S', 'T', 'V', 'W', 'Y'};
  std::map<char, int> letter_int_map;
  for (int i=0;i<20;i++){
    letter_int_map[kAlphabtets[i]] = i;
  }
  
  std::vector<std::string> vecOfStrs = read_file(filename);
  std::string current_header;
  std::vector<int> current_seq;
  
  for (std::string &line: vecOfStrs) {
    if (line[0] == '>') {
//      Header line
      if (current_seq.empty()) {
        current_header = line;
      } else {
        current_seq.shrink_to_fit();
        current_header.shrink_to_fit();
        sequences.push_back(current_seq);
        headers.push_back(current_header);
        current_header = line;
        current_seq.clear();
      }
    } else {
//      Sequence line
      for (char letter: line) {
        bool is_in = kAlphabtets_set.find(letter) != kAlphabtets_set.end();
        if (is_in) {
          current_seq.push_back(letter_int_map.at(letter));
        }
      }
    }
  }
  if (!current_seq.empty()) {
    sequences.push_back(current_seq);
    headers.push_back(current_header);
  }
  sequences.shrink_to_fit();
  headers.shrink_to_fit();
  assert (sequences.size() == headers.size());
}


void split_seq(std::vector<int> &seq, int denom, int length,
  std::vector<std::vector<int>> seed_seqs) {
  int num_full_cycles = (((int) seq.size())-length) / denom;
  for (int i = 0; i < num_full_cycles; i++) {
    auto first = seq.begin() + i * denom;
    auto last = seq.begin() + i * denom + length;
    std::vector<int> sub_vector(first, last);
    seed_seqs.push_back(sub_vector);
  }
  if (seq.size() > num_full_cycles * denom) {
    auto first = seq.end() - length;
    auto last = seq.end();
    std::vector<int> sub_vector(first, last);
    seed_seqs.push_back(sub_vector);
  }
}


std::vector<std::vector<int>> load_seed_seq(const std::string& filename) {
  std::vector<std::string> headers;
  std::vector<std::vector<int>> seed_seq_rawsplit;
  load_fasta_sequences(filename, headers, seed_seq_rawsplit);
  std::vector<std::vector<int>> seed_seqs;
  for (std::vector<int> seed_seq_wronglen: seed_seq_rawsplit) {
    std::vector<std::vector<int>> seed_seq_split;
    split_seq(seed_seq_wronglen, 10, 30, seed_seq_split);
    for (std::vector<int> seed_seq: seed_seq_split) {
      seed_seq.shrink_to_fit();
      seed_seqs.push_back(seed_seq);
    }
  }
  return seed_seqs;
}


std::vector<std::vector<double>> read_blosum(const std::string &filename) {
  std::vector<std::vector<double>> kBlosum(20, std::vector<double>(20));
  
  std::array<char, 20> kAlphabets = {'A', 'C', 'D', 'E', 'F', 'G',
                                     'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                                     'R', 'S', 'T', 'V', 'W', 'Y'};
  std::unordered_set<char> kAlphabtets_set = {'A', 'C', 'D', 'E', 'F', 'G',
                                              'H', 'I', 'K', 'L', 'M', 'N',
                                              'P', 'Q', 'R', 'S', 'T', 'V',
                                              'W', 'Y'};
  
  std::map<char, int> letter_int_map;
  for (int i=0;i<20;i++){
    letter_int_map[kAlphabets[i]] = i;
  }
  
  std::vector<std::string> vecOfStrs = read_file(filename);
  std::vector<char> local_alphabets;
  
  int row_i = 0;
  bool is_in;
  for (const std::string &line: vecOfStrs) {
    if (line[0] == '#') {
      continue;
    }
    if (line[0] == '*') {
      break;
    }
    if (line[0] == ' ') {
      for (int i=3;i<3*23+3;i+=3) {
        assert(line[i] != ' ');
        local_alphabets.push_back(line[i]);
      };
      continue;
    } else {
      char current_row_char = local_alphabets[row_i];
      is_in = kAlphabtets_set.find(current_row_char) != kAlphabtets_set.end();
      if (!is_in) {
        continue;
      }
      for (int i=2; i<3*23+2; i+=3) {
        bool is_negative = false;
        if (line[i] == '-') {
          is_negative = true;
        }
        auto digit = (double) (((int) line[i + 1]) - '0');
        
        if (is_negative and line[i + 1] != '0') {
          digit *= -1;
        }
        char current_col_char = local_alphabets[i/3];
        is_in = kAlphabtets_set.find(current_col_char) != kAlphabtets_set.end();
        if (is_in) {
          int blosum_row = letter_int_map.at(current_row_char);
          int blosum_col = letter_int_map.at(current_col_char);
          kBlosum[blosum_row][blosum_col] = digit;
        }
      }
      row_i++;
    }
  }
  return kBlosum;
}


template <typename... Args>
void save(const std::string& filename, const Args&... saves){
  std::ofstream file;
  file.open(filename, std::ios_base::binary);
  {
    cereal::BinaryOutputArchive oarchive(file); // Create an output archive
    oarchive(saves...);
  }
  file.close();
}


template <typename... Args>
void load(const std::string& filename, Args&... outputs){
  std::ifstream file;
  file.open(filename, std::ios_base::binary);
  {
    cereal::BinaryInputArchive iarchive(file); // Create an output archive
    iarchive(outputs...);
  }
}


int main(){
  //  Encode proteome into vector<pair<string, string>>
  // save fasta seq and names separately.
  std::string proteome_filename = "proteome.fasta";
  std::vector<std::string> headers;
  std::vector<std::vector<int>> sequences;
  load_fasta_sequences(proteome_filename, headers, sequences);
  save("proteome_binary", headers, sequences);

// Encode seed into vector<string>, each of spacing 10
  std::string seed_filename = "initial.fasta";
  std::vector<std::vector<int>> seed_seqs = load_seed_seq(seed_filename);
  save("seed_seq_binary", seed_seqs);
  
// Encode blosum
  std::string blosum_filename = "BLOSUM62";
  std::vector<std::vector<double>> kBlosum = read_blosum(blosum_filename);
  save("blosum_binary", kBlosum);
}
