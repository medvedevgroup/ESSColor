#pragma once

#include "../dictionary.hpp"
#include "../util.hpp"

#include "../gz/zip_stream.hpp"
#include "streaming_query_canonical_parsing.hpp"
#include "streaming_query_regular_parsing.hpp"

namespace sshash {

template <typename Query>
void streaming_query_from_fasta_file(dictionary const* dict, std::istream& is) {
    std::string line;
    uint64_t k = dict->k();
    Query query(dict);
    while (!is.eof()) {
        query.start();
        std::getline(is, line);  // skip first header line
        std::getline(is, line);
        if (line.size() < k) continue;
        for (uint64_t i = 0; i != line.size() - k + 1; ++i) {
            char const* kmer = line.data() + i;
            auto answer = query.lookup_advanced(kmer);
            if( answer.kmer_id != constants::invalid_uint64){
                std::cout<<kmer<<" "<<answer.kmer_id<<std::endl;
            }else{
                std::cout<<kmer<<" "<<-1<<std::endl;
            }
        }
    }
}


void dictionary::streaming_query_from_file(std::string const& filename) const {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    streaming_query_from_fasta_file<streaming_query_canonical_parsing>(this, is);
    is.close();
    return;
}

}  // namespace sshash