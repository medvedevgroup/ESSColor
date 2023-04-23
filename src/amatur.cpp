#include <iostream>

#include "lib/sshash/src/common.hpp"
#include "lib/sshash/src/bench_utils.hpp"
#include "lib/sshash/src/check_utils.hpp"
#include "lib/sshash/src/build.cpp"
#include "lib/sshash/src/query.cpp"
#include "lib/sshash/src/permute.cpp"


// using namespace sshash;

int main(int argc, char** argv) {
    std::string input_filename = "mega.essd";
    std::cout<<"hellp"<<std::endl;
    auto k = 31;
    auto m = k/2;

    dictionary dict;

    build_configuration build_config;
    build_config.k = k;
    build_config.m = m;

    build_config.canonical_parsing = true;
    build_config.verbose = true;
    build_config.print();

    dict.build(input_filename, build_config);
    assert(dict.k() == k);
    std::cout<<"dict built complete";

    //     // auto output_filename = "amatur_out";
    //     // // essentials::logger("saving data structure to disk...");
    //     // // essentials::save(dict, output_filename.c_str());
    //     // // essentials::logger("DONE");
    // }

    // return 0;
}