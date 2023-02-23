rm bin/ess_color_kmerlist
rm bin/ess_color_compress
rm bin/ess_color_decompress

g++  -g -std=c++11  -DNDEBUG  src/kmerlist.cpp -o bin/ess_color_kmerlist
g++  -g -std=c++11  -DNDEBUG  src/compress.cpp -o bin/ess_color_compress -lcmph -lsdsl -ldivsufsort -ldivsufsort64  -lpthread  -Wno-unused-result -Wno-format
g++  -g -std=c++11  -DNDEBUG  src/decompress.cpp -o bin/ess_color_decompress -lcmph -lsdsl -ldivsufsort -ldivsufsort64  -lpthread  -Wno-unused-result -Wno-format

