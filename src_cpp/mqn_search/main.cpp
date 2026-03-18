#include <iostream>
#include <fstream>
#include <algorithm>

#include "reader.h"

typedef unsigned int uint;

using namespace std;

static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

int main(int argc, char *argv[])
{
    std::string in_path = argv[1];
    std::string out_path = argv[2];

    std::ifstream in_file(in_path);
    std::ofstream out_file(out_path);

    if( !in_file.is_open() || !out_file.is_open() ){
        std::cout << "Problem with file I/O." << std::endl;
        return 1;
    }

    reader rd;

    std::string line;

    while( std::getline( in_file, line ) ){

        rtrim(line);

        std::size_t npos = line.find('\t');
        std::string SMILES = line.substr(0,npos);

        rd.read( SMILES );

        out_file << SMILES << ';';

        const std::vector<uint>& mqn = rd.get_mqn();

        for( uint i = 0; i <= 41; ++i ){ // mqn_sum to fluorine_count
            out_file << mqn[i] << ';';
        }
        out_file << mqn[42] << '\n'; // phosphorous_count


        rd.cleanup();
    }
    return 0;
}
