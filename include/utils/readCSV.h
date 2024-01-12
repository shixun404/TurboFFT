#include "utils.h"
namespace utils{
std::vector<long long int> split(const std::string& s, char delimiter) {
    std::vector<long long int> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(std::stoi(token));
    }
    return tokens;
}

std::vector<std::vector<long long int>> load_parameters(
    std::ifstream &file){
    std::vector<std::vector<long long int>> data;
    std::string line;

    while (std::getline(file, line)) {
        std::vector<long long int> row = split(line, ',');
        printf("%s\n", line.c_str());
        data.push_back(row);
    }
    return data;
}
}