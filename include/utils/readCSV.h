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
    std::string &file_name, bool if_print=false){
    // std::ifstream &file, bool if_print=false){
    std::vector<std::vector<long long int>> data;
    std::ifstream file(file_name);
    
    if (file.is_open()) {
        // std::cout << "File opened successfully." << std::endl;
        std::string line;

        while (std::getline(file, line)) {
            std::vector<long long int> row = split(line, ',');
            if(if_print)
            printf("%s\n", line.c_str());
            data.push_back(row);
        }
    } else {
        std::cout << "Failed to open file." << std::endl;
    }    
    return data;
}
}