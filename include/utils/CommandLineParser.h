#include <iostream>
#include <cstdlib> 
#include <cstring> 
class ProgramConfig {
public:
    long long N = 0;
    long long logN = 0;
    long long bs = 0;
    long long bs_end = 0;
    long long bs_gap = 1;
    bool if_profile = false;
    bool if_verify = false;
    bool if_bench = true;
    bool if_ft = false;
    bool if_err = false;
    int datatype = 0;
    int thread_bs = 1;
    std::string param_file_path = "../include/param/A100/param_float2.csv";

    static void displayHelp() {
        std::cout << "Usage: program [options]\n"
                  << "Options:\n"
                  << "  --logN <value>       Set logN to <value>, which determines N as 2^<value>.\n"
                  << "  --bs <value>         Set block size to <value>.\n"
                  << "  --bs_end <value>     Set block size end to <value> (for iterative tests).\n"
                  << "  --bs_gap <value>     Set block size gap to <value> (for iterative tests).\n"
                  << "  --if_profile <0|1>   Enable (1) or disable (0) profiling.\n"
                  << "  --if_verify <0|1>    Enable (1) or disable (0) verification.\n"
                  << "  --if_bench <0|1>     Enable (1) or disable (0) benchmarking.\n"
                  << "  --if_ft <0|1>        Enable (1) or disable (0) fault tolerance.\n"
                  << "  --if_err <0|1>       Enable (1) or disable (0) error injection.\n"
                  << "  --datatype <type>    0 for FP32, 1 for FP64.\n"
                  << "  --thread_bs <value>  Set batches per block to <value>.\n"
                  << "  -h, --help               Display this help message and exit.\n";
    }
    
    bool parseParameter(const char* parameterName, const char* value) {
        if (strcmp(parameterName, "--logN") == 0) {
            logN = std::atol(value);
            N = 1LL << std::atol(value);
        } else if (strcmp(parameterName, "--bs") == 0) {
            bs = std::atol(value);
        } else if (strcmp(parameterName, "--bs_end") == 0) {
            bs_end = std::atol(value);
        } else if (strcmp(parameterName, "--bs_gap") == 0) {
            bs_gap = std::atol(value);
        } else if (strcmp(parameterName, "--if_profile") == 0) {
            if_profile = std::atol(value);
        } else if (strcmp(parameterName, "--if_verify") == 0) {
            if_verify = std::atol(value);
        } else if (strcmp(parameterName, "--if_ft") == 0) {
            if_ft = std::atol(value);
        } else if (strcmp(parameterName, "--if_err") == 0) {
            if_err = std::atol(value);
            if(if_err != 0) if_ft = 1;
        } else if (strcmp(parameterName, "--if_bench") == 0) {
            if_bench = std::atol(value);
        } else if (strcmp(parameterName, "--datatype") == 0) {
            datatype = std::atoi(value);
            if (datatype == 1) param_file_path = "../include/param/A100/param_double2.csv";
        } else if (strcmp(parameterName, "--thread_bs") == 0) {
            thread_bs = std::atoi(value);
        } else {
            return false; // 未知的参数名
        }
        return true; // 参数解析成功
    }

    
    void parseCommandLine(int argc, char* argv[]) {
        for (int i = 1; i < argc; i += 1){
            if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
                displayHelp();
                exit(EXIT_FAILURE); // Return false to indicate the program should exit after showing help
            }
        }
        for (int i = 1; i < argc; i += 2) {
            if (i + 1 >= argc || !parseParameter(argv[i], argv[i + 1])) {
                std::cerr << "Invalid or incomplete argument: " << argv[i] << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    void displayConfig() const {
        std::cout << "N: " << N << ", bs: " << bs << ", bs_end: " << bs_end << ", bs_gap: " << bs_gap << std::endl;
        std::cout << "if_profile: " << if_profile << ", if_verify: " << if_verify << ", if_bench: " << if_bench << std::endl;
        std::cout << "if_ft: " << if_ft << ", if_err: " << if_err << std::endl;
        std::cout << "datatype: " << datatype << ", thread_bs: " << thread_bs << std::endl;
    }
};