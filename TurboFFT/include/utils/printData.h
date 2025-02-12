#include "utils.h"
namespace utils{
template <typename DataType>
void printData(DataType* res, long long int N){
    double rel_error = 0.;
    for(int i = 0; i < N; ++i){
        printf("res[%d] = %f + %f j\n", i, res[i].x, res[i].y);
    }
}

}