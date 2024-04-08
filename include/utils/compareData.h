#include "utils.h"
namespace utils{
    
template<typename DataType>
void compareData(DataType* res, DataType *res_ref, long long int N, 
                double error_bound){
    printf("Compare begin.\n");
    double rel_error = 0.;
    for(int i = 0; i < N; ++i){
        rel_error = abs((res[i].x - res_ref[i].x) / (1 + abs(res_ref[i].x)));
        // if(i % 10000 == 0){
        //     printf("res[%d].x=%f, res_ref[%d].x=%f, rel_error=%f\n", 
        //     i, res[i].x, i, res_ref[i].x, rel_error);
        // }
        if(rel_error > error_bound || res[i].x != res[i].x){
        // if(true){
            printf("Error detected: res[%d].x=%f, res_ref[%d].x=%f, rel_error=%f\n", 
            i, res[i].x, i, res_ref[i].x, rel_error);
            return;
        }
        rel_error = abs((res[i].y - res_ref[i].y) / (1 + abs(res_ref[i].y)));
        // if(i % 10000 == 0){
        //     printf("res[%d].y=%f, res_ref[%d].y=%f, rel_error=%f\n", 
        //     i, res[i].y, i, res_ref[i].y, rel_error);
        // }
        if(rel_error > error_bound || res[i].y != res[i].y){
        // if(true){
            // if(i % 10000 == 0){
            printf("Error detected: res[%d].y=%f, res_ref[%d].y=%f, rel_error=%f\n", 
            i, res[i].y, i, res_ref[i].y, rel_error);
            return;
        }
    }
    printf("Compare finished.\n");
}

}
