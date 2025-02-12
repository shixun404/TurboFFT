#include "utils.h"
namespace utils{
    void getDFTMatrix(double2* dest, long long int N){
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < N; ++j){
                dest[i * N + j].x = cos(- 2 * M_PI * i * j / N);
                dest[i * N + j].y = sin(- 2 * M_PI * i * j / N);
            }
        }
    }

    void getDFTMatrix(float2* dest, long long int N){
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < N; ++j){
                dest[i * N + j].x = cosf(- (float)(2 * M_PI * i * j) / float(N));
                dest[i * N + j].y = sinf(- (float)(2 * M_PI * i * j) / float(N));
            }
        }
    }

    void getDFTMatrixChecksum(double2* dest, long long int N){
        double2 r[3];
        r[0].x = 1.0;
        r[0].y = 0.0;
        r[1].x = -0.5;
        r[1].y = -0.8660253882408142;
        r[2].x = -0.5;
        r[2].y = 0.8660253882408142;
        
        for(int i = 0; i < N; ++i){
            dest[i].x = 0;
            dest[i].y = 0;
            for(int j = 0; j < N; ++j){
                double x = cos(- (float)(2 * M_PI * i * j) / float(N));
                double y = sin(- (float)(2 * M_PI * i * j) / float(N));
                dest[i].x += x * r[j % 3].x - y * r[j % 3].y;
                dest[i].y += x * r[j % 3].y + y * r[j % 3].x;
            }
        }
    }

    void getDFTMatrixChecksum(float2* dest, long long int N){
        float2 r[3];
        r[0].x = 1.0f;
        r[0].y = 0.0f;
        r[1].x = -0.5f;
        r[1].y = -0.8660253882408142f;
        r[2].x = -0.5f;
        r[2].y = 0.8660253882408142f;
        
        for(int i = 0; i < N; ++i){
            dest[i].x = 0;
            dest[i].y = 0;
            for(int j = 0; j < N; ++j){
                float x = cosf(- (float)(2 * M_PI * i * j) / float(N));
                float y = sinf(- (float)(2 * M_PI * i * j) / float(N));
                dest[i].x += x * r[j % 3].x - y * r[j % 3].y;
                dest[i].y += x * r[j % 3].y + y * r[j % 3].x;
            }
        }
    }
}