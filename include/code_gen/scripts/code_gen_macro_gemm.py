from math import *
import numpy as np

def code_gen_macro_gemm(radix=3):
    M_PI = 3.141592653589793
    macro_gemm = ''''''
    for i in range(radix):
        for j in range(radix):
            macro_gemm += f'''
#define radix3_a{i}{j}_x {cos(- 2.0 * M_PI * (float)(i * j) / float(radix))}f;
#define radix3_a{i}{j}_y {sin(- 2.0 * M_PI * (float)(i * j) / float(radix))}f;
'''
    macro_gemm += f'''

#define GEMM_radix{radix} (b1, b2, b3, c1, c2, c3) \\
'''
    for i in range(radix):
        for j in range(radix):
            macro_gemm += f'''c{i}.x += radix{radix}_a{i}{j}_x * b{j}.x - radix{radix}_a{i}{j}_y * b{j}.y; c{i}.y = radix{radix}_a{i}{j}_y * b{j}.x + radix{radix}_a{i}{j}_x * b{j}.y;\    
'''
    print(macro_gemm)
    
    
if __name__ == "__main__":
    code_gen_macro_gemm(3)