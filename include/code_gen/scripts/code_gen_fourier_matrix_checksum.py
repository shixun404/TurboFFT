import torch as th
from math import cos, sin, pi
def code_gen_fourier_matrix_checksum(radix=2):
    v = th.ones(radix, dtype=th.complex64)
    w = th.ones((radix, radix), dtype=th.complex64)
    a = cos(-2 * pi / radix) + sin(-2 * pi / radix) * 1.j
    a_3 = cos(-2 * pi / 3) + sin(-2 * pi / 3) * 1.j
    sum = 1.0 + 0.j
    
    for i in range(radix):    
        for j in range(radix):
            w[i, j] = cos(-2 * pi * i * j / radix) + sin(-2 * pi * i * j / radix) * 1.j
            
    for i in range(0, radix):
        v[i] = cos(-2 * pi * (i % 3) / 3) + sin(-2 * pi * (i % 3) / 3) * 1.j
    print(v)
    # for i in range(3):
    #     print(i, v[i].item().real, v[i].item().imag)
    # assert 0
    # for i in range(radix):
    #     print(w[i])
    # print(w)
    # print(th.matmul(w, v))
    checksum = th.matmul(w, v)
    for i in range(radix):
        print(i, checksum[i].item().real, checksum[i].item().imag)
    # print(checksum)
    assert 0
    for i in range(radix):
        print(f"#define A_radix{radix}_{i}_x {checksum[i].item().real}f")
        print(f"#define A_radix{radix}_{i}_y {checksum[i].item().imag}f")

if __name__ == "__main__":
    code_gen_fourier_matrix_checksum(256)
    # print(cos(pi))
    
            
    