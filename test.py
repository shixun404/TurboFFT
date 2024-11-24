
func_entry_float = ""
func_entry_double = ""
A100_float_bd = [14, 23]
A100_double_bd = [14, 23]
T4_float_bd = [13, 23]
T4_double_bd = [13, 23]

float_bd =  {"T4":T4_float_bd, "A100":A100_float_bd}
double_bd =  {"T4":T4_double_bd, "A100":A100_double_bd}

gpus = ["T4", "A100"]
spec = {"T4":75, "A100":80}
def helper_include(dtype, if_ft, if_err_inj, bd, gpu, spec, max_logN=26):
    func_entry = f'''
template<> struct TurboFFT_Kernel_Entry<{dtype}, {if_ft}, {if_err_inj}, {spec}>
'''
    name = f"void (*turboFFTArr [{max_logN}][3])({dtype} *, {dtype} *, {dtype} *, {dtype} *, int, int) ="
    func_entry += '{\n' + name + '{\n {NULL, NULL, NULL},\n'
    for i in range(1, max_logN):

        func_0 = f'''fft_radix_2<{dtype}, {i}, 0, {if_ft}, {if_err_inj}>'''
        func_1 = f'''fft_radix_2<{dtype}, {i}, 1, {if_ft}, {if_err_inj}>''' if i >= bd[0] else "NULL"
        func_2 = f'''fft_radix_2<{dtype}, {i}, 2, {if_ft}, {if_err_inj}>''' if i >= bd[1] else "NULL"
        func_entry += "{" + func_0 + ", " +  func_1 + ", " + func_2 + "},\n"
    func_entry += '''
};
};
'''
    print(func_entry)


def helper_filename(dtype, bd, max_logN=26):
    include_file = ''
    for i in range(1, max_logN):
        include_file += f'''
        #include "code_gen/generated/{dtype}/fft_radix_2_logN_{i}_upload_0.cuh"'''
        include_file += f'''
        #include "code_gen/generated/{dtype}/fft_radix_2_logN_{i}_upload_1.cuh" ''' if i >= bd[0] else ''
        include_file += f'''
        #include "code_gen/generated/{dtype}/fft_radix_2_logN_{i}_upload_2.cuh" ''' if i >= bd[1] else ''
    return include_file


for gpu in gpus:
    print(f"#if ARCH_SM == {spec[gpu]}")
    print(helper_filename('float2', bd=float_bd[gpu]))
    helper_include("float2", if_ft=0, if_err_inj=0, bd=float_bd[gpu], gpu=gpu, spec=spec[gpu])
    helper_include("float2", if_ft=1, if_err_inj=0, bd=float_bd[gpu], gpu=gpu, spec=spec[gpu])
    helper_include("float2", if_ft=1, if_err_inj=1, bd=float_bd[gpu], gpu=gpu, spec=spec[gpu])

    print(helper_filename('double2', bd=float_bd[gpu]))
    helper_include("double2", if_ft=0, if_err_inj=0, bd=double_bd[gpu], gpu=gpu, spec=spec[gpu])
    helper_include("double2", if_ft=1, if_err_inj=0, bd=double_bd[gpu], gpu=gpu, spec=spec[gpu])
    helper_include("double2", if_ft=1, if_err_inj=1, bd=double_bd[gpu], gpu=gpu, spec=spec[gpu])
    print("#endif")