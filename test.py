template_function  = ''
for i in range(30):
    for j in range(3):
        template_function += f"tempate<typename DataType> void __global__ fft_radix_2_logN_{i}_dim_{j}(DataType *, DataType *, DataType *, DataType*, int, int);\n"

print(template_function)
func_entry_float = ""
func_entry_double = ""
float_bd = [14, 23]
double_bd = [14, 23]
def helper_include(dtype, if_ft, if_err_inj, bd):
    func_entry = f'''
template<> struct TurboFFT_Kernel_Entry<{dtype}, {if_ft}, {if_err_inj}>
'''
    name = f"void (*turboFFTArr [26][3])({dtype} *, {dtype} *, {dtype} *, {dtype} *, int, int) ="
    func_entry += '{\n' + name + '{\n {NULL, NULL, NULL},\n'
    for i in range(1, 26):
        func_0 = f"fft_radix_2<{dtype}, {i}, 0, {if_ft}, {if_err_inj}>"
        func_1 = f"fft_radix_2<{dtype}, {i}, 1, {if_ft}, {if_err_inj}>" if i >= bd[0] else "NULL"
        func_2 = f"fft_radix_2<{dtype}, {i}, 2, {if_ft}, {if_err_inj}>" if i >= bd[1] else "NULL"
        func_entry += "{" + func_0 + ", " +  func_1 + ", " + func_2 + "},\n"
    func_entry += '''
};
};
'''
    print(func_entry)

helper_include("float2", if_ft=0, if_err_inj=0, bd=float_bd)
helper_include("float2", if_ft=1, if_err_inj=0, bd=float_bd)
helper_include("float2", if_ft=1, if_err_inj=1, bd=float_bd)


helper_include("double2", if_ft=0, if_err_inj=0, bd=double_bd)
helper_include("double2", if_ft=1, if_err_inj=0, bd=double_bd)
helper_include("double2", if_ft=1, if_err_inj=1, bd=double_bd)