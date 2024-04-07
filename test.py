template_function  = ''
for i in range(30):
    for j in range(3):
        template_function += f"tempate<typename DataType> void __global__ fft_radix_2_logN_{i}_dim_{j}(DataType *, DataType *, DataType *, DataType*, int, int);\n"

print(template_function)
func_entry_float = ""
func_entry_double = ""
float_bd = [14, 22]
double_bd = [14, 23]
for i in range(26):
    func_0 = f"fft_radix_2<float2, {i}, 0>"
    func_1 = f"fft_radix_2<float2, {i}, 1>" if i >= float_bd[0] else "NULL"
    func_2 = f"fft_radix_2<float2, {i}, 2>" if i >= float_bd[1] else "NULL"
    func_entry_float += "{" + func_0 + ", " +  func_1 + ", " + func_2 + "},\n"


for i in range(26):
    func_0 = f"fft_radix_2<double2, {i}, 0>"
    func_1 = f"fft_radix_2<double2, {i}, 1>" if i >= double_bd[0] else "NULL"
    func_2 = f"fft_radix_2<double2, {i}, 2>" if i >= double_bd[1] else "NULL"
    func_entry_double += "{" + func_0 + ", " + func_1 + ", " + func_2 + "},\n"

print(func_entry_float)
print(func_entry_double)