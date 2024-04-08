cd include/code_gen/scripts
python fft_codegen.py --datatype float2 
python fft_codegen.py --datatype float2  --if_ft 1
python fft_codegen.py --datatype float2 --if_ft 1 --if_err_injection 1
python fft_codegen.py --datatype double2
python fft_codegen.py --datatype double2  --if_ft 1
python fft_codegen.py --datatype double2 --if_ft 1 --if_err_injection 1

# python fft_codegen.py --datatype float2 #--if_ft 1 --if_err_injection 1
# python fft_codegen.py --datatype float2 --if_ft 1 --if_err_injection 1


# python fft_codegen.py --datatype  double2  --if_ft 1 # --if_err_injection 1
# python fft_codegen.py --datatype  double2  --if_ft 1  --if_err_injection 1
cd -

# cd build_no_ft_fp32
# cd build_ft_fp32
# cd build_err_fp32

cd build

# cd build_err_fp64
# cd build_no_ft_fp64
# cd build_no_ft_fp64
# cd build_cufft_fp64
make -j
