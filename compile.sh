cd include/code_gen/scripts
python fft_codegen.py --gpu T4 --datatype float2
# python fft_codegen.py --gpu T4 --datatype float2  --if_ft 1
# python fft_codegen.py --gpu T4 --datatype float2 --if_ft 1 --if_err_injection 1
# python fft_codegen.py --gpu T4 --datatype double2
# python fft_codegen.py --gpu T4 --datatype double2  --if_ft 1
# python fft_codegen.py --gpu T4 --datatype double2 --if_ft 1 --if_err_injection 1
# cd -
# cd build
# make -j
