gpu=$1
cd include/code_gen/scripts
python fft_codegen.py --gpu $gpu --datatype float2
python fft_codegen.py --gpu $gpu --datatype float2  --if_ft 1
python fft_codegen.py --gpu $gpu --datatype float2 --if_ft 1 --if_err_injection 1
python fft_codegen.py --gpu $gpu --datatype double2
python fft_codegen.py --gpu $gpu --datatype double2  --if_ft 1
python fft_codegen.py --gpu $gpu --datatype double2 --if_ft 1 --if_err_injection 1
cd -
cd build
make -j
./turbofft --if_bench 2  --thread_bs 1 --gpu A100 --datatype 0 --if_ft  0  --if_err 0
./turbofft --if_bench 2  --thread_bs 32 --gpu A100 --datatype 0 --if_ft  1  --if_err 0
./turbofft --if_bench 2  --thread_bs 2 --gpu A100 --datatype 0 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 4 --gpu A100 --datatype 0 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 8 --gpu A100 --datatype 0 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 16 --gpu A100 --datatype 0 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 32 --gpu A100 --datatype 0 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 1 --gpu A100 --datatype 0 --if_ft  1  --if_err 1


./turbofft --if_bench 2  --thread_bs 1 --gpu A100 --datatype 1 --if_ft  0  --if_err 0
./turbofft --if_bench 2  --thread_bs 32 --gpu A100 --datatype 1 --if_ft  1  --if_err 0
./turbofft --if_bench 2  --thread_bs 2 --gpu A100 --datatype 1 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 4 --gpu A100 --datatype 1 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 8 --gpu A100 --datatype 1 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 16 --gpu A100 --datatype 1 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 32 --gpu A100 --datatype 1 --if_ft  1  --if_err 1
./turbofft --if_bench 2  --thread_bs 1 --gpu A100 --datatype 1 --if_ft  1  --if_err 1
