export CUDA_SAMPLES_COMMON=$(pwd)/Common
gpu=T4
cd TurboFFT
cd include/code_gen/scripts
python fft_codegen.py --gpu $gpu --datatype float2
python fft_codegen.py --gpu $gpu --datatype float2 --if_thread_ft 1
python fft_codegen.py --gpu $gpu --datatype float2  --if_ft 1
python fft_codegen.py --gpu $gpu --datatype float2 --if_ft 1 --if_err_injection 1
python fft_codegen.py --gpu $gpu --datatype double2
python fft_codegen.py --gpu $gpu --datatype double2 --if_thread_ft 1
python fft_codegen.py --gpu $gpu --datatype double2  --if_ft 1
python fft_codegen.py --gpu $gpu --datatype double2 --if_ft 1 --if_err_injection 1
cd -
mkdir build
cd build
cmake .. -DARCH_SM=75
make -j


echo "Benchmarking 1/7"
./turbofft --if_bench 1  --thread_bs 1 --gpu T4 --datatype 0 --if_thread_ft 0 --if_ft  0  --if_err 0 > Benchmark=1_TurboFFT_FP32.csv
echo "Benchmarking 2/7"
./turbofft --if_bench 2  --thread_bs 1 --gpu T4 --datatype 0 --if_thread_ft 0 --if_ft  0  --if_err 0 > Benchmark=2_TurboFFT_FP32.csv
echo "Benchmarking 3/7"
./turbofft --if_bench 2  --thread_bs 4 --gpu T4 --datatype 0 --if_thread_ft 0 --if_ft  1  --if_err 0 > Benchmark=2_TurboFFT_FP32_FT_BS=1.csv
echo "Benchmarking 4/7"
./turbofft --if_bench 2  --thread_bs 4 --gpu T4 --datatype 0 --if_thread_ft 0 --if_ft  1  --if_err 1 > Benchmark=2_TurboFFT_FP32_ERR_BS=32.csv

echo "Benchmarking 5/7"
./turbofft --if_bench 11  --gpu T4 --datatype 0 > Benchmark=1_cuFFT_FP32.csv

echo "Benchmarking 6/7"
./turbofft --if_bench 12  --gpu T4 --datatype 0 > Benchmark=2_cuFFT_FP32.csv

echo "Benchmarking 7/7"
./turbofft --if_bench 32  --gpu T4 --datatype 0 > Benchmark=2_Offline_FP32_ERR.csv

echo "Benchmarking Finished"



cp Benchmark=1_cuFFT_FP32.csv  ../artifact_data/cuFFT_data
cp Benchmark=2_cuFFT_FP32.csv  ../artifact_data/cuFFT_data
cp Benchmark=2_Offline_FP32_ERR.csv  ../artifact_data/cuFFT_data


cp Benchmark=1_TurboFFT_FP32.csv ../artifact_data/TurboFFT_data
cp Benchmark=2_TurboFFT_FP32_ERR_BS=32.csv  ../artifact_data/TurboFFT_data
cp Benchmark=2_TurboFFT_FP32_FT_BS=1.csv ../artifact_data/TurboFFT_data
cp Benchmark=2_TurboFFT_FP32.csv  ../artifact_data/TurboFFT_data


cd ../../
git clone https://github.com/DTolm/VkFFT.git
cp CMakeLists_VkFFT.txt VkFFT/CMakeLists.txt
cd VkFFT
mkdir build
cd build
cmake -DVKFFT_BACKEND=1 ..
make -j

cp ../../test_VkFFT_T4.py .
python test_VkFFT_T4.py
cp vkFFT.pt ../../TurboFFT/artifact_data/VkFFT_data/vkFFT_T4.pt


cd ../../TurboFFT/plot_scripts
python3 plot_fig21.py
python3 plot_fig22.py