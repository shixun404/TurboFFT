export CUDA_SAMPLES_COMMON=$(pwd)/Common
gpu=A100
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
cmake ..
make -j


echo "Benchmarking 1/24"
./turbofft --if_bench 1  --thread_bs 1 --gpu A100 --datatype 0 --if_thread_ft 0 --if_ft  0  --if_err 0 > Benchmark=1_TurboFFT_FP32.csv
echo "Benchmarking 2/24"
./turbofft --if_bench 1  --thread_bs 1 --gpu A100 --datatype 0 --if_thread_ft 1 --if_ft  0  --if_err 0 > Benchmark=1_TurboFFT_FP32_FT_Thread.csv
echo "Benchmarking 3/24"
./turbofft --if_bench 1  --thread_bs 1 --gpu A100 --datatype 0 --if_thread_ft 0 --if_ft  1  --if_err 0 > Benchmark=1_TurboFFT_FP32_FT_BS=1.csv
echo "Benchmarking 4/24"
./turbofft --if_bench 1  --thread_bs 2 --gpu A100 --datatype 0 --if_thread_ft 0 --if_ft  1  --if_err 0 > Benchmark=1_TurboFFT_FP32_FT_BS=2.csv
echo "Benchmarking 5/24"
./turbofft --if_bench 1  --thread_bs 4 --gpu A100 --datatype 0 --if_thread_ft 0 --if_ft  1  --if_err 0 > Benchmark=1_TurboFFT_FP32_FT_BS=4.csv

echo "Benchmarking 6/24"
./turbofft --if_bench 2  --thread_bs 1 --gpu A100 --datatype 0 --if_thread_ft 0 --if_ft  0  --if_err 0 > Benchmark=2_TurboFFT_FP32.csv
echo "Benchmarking 7/24"
./turbofft --if_bench 2  --thread_bs 4 --gpu A100 --datatype 0 --if_thread_ft 0 --if_ft  1  --if_err 0 > Benchmark=2_TurboFFT_FP32_FT_BS=1.csv
echo "Benchmarking 8/24"
./turbofft --if_bench 2  --thread_bs 4 --gpu A100 --datatype 0 --if_thread_ft 0 --if_ft  1  --if_err 1 > Benchmark=2_TurboFFT_FP32_ERR_BS=32.csv


echo "Benchmarking 9/24"
./turbofft --if_bench 1  --thread_bs 1 --gpu A100 --datatype 1 --if_thread_ft 0 --if_ft  0  --if_err 0 > Benchmark=1_TurboFFT_FP64.csv
echo "Benchmarking 10/24"
./turbofft --if_bench 1  --thread_bs 1 --gpu A100 --datatype 1 --if_thread_ft 1 --if_ft  0  --if_err 0 > Benchmark=1_TurboFFT_FP64_FT_Thread.csv
echo "Benchmarking 11/24"
./turbofft --if_bench 1  --thread_bs 1 --gpu A100 --datatype 1 --if_thread_ft 0 --if_ft  1  --if_err 0 > Benchmark=1_TurboFFT_FP64_FT_BS=1.csv
echo "Benchmarking 12/24"
./turbofft --if_bench 1  --thread_bs 2 --gpu A100 --datatype 1 --if_thread_ft 0 --if_ft  1  --if_err 0 > Benchmark=1_TurboFFT_FP64_FT_BS=2.csv
echo "Benchmarking 13/24"
./turbofft --if_bench 1  --thread_bs 4 --gpu A100 --datatype 1 --if_thread_ft 0 --if_ft  1  --if_err 0 > Benchmark=1_TurboFFT_FP64_FT_BS=4.csv

echo "Benchmarking 14/24"
./turbofft --if_bench 2  --thread_bs 1 --gpu A100 --datatype 1 --if_thread_ft 0 --if_ft  0  --if_err 0 > Benchmark=2_TurboFFT_FP64.csv
echo "Benchmarking 15/24"
./turbofft --if_bench 2  --thread_bs 4 --gpu A100 --datatype 1 --if_thread_ft 0 --if_ft  1  --if_err 0 > Benchmark=2_TurboFFT_FP64_FT_BS=1.csv
echo "Benchmarking 16/24"
./turbofft --if_bench 2  --thread_bs 4 --gpu A100 --datatype 1 --if_thread_ft 0 --if_ft  1  --if_err 1 > Benchmark=2_TurboFFT_FP64_ERR_BS=32.csv

echo "Benchmarking 17/24"
./turbofft --if_bench 11  --gpu A100 --datatype 0 > Benchmark=1_cuFFT_FP32.csv
echo "Benchmarking 18/24"
./turbofft --if_bench 11  --gpu A100 --datatype 1 > Benchmark=1_cuFFT_FP64.csv
echo "Benchmarking 19/24"
./turbofft --if_bench 12  --gpu A100 --datatype 0 > Benchmark=2_cuFFT_FP32.csv
echo "Benchmarking 20/24"
./turbofft --if_bench 12  --gpu A100 --datatype 1 > Benchmark=2_cuFFT_FP64.csv
echo "Benchmarking 21/24"
./turbofft --if_bench 21  --gpu A100 --datatype 0 > Benchmark=1_cuFFT_FP32_FT_OFFLINE.csv
echo "Benchmarking 22/24"
./turbofft --if_bench 21  --gpu A100 --datatype 1 > Benchmark=1_cuFFT_FP64_FT_OFFLINE.csv
echo "Benchmarking 23/24"
./turbofft --if_bench 32  --gpu A100 --datatype 0 > Benchmark=2_Offline_FP32_ERR.csv
echo "Benchmarking 24/24"
./turbofft --if_bench 32  --gpu A100 --datatype 1 > Benchmark=2_Offline_FP64_ERR.csv
echo "Benchmarking Finished"


cp Benchmark=1_cuFFT_FP32_FT_OFFLINE.csv ../artifact_data/cuFFT_data
cp Benchmark=1_cuFFT_FP32.csv  ../artifact_data/cuFFT_data
cp Benchmark=1_cuFFT_FP64_FT_OFFLINE.csv  ../artifact_data/cuFFT_data
cp Benchmark=1_cuFFT_FP64.csv  ../artifact_data/cuFFT_data
cp Benchmark=2_cuFFT_FP32.csv  ../artifact_data/cuFFT_data
cp Benchmark=2_cuFFT_FP64.csv  ../artifact_data/cuFFT_data
cp Benchmark=2_Offline_FP32_ERR.csv  ../artifact_data/cuFFT_data
cp Benchmark=2_Offline_FP64_ERR.csv ../artifact_data/cuFFT_data

cp Benchmark=1_TurboFFT_FP32_FT_BS=1.csv ../artifact_data/TurboFFT_data
cp Benchmark=1_TurboFFT_FP32_FT_BS=2.csv ../artifact_data/TurboFFT_data
cp Benchmark=1_TurboFFT_FP32_FT_BS=4.csv ../artifact_data/TurboFFT_data
cp Benchmark=1_TurboFFT_FP32_FT_Thread.csv ../artifact_data/TurboFFT_data
cp Benchmark=1_TurboFFT_FP32.csv ../artifact_data/TurboFFT_data
cp Benchmark=1_TurboFFT_FP64_FT_BS=1.csv ../artifact_data/TurboFFT_data
cp Benchmark=1_TurboFFT_FP64_FT_BS=2.csv ../artifact_data/TurboFFT_data
cp Benchmark=1_TurboFFT_FP64_FT_BS=4.csv ../artifact_data/TurboFFT_data
cp Benchmark=1_TurboFFT_FP64_FT_Thread.csv ../artifact_data/TurboFFT_data
cp Benchmark=1_TurboFFT_FP64.csv ../artifact_data/TurboFFT_data
cp Benchmark=2_TurboFFT_FP32_ERR_BS=32.csv  ../artifact_data/TurboFFT_data
cp Benchmark=2_TurboFFT_FP32_FT_BS=1.csv ../artifact_data/TurboFFT_data
cp Benchmark=2_TurboFFT_FP32.csv  ../artifact_data/TurboFFT_data
cp Benchmark=2_TurboFFT_FP64_ERR_BS=32.csv  ../artifact_data/TurboFFT_data
cp Benchmark=2_TurboFFT_FP64_FT_BS=1.csv ../artifact_data/TurboFFT_data
cp Benchmark=2_TurboFFT_FP64.csv ../artifact_data/TurboFFT_data

cd ../../
git clone https://github.com/DTolm/VkFFT.git
cp CMakeLists_VkFFT.txt VkFFT/CMakeLists.txt
cd VkFFT
mkdir build
cd build
cmake -DVKFFT_BACKEND=1 ..
make -j
cp ../../test_VkFFT.py .
python test_VkFFT.py
cp vkFFT.pt ../../TurboFFT/artifact_data/VkFFT_data/vkFFT.pt



cd ../../TurboFFT/plot_scripts
python3 plot.py
python3 plot_fig1.py
python3 plot_fig10.py
python3 plot_fig11.py
python3 plot_fig12.py
python3 plot_fig13.py
python3 plot_fig14.py
python3 plot_fig16.py
python3 plot_fig17.py
python3 plot_fig18.py
python3 plot_fig19.py
python3 plot_fig20.py