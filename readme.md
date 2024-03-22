# TurboFFT

## Run the following command to generate FFT kernel and compile
```
mkdir build
cd build
cmake ..
cd ..
./compile.sh
```
## Run the following command to profile/verify/benchmark
turboFFT takes in 5 arguments, which are `logN`, `batch_size`, `if_profile`, `if_verify`, and `if_bench`.

### Profile a $16 \times 2^{10}$ FFT computation
```
cd build
./turbofft 10 16 1 0 0 
```
### Verify a $16 \times 2^{10}$ FFT computation
```
cd build
./turbofft 10 16 0 1 0 
```
### Benchmark turboFFT
```
cd build
./turbofft 0 0 0 0 1 
```

## Set your kernel parameters!

|logN | Kernel Launch Times | logN1 |logN2 |logN3 | BS1| BS2| BS3 | n1| n2| n3|
|--|--|-|--|--|--|--|--|--|--|--|
|7| 1| 7| 0| 0| 1| 0| 0| 8| 0| 0|
|8| 1| 8| 0| 0| 1| 0| 0| 16| 0| 0|
|9| 1| 9| 0| 0| 1| 0| 0| 8| 0| 0|
|10| 1| 10| 0| 0| 1| 0| 0| 8| 0| 0|
|11| 1| 11| 0| 0| 1| 0| 0| 8| 0| 0|
|12| 1| 12| 0| 0| 1| 0| 0| 16| 0| 0|
|13 |1 |13| 0| 0| 1| 0| 0| 32| 0| 0|
|14 |2 |7 |7| 0| 16| 16| 0| 8| 8| 0|
|15 |2 |7 |8| 0| 16| 8| 0| 8| 8| 0|
|23| 3| 7| 8| 8| 32| 16| 16|  16| 16| 16
|24| 3| 8| 8| 8| 16| 16| 16|  16| 16| 16

Some parameters combination may not work, for example you should make sure $BS_i \times logN_i \leq max\_shared\_memory\_size$
1. You can config your own kernel parameters by edit the parameter file in `include/param/param_A100_double2.csv`
2. Run the following command to regenerate the kernel code and compile again
```
./compile.sh
```

## Fault Tolerance mode and Error Injection Mode

1. Modify the code generator `fft_codegen.py` arguments

| Parameter        | Type  | Default Value | Description                                              |
|------------------|-------|---------------|----------------------------------------------------------|
| `if_ft`          | int   | 0             | Flag to indicate fault detection (0 for False, 1 for True) |
| `if_err_injection` | int   | 0             | Flag to indicate error injection (0 for False, 1 for True) |
| `err_smoothing`  | int   | 1000          | Error smoothing parameter                                |
| `err_inj`        | int   | 100           | Error injection value                                |
| `err_threshold`  | float | 1e-3          | Error threshold                                          |
| `datatype`       | str   | "double2"     | Data type with a default value of "double2"              |

1. Run the following command to 
```
cd include/code_gen/scripts
python fft_codegen.py --datatype="double2" --if_ft=1  --if_err_injection=1 
```

## We currently only support radix-2