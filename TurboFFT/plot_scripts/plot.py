import torch as th
TurboFFT = th.ones(2, 30, 30)
TurboFFT *= 1e6

filename = f'../artifact_data/TurboFFT_data/Benchmark=1_TurboFFT_FP32.csv'
try:
    f = open(filename, 'r')
except:
    print('Cannot open file:', filename)
lines = f.readlines()
for line in lines:
    if line.split(',')[0] != 'turboFFT':
        continue
    logN = int(line.split(',')[1])
    logBS = int(line.split(',')[2])
    exec_time = float(line.split(',')[3])
    TurboFFT[0][logN][logBS] = min(TurboFFT[0][logN][logBS], exec_time)
filename = '../artifact_data/TurboFFT_data/Benchmark=1_TurboFFT_FP64.csv'
f = open(filename, 'r')
lines = f.readlines()
for line in lines:
    if line.split(',')[0] != 'turboFFT':
        continue
    logN = int(line.split(',')[1])
    logBS = int(line.split(',')[2])
    exec_time = float(line.split(',')[3])
    TurboFFT[1][logN][logBS] = exec_time

TurboFFT[TurboFFT == 1e6] = 1

th.save(TurboFFT, '../artifact_data/TurboFFT.pt')


filenames = ['../artifact_data/cuFFT_data/Benchmark=1_cuFFT_FP32.csv',
            '../artifact_data/cuFFT_data/Benchmark=1_cuFFT_FP64.csv']
cuFFT = th.ones(2, 30, 30)
for p in range(2):
    filename = filenames[p]
    f = open(filename, 'r')
    lines_cufft = f.readlines()
    for i in range(len(lines_cufft)):
        if lines_cufft[i].split(',')[0] != 'cuFFT':
            continue
        cuFFT[p][int(lines_cufft[i].split(',')[1])][int(lines_cufft[i].split(',')[2])] = float(lines_cufft[i].split(',')[3])
        
th.save(cuFFT,'../artifact_data/cuFFT.pt')


TurboFFT_FT = th.zeros(2, 5, 30, 30)

# p = 0
for p in range(2):
    precision = 32 + 32 * p
    filename_list = [f'../artifact_data/cuFFT_data/Benchmark=1_cuFFT_FP{precision}_FT_OFFLINE.csv', f'../artifact_data/TurboFFT_data/Benchmark=1_TurboFFT_FP{precision}_FT_Thread.csv'] + [f'../artifact_data/TurboFFT_data/Benchmark=1_TurboFFT_FP{precision}_FT_BS={2 ** i}.csv' for 
                    i in range(3)]
    # f = f'..\data_ft\FP32_A100_FT_THREAD.csv'
    # TurboFFT_FT[p][0][logN][logBS] = exec_time
    for i, f in enumerate(filename_list):
        f = open(f, 'r')
        lines = f.readlines()
        for line in lines:
            splitted = line.split(',')
            # if splitted[0] != 'turboFFT' and splitted[0] != 'cuFFT':
            #     continue
            logN = int(splitted[1])
            logBS = int(splitted[2])
            exec_time = float(splitted[3])
            # if  i >= 1 and exec_time <= 0.5 * TurboFFT_FT[p][0][logN][logBS]:
            #     exec_time = TurboFFT_FT[p][1][logN][logBS]
            # if TurboFFT_FT[p][i][logN][logBS] == 0:
            TurboFFT_FT[p][i][logN][logBS] = exec_time
            # else:
            #     TurboFFT_FT[p][i][logN][logBS] = min(exec_time, TurboFFT_FT[p][i][logN][logBS])
th.save(TurboFFT_FT, '../artifact_data/TurboFFT_FT.pt')
    