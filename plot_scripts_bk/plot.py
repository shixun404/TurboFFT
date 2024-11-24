import torch as th
TurboFFT = th.ones(2, 30, 30)
TurboFFT *= 1e6

filename = f'../build/Benchmark=1_TurboFFT_FP32.csv'
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
filename = '../build/Benchmark=1_TurboFFT_FP64.csv'
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

th.save(TurboFFT, '../artifact_data/TurboFFT_bk.pt')


filenames = ['../artifact_data/Benchmark=1_cuFFT_FP32.csv',
            '../artifact_data/Benchmark=1_cuFFT_FP64.csv']
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