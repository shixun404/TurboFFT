import torch as th
from math import *
import numpy as np
class TurboFFT:
    def __init__(self, global_tensor_shape=[256, 1], radix=2, WorkerFFTSize = 8, data_type='double'):
        self.fft = ''''''
        self.exponent = int(log(N, radix))
        self.data_type = "double"
        self.gPtr = "gPtr"
        self.rPtr = "rPtr"
        self.shPtr = "shPtr"
        self.WorkerFFTSizes = WorkerFFTSizes
        self.global_tensor_shape = global_tensor_shape
        self.logWorkerFFTSize = int(log(self.WorkerFFTSize, 2))
        # self.coalesced = 128 / 64 if data_type == "float" else 128 / 128
        self.state_vec = th.zeros(self.WorkerFFTSize, int(log(self.WorkerFFTSize, 2)))
        self.data = th.rand(self.WorkerFFTSize, dtype=th.cfloat)
        self.output = th.zeros(self.WorkerFFTSize, dtype=th.cfloat)
        for i in range(32):
            for j in range(5):
                self.state_vec[i, j] = (i // (2 ** j)) % 2

        self.threadblock_tensor_shape = []
        for size, N_tmp in zip(WorkerFFTSizes, self.global_tensor_shape):
            threadblock_tensor_shape = []
            while N_tmp > 1:
                threadblock_tensor_shape.append(size if N_tmp >= size else N_tmp)
                N_tmp /= self.WorkerFFTSize
            threadblock_tensor_shape.reverse()
            self.threadblock_tensor_shape.append(threadblock_tensor_shape)

        print(self.threadblock_tensor_shape)
        assert 0

    def codegen(self,):
        exponent = int(log(N, radix))
        N__ = N
        # print(f"########################### TRANSPOSE N = 2 ** {exponent} ###########################################################")
        # print(f"N={N}, radix={radix}, N / radix = {N / radix}, WorkerFFTSize={WorkerFFTSize}")
        
        NumWorkersPerSignal = int((N // WorkerFFTSize))
        ThreadblockFFTBatchSize = int(blockdim // (N // WorkerFFTSize))
        

            

        order = []
        for i in range(WorkerFFTSize * 2):
            order.append(i)
        offset = WorkerFFTSize

        self.fft = f'''extern __shared__ {data_type} shared[];
        __global__ void __launch_bounds__({blockdim}) fft_radix{radix}_logN{exponent}''' \
        + f'''({data_type}2* inputs, {data_type}2* outputs, {data_type}2* r_2)''' + ''' {
        '''

        self.fft += self.global2reg(self.tensor_shape, BS, strdBS, strdTensor, strdTB, if_shared=False, bs_major=False)

        for i in len(tensor_shape) - 1:
            next_tensor_shape = tensor_reshape(self.tensor_shape, )
            self.fft += fft_reg()
            self.fft += tensor_reshape_shared(self.tensor_shape, next_tensor_shape)
            self.fft += shared2reg(next_tensor_shape)
            self.tensor_shape = next_tensor_shape
        self.fft += fft_reg()
        self.fft += reg2global()
        

        fft += epilogue()
        return self.fft

    def global2reg(self, bs_dim, data_dim, ):
        pass
            
    def fft_reg(self, tensor_shape, WorkerFFTSize, cur_stage):
        print(self.tensor_shape, self.tensor_shape[-1])
        bs = self.WorkerFFTSize // self.tensor_shape[-1]
        logbs = int(log(bs, 2))
        reg_tensor_stride = th.zeros(self.logWorkerFFTSize)
        strd = 1
        i = 0
        while strd < self.WorkerFFTSize:
            reg_tensor_stride[i] = strd
            i += 1
            strd *= 2
        # st = self.logWorkerFFTSize - 1 if self.WorkerFFTSize == self.tensor_shape[-1] else self.logWorkerFFTSize - 1 - logbs
        st = self.logWorkerFFTSize - 1
        

        tmp = 1
        stride = 1
        cur_stage_stride = 1
        tensor_shape_bs = [self.BS] + self.tensor_shape
        cur_stage_stride = int(self.BS * th.prod(th.as_tensor(self.tensor_shape)) / self.WorkerFFTSize)
        self.fft += f'''
        offset += tid;
        '''
        
        self.fft += '''
            __syncthreads();
        '''
        for i in range(self.WorkerFFTSize):
            self.fft += f'''
            {self.rPtr}[{i}] = {self.shPtr}[offset + {cur_stage_stride * i}];
            '''
    
        print(st, reg_tensor_stride)
        
        # for i in range(st, -1, -1):
        for i in range(st, logbs - 1, -1):
            print("*************************")
            print(reg_tensor_stride[i])
            for j in range(self.WorkerFFTSize):
                if self.state_vec[j, i] == 1:
                    continue
                id_j1 = int(th.dot(self.state_vec[j], reg_tensor_stride))
                id_j2 = int(id_j1 + reg_tensor_stride[i])
                # id_k = int(th.dot(self.state_vec[j, :i], reg_tensor_stride[:i]))
                id_k = int(th.dot(self.state_vec[j, logbs:i], reg_tensor_stride[logbs:i]))
                
                print(id_j1, id_j2, id_k, i,  (2 ** (i + 1 - logbs)), st, logbs)
                # tmp_angle = -2 * pi * id_k * 1 / (2 ** (i + 1))
                tmp_angle = (-2 * id_k * 1 / (2 ** (i + 1 - logbs))) * pi
                rel_bounds = 1e-8
                abs_bounds = 1e-5
                self.fft += f'''
                tmp = {self.rPtr}[{id_j1}] 
                {self.rPtr}[{id_j1}] = turboFFT_ZADD({self.rPtr}[{id_j1}], tmp, {self.rPtr}[{id_j2}]);
                {self.rPtr}[{id_j2}] = turboFFT_ZSUB({self.rPtr}[{id_j2}], tmp, {self.rPtr}[{id_j2}]);
                angle.x = {0 if np.allclose(0, cos(tmp_angle), rel_bounds, abs_bounds) else cos(tmp_angle)};
                angle.y = {0 if np.allclose(0, sin(tmp_angle), rel_bounds, abs_bounds) else sin(tmp_angle)};
                tmp = {self.rPtr}[{id_j2}];
                turboFFT_ZMUL({self.rPtr}[{id_j2}], tmp, angle);
                '''

                tmp = self.data[id_j1].item()
                self.data[id_j1] = tmp + self.data[id_j2]
                self.data[id_j2] = tmp - self.data[id_j2]
                # angle = cos(-2 * pi * id_k * 1 / (2 ** (i + 1))) + sin(-2 * pi * id_k * 1 / (2 ** (i + 1))) * 1.j
                angle = cos(-2 * pi * id_k * 1 / (2 ** (i + 1 - logbs))) + sin(-2 * pi * id_k * 1 / (2 ** (i + 1 - logbs))) * 1.j
                self.data[id_j2] = self.data[id_j2] * angle
        # reg_tensor_stride_reverse = th.cat((reg_tensor_stride[:(st + 1)].flip(0), reg_tensor_stride[(st + 1):]), dim=0)
        reg_tensor_stride_reverse = th.cat((reg_tensor_stride[:logbs], reg_tensor_stride[logbs:].flip(0)), dim=0)

        for i in range(self.WorkerFFTSize):
            output_id = int(th.dot(self.state_vec[i], reg_tensor_stride_reverse))
            input_id = int(th.dot(self.state_vec[i], reg_tensor_stride))
            print("output, input id:", output_id, input_id)
            self.output[output_id] = self.data[input_id]

        
        
        # threadId to tensor coordinates
        tmp = 1
        stride = 1
        cur_stage_stride = 1
        self.tensor_shape = self.tensor_shape[:cur_stage] + [self.tensor_shape[-1]] + self.tensor_shape[cur_stage:-1]
        tensor_shape_bs = [self.BS] + self.tensor_shape
        print(tensor_shape_bs)
        for i in range(len(tensor_shape_bs)):
            stride *= tensor_shape_bs[i]
            if i == cur_stage + 1:
                cur_stage_stride = int(stride / tensor_shape_bs[i])
                continue
            self.fft += f'''
            offset += ((tid / {tmp}) % {tensor_shape_bs[i]}) * {stride};
            '''
            tmp *= tensor_shape_bs[i]
        if cur_stage == len(self.tensor_shape) - 1:
            cur_stage_stride = int(self.BS * th.prod(th.as_tensor(self.tensor_shape)) / self.WorkerFFTSize)
        # assert 0
        self.fft += '''
        __syncthreads();
        '''
        dict_output = {}
        for i in range(self.WorkerFFTSize):
            output_id = int(th.dot(self.state_vec[i], reg_tensor_stride_reverse))
            input_id = int(th.dot(self.state_vec[i], reg_tensor_stride))
            dict_output[output_id] = input_id
        print(dict_output)
        for output_id in range(self.WorkerFFTSize): 
            print(output_id, dict_output[output_id])
            self.fft += f'''
            {self.shPtr}[offset + {cur_stage_stride * output_id}] = {self.rPtr}[{dict_output[output_id]}];
            '''
        print(cur_stage_stride)
        # assert 0

if __name__ == '__main__':
    global_tensor_shape = [256,1]
    fft = TurboFFT(global_tensor_shape=global_tensor_shape)
    fft.fft_reg(0)
    fft.data = fft.output.clone()
    fft.fft_reg(1)
    fft.data = fft.output.clone()
    # torch_fft_result = th.fft.fft(fft.data).flatten()
    torch_fft_result = th.fft.fft(fft.data.reshape(fft.WorkerFFTSize // 4, 4), dim=0).flatten()
    fft.fft_reg(2)
    print(fft.output)
    print(torch_fft_result)
    rel_err = th.norm(fft.output - torch_fft_result) / th.norm(torch_fft_result)
    print(fft.fft)
    print(f"Rel_ERR = {rel_err}")
    if rel_err > 1e-3:
        print("Error!\n")
    else:
        print("Passed!\n")