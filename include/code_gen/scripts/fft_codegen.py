import torch as th
from math import *
import numpy as np
class TurboFFT:
    def __init__(self, N=256, radix=2, WorkerFFTSize = 8, data_type='double'):
        self.fft = ''''''
        self.N = N
        self.BS = 4
        self.strdBS = 1
        self.strdTensor = 1
        self.strdTB = 1
        self.exponent = int(log(N, radix))
        self.data_type = "double"
        self.gPtr = "gPtr"
        self.rPtr = "rPtr"
        self.shPtr = "shPtr"
        self.WorkerFFTSize = WorkerFFTSize
        self.logWorkerFFTSize = int(log(self.WorkerFFTSize, 2))
        self.coalesced = 128 / 64 if data_type == "float" else 128 / 128
        self.state_vec = th.zeros(self.WorkerFFTSize, int(log(self.WorkerFFTSize, 2)))
        self.data = th.rand(self.WorkerFFTSize, dtype=th.cfloat)
        print(self.data)
        self.output = th.zeros(self.WorkerFFTSize, dtype=th.cfloat)
        for i in range(self.WorkerFFTSize):
            for j in range(int(log(self.WorkerFFTSize, 2))):
                self.state_vec[i, j] = (i // (2 ** j)) % 2

        N_tmp = N
        plan = []
        self.tensor_shape = []
        while N_tmp > 1:
            plan.append(self.WorkerFFTSize if N_tmp >= self.WorkerFFTSize else N_tmp)
            self.tensor_shape.append(int(plan[-1]))
            N_tmp /= self.WorkerFFTSize
            
        
        self.tensor_shape.reverse()
        print(plan)
        print(self.state_vec)
        print(self.tensor_shape)
        # assert 0

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

    def global2reg(self, if_shared):
        
        self.fft = ''''''
        if not if_shared:
            self.fft += f'''
            {self.data_type}* gPtr = {self.gPtr};
            gPtr += tid / {self.bdim / self.BS} * {self.strdBS}
            '''
            for i in range(self.WorkerFFTSize):
                self.fft += f'''
                gPtr += {self.bdim / self.BS * self.strdTensor};
                {self.rPtr}[{i}] = *gPtr;
                ''' 
        else:
            if strdBS == 1:
                self.fft += f'''
            {self.data_type}* gPtr = {self.gPtr} + tid % {self.BS // self.coalesced};
            gPtr +=  tid / {self.BS // self.coalesced} * {self.strdTensor};
            {self.data_type}* shPtr = {self.shPtr} + tid % {self.BS // self.coalesced} * {self.N};
            shPtr +=  tid / {self.BS // self.coalesced};
            int offset = 0;
                '''
                for i in range(self.WorkerFFTSize / self.coalesced):
                    self.fft += f'''
            offset = {i} * {(self.bdim // (self.BS // self.coalesced)) * self.strdTensor};
            (({self.data_type}{self.coalesced * 2}){self.rPtr})[{i}] = *({self.data_type}{self.coalesced * 2})(gPtr + offset);
            ''' 
            
                for i in range(self.WorkerFFTSize / self.coalesced):
                    self.fft += f'''
            offset = {i} * {(self.bdim // (self.BS // self.coalesced))};
            *({self.data_type}{self.coalesced * 2})(shPtr + offset) = (({self.data_type}{self.coalesced * 2}){self.rPtr})[{i}];
            ''' 
            else:
                self.fft += f'''
            {self.data_type}* gPtr = {self.gPtr} + tid % {min(self.N // self.coalesced, self.bdim)};
            gPtr +=  tid / {min(self.N // self.coalesced, self.bdim)} * {self.strdBS};
            {self.data_type}* shPtr = {self.shPtr} + tid % {min(self.N // self.coalesced, self.bdim)};
            shPtr +=  tid / {min(self.N // self.coalesced, self.bdim)} * {self.N};
            int offset = 0;
                '''
                for i in range(self.WorkerFFTSize / self.coalesced):
                    self.fft += f'''
            offset = {i // (self.N // min(self.N // self.coalesced, self.bdim))} * {self.strdBS};
            offset += {i % (self.N // min(self.N // self.coalesced, self.bdim))} * {min(self.N // self.coalesced, self.bdim)};
            (({self.data_type}{self.coalesced * 2}){self.rPtr})[{i}] = *({self.data_type}{self.coalesced * 2})(gPtr + offset);
            ''' 
                for i in range(self.WorkerFFTSize / self.coalesced):
                    self.fft += f'''
            offset = {i // (self.N // min(self.N // self.coalesced, self.bdim))} * {self.N};
            offset += {i % (self.N // min(self.N // self.coalesced, self.bdim))} * {min(self.N // self.coalesced, self.bdim)};
            *({self.data_type}{self.coalesced * 2})(shPtr + offset) = (({self.data_type}{self.coalesced * 2}){self.rPtr})[{i}];
            ''' 
                self.fft += f'''
                syncthreads();
                {self.data_type}* shPtr = {self.shPtr} + tid % {self.N // self.WorkerFFTSize};
                shPtr +=  tid / {self.N // self.WorkerFFTSize} * {self.N};
                '''
                for i in range(self.WorkerFFTSize):
                    self.fft += f'''
            offset = i * {self.N // self.WorkerFFTSize};
            (({self.data_type}2){self.rPtr})[{i}] = *({self.data_type}2)(shPtr + offset);
            ''' 
            
    def fft_reg(self, cur_stage):
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
    N = 128
    fft = TurboFFT(N=N)
    # print(fft.data)
    # torch_fft_result = th.fft.fft(fft.data.reshape(fft.WorkerFFTSize // 4, 4), dim=1).flatten()
    fft.fft_reg(0)
    fft.data = fft.output.clone()
    fft.fft_reg(1)
    fft.data = fft.output.clone()
    # fft.fft_reg(1)
    # fft.fft_reg(2)
    # torch_fft_result = th.fft.fft(fft.data).flatten()
    torch_fft_result = th.fft.fft(fft.data.reshape(fft.WorkerFFTSize // 4, 4), dim=0).flatten()
    print(fft.tensor_shape)
    # print(fft.data)
    # print(fft.data.reshape(4, -1)[:, 0])
    # print(torch_fft_result)
    # print( th.fft.fft(fft.data.reshape(4, -1)[:, 0]))
    # assert 0
    # print(fft.data)
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