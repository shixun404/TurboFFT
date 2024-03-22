import torch as th
from math import *
import argparse
import numpy as np
from main_codegen import main_codegen
import sys
class TurboFFT:
    def __init__(self, global_tensor_shape=[256, 1], radix=2, WorkerFFTSizes = [8],
                        threadblock_bs=[1], threadblock_bs_dim=[0], shared_mem_size=[0], data_type='double2',
                        if_special=False, if_ft=0, if_err_injection=0,  err_smoothing=1000, err_inj=100, err_threshold=1e-3):
        self.fft_code = []
        self.data_type = data_type
        self.gPtr = "gPtr"
        self.rPtr = "rPtr"
        self.rPtr_2 = "rPtr_2"
        self.rPtr_3 = "rPtr_3"
        self.rPtr_4 = "rPtr_4"
        self.shPtr = "shPtr"
        self.if_err_injection = if_err_injection
        self.err_inj = err_inj
        self.err_smoothing = err_smoothing
        self.err_threshold = err_threshold
        self.shared_mem_size=shared_mem_size
        self.thread_bs = 1
        self.ft = 0
        self.WorkerFFTSizes = WorkerFFTSizes
        self.threadblock_bs = threadblock_bs
        self.threadblock_bs_dim = threadblock_bs_dim
        self.global_tensor_shape = global_tensor_shape
        self.radix = radix
        self.state_vec = th.zeros(64, 6)
        self.if_special = if_special
        for i in range(64):
            for j in range(6):
                self.state_vec[i, j] = int((i // (2 ** j))) % 2

        self.threadblock_tensor_shape = []
        for size, N_tmp in zip(WorkerFFTSizes, self.global_tensor_shape[:-1]):
            threadblock_tensor_shape = []
            while N_tmp > 1:
                threadblock_tensor_shape.append(size if N_tmp >= size else int(N_tmp))
                N_tmp /= size
            threadblock_tensor_shape.reverse()
            self.threadblock_tensor_shape.append(threadblock_tensor_shape)
    def init(self, dim=0):
        self.local_variable = {
            "j" : ("int", "0"),
            "k" : ("int", "-1"),
            "global_j" : ("int", "0"),
            "global_k" : ("int", "0"),
            "data_id" : ("int", "0"),
            "bs_id" : ("int", "0"),
            "shared_offset_bs" : ("int", "0"),
            "shared_offset_data" : ("int", "0"),
            "bx": ("int", "gridDim.x - blockIdx.x - 1") if dim == 0 and len(self.WorkerFFTSizes) == 3 else ("int", "blockIdx.x"),
            "tx": ("int", "threadIdx.x"),
            "offset": ("int", "0"),
            self.gPtr: (f"{self.data_type}*", "inputs"),
            self.shPtr: (f"{self.data_type}*", "shared"),
            f"{self.rPtr}[{self.WorkerFFTSizes[dim]}]": (self.data_type, None),
            f"{self.rPtr_2}[{self.WorkerFFTSizes[dim] }]": (self.data_type, None),
            f"{self.rPtr_3}[{self.WorkerFFTSizes[dim] }]": (self.data_type, None),
            f"{self.rPtr_4}[{self.WorkerFFTSizes[dim] }]": (self.data_type, None),
            "tmp": (self.data_type, None),
            "tmp_1": (self.data_type, None),
            "tmp_2": (self.data_type, None),
            "angle": (self.data_type, None),
            "delta_angle": (self.data_type, None),
        }


    def save_generated_code(self, ):
        if not self.if_special:
            N = th.prod(th.as_tensor(self.global_tensor_shape[:-1]))
            for i in range(len(self.global_tensor_shape) - 1):
                file_name = f"../generated/{self.data_type}/fft_radix_{self.radix}_logN_{int(log(N, 2))}_upload_{i}.cuh"
                with open(file_name, 'w') as f:
                    f.write(self.fft_code[i])
        else:
            N = th.prod(th.as_tensor(self.global_tensor_shape[:-2]))
            file_name = f"../generated/{self.data_type}/fft_radix_{self.radix}_logN_{int(log(N, 2))}_upload_{0}.cuh"
            with open(file_name, 'w') as f:
                f.write(self.fft_code[1])

    def codegen(self,):

        reg_tensor_stride = th.as_tensor([1, 2, 4, 8, 16, 32, 64], dtype=th.float)
        state_vec = self.state_vec.clone()
        for dim in range(len(self.global_tensor_shape) - 2, -1, -1):
            self.init(dim)
            
            threadblock_tensor_shape =  self.threadblock_tensor_shape[dim]
            threadblock_bs = self.threadblock_bs[dim]
            threadblock_bs_dim = self.threadblock_bs_dim[dim]
            WorkerFFTSize = self.WorkerFFTSizes[dim]
            print( self.shared_mem_size,self.global_tensor_shape,self.threadblock_bs, dim)
            smem_size = self.shared_mem_size[dim - 1]
            logWorkerFFTSize = int(log(WorkerFFTSize, 2))
            global_tensor_shape = self.global_tensor_shape
            blockorder = [i for i in range(len(global_tensor_shape))]
            self.dim_ = dim
            fft_code = self.head(len(self.global_tensor_shape) - 2 - dim, global_tensor_shape,threadblock_bs, WorkerFFTSize,dim,smem_size)
            fft_code += self.globalAccess(dim, global_tensor_shape, 
                                    threadblock_bs_dim, threadblock_bs, WorkerFFTSize,
                                     blockorder)
            for threadblock_dim in range(len(threadblock_tensor_shape)):
                self.state_vec = state_vec[:, :logWorkerFFTSize]
                if threadblock_dim != 0:
                    fft_code += self.shared2reg(threadblock_bs, threadblock_tensor_shape,
                                            WorkerFFTSize, dim=threadblock_dim)
                fft_code += self.fft_reg(threadblock_bs, threadblock_tensor_shape, 
                                    WorkerFFTSize, threadblock_dim, reg_tensor_stride[:logWorkerFFTSize])
                dict_output = self.reg_output_remap(WorkerFFTSize // threadblock_tensor_shape[-1], 
                                                    reg_tensor_stride[:logWorkerFFTSize], WorkerFFTSize)
                threadblock_tensor_shape = threadblock_tensor_shape[:threadblock_dim] + \
                                            [threadblock_tensor_shape[-1]] + \
                                            threadblock_tensor_shape[threadblock_dim:-1]
                
                # if True:
                if threadblock_dim != len(threadblock_tensor_shape) - 1:
                
                    fft_code += self.reg2shared(threadblock_bs, threadblock_tensor_shape, 
                                        WorkerFFTSize, threadblock_dim, dict_output)
            if dim == 0 and not if_special:
                blockorder = self.list_reverse(blockorder, 0, -1)
                global_tensor_shape = self.list_reverse(global_tensor_shape, 0, -1)
            fft_code += self.globalAccess(blockorder[dim], global_tensor_shape, 
                        blockorder[threadblock_bs_dim], threadblock_bs, WorkerFFTSize,
                            blockorder, if_output=True, dict_output=dict_output, if_twiddle=(dim!=0))
            if self.ft == 1 and self.if_err_injection == 1:
                fft_code += '''
                }
                 if(k != -1){
                //if(false){
                '''
                fft_code += f'''
                // if(threadIdx.x == 0)printf("bid=%d, upload=%d, bx=%d, tx=%d, k = %d\\n", bid, {len(self.fft_code)}, blockIdx.x, threadIdx.x, k);
                bid = k;
                '''
                threadblock_tensor_shape =  self.threadblock_tensor_shape[dim]
                threadblock_bs = self.threadblock_bs[dim]
                threadblock_bs_dim = self.threadblock_bs_dim[dim]
                WorkerFFTSize = self.WorkerFFTSizes[dim]
                print( self.shared_mem_size,self.global_tensor_shape,self.threadblock_bs, dim)
                smem_size = self.shared_mem_size[dim - 1]
                logWorkerFFTSize = int(log(WorkerFFTSize, 2))
                global_tensor_shape = self.global_tensor_shape
                blockorder = [i for i in range(len(global_tensor_shape))]
                fft_code += self.globalAccess(dim, global_tensor_shape, 
                                    threadblock_bs_dim, threadblock_bs, WorkerFFTSize,
                                     blockorder, if_correction=self.if_err_injection)
                for threadblock_dim in range(len(threadblock_tensor_shape)):
                    self.state_vec = state_vec[:, :logWorkerFFTSize]
                    if threadblock_dim != 0:
                        fft_code += self.shared2reg(threadblock_bs, threadblock_tensor_shape,
                                                WorkerFFTSize, dim=threadblock_dim)
                    fft_code += self.fft_reg(threadblock_bs, threadblock_tensor_shape, 
                                        WorkerFFTSize, threadblock_dim, reg_tensor_stride[:logWorkerFFTSize])
                    dict_output = self.reg_output_remap(WorkerFFTSize // threadblock_tensor_shape[-1], 
                                                        reg_tensor_stride[:logWorkerFFTSize], WorkerFFTSize)
                    threadblock_tensor_shape = threadblock_tensor_shape[:threadblock_dim] + \
                                                [threadblock_tensor_shape[-1]] + \
                                                threadblock_tensor_shape[threadblock_dim:-1]
                    if threadblock_dim != len(threadblock_tensor_shape) - 1:
                    
                        fft_code += self.reg2shared(threadblock_bs, threadblock_tensor_shape, 
                                            WorkerFFTSize, threadblock_dim, dict_output)
                if dim == 0 and not if_special:
                    blockorder = self.list_reverse(blockorder, 0, -1)
                    global_tensor_shape = self.list_reverse(global_tensor_shape, 0, -1)
                fft_code += self.globalAccess(blockorder[dim], global_tensor_shape, 
                            blockorder[threadblock_bs_dim], threadblock_bs, WorkerFFTSize,
                                blockorder, if_output=True, dict_output=dict_output, if_twiddle=(dim!=0), if_correction=self.if_err_injection)
                fft_code += '''}}'''
            else:
                fft_code += self.epilogue()
            self.fft_code.append(fft_code)

    def head(self,  dim, global_tensor_shape, threadblock_bs, WorkerFFTSize, dim_, smem_size):
        N = th.prod(th.as_tensor(self.global_tensor_shape[:-1]))
        global_tensor_shape = th.as_tensor(global_tensor_shape)
        num_thread = (global_tensor_shape[dim_] / WorkerFFTSize * threadblock_bs)
        Ni = self.global_tensor_shape[len(self.global_tensor_shape) - 2 -dim]
        threadblock_bs = self.threadblock_bs[len(self.global_tensor_shape) - 2 -dim]
        if self.if_special:
            N = th.prod(th.as_tensor(self.global_tensor_shape[:-2]))
            dim = 0
        head = f'''extern __shared__ {self.data_type} shared[];
__global__ void fft_radix_{self.radix}_logN_{int(log(N, self.radix))}_dim_{dim}''' \
        + f'''({self.data_type}* inputs, {self.data_type}* outputs, {self.data_type}* twiddle, {self.data_type}* checksum_DFT, int BS)''' + ''' {
    int bid_cnt = 0;
    '''
        head += f'''
    int threadblock_per_SM = {int(128 * 1024 / (smem_size * 16 if self.data_type == "double2" else smem_size * 8))};
    int tb_gap = threadblock_per_SM * 108;
    int delta_bid = ((blockIdx.x / tb_gap) ==  (gridDim.x / tb_gap)) ? (gridDim.x % tb_gap) : tb_gap;
    {self.data_type} r[3];
    r[0].x = 1.0;
    r[0].y = 0.0;
    r[1].x = -0.5;
    r[1].y = -0.8660253882408142;
    r[2].x = -0.5;
    r[2].y = 0.8660253882408142;
    '''
        for key in self.local_variable.keys():
            head += f'''{self.local_variable[key][0]} {key};
    '''
        for key in self.local_variable.keys():
            if self.local_variable[key][1] is not None:
                head += f'''{key} = {self.local_variable[key][1]};
    '''
        if self.ft == 1:
            for i in range(max(1, global_tensor_shape[dim_] / (global_tensor_shape[dim_] / WorkerFFTSize * threadblock_bs))):
                head += f'''
    {self.rPtr_2}[{i}] = *(checksum_DFT + {global_tensor_shape[dim_]} - 2 + tx + {i * num_thread});
    {self.shPtr}[tx + {i * num_thread}] = {self.rPtr_2}[{i}];
    '''
            head += f'''
    __syncthreads();
    tmp_1.x = 0;
    tmp_1.y = 0;
    '''
            for i in range(WorkerFFTSize):
                head += f'''
    {self.rPtr_2}[{i}] = *({self.shPtr} +  tx / {threadblock_bs} + {i * (global_tensor_shape[dim_] // WorkerFFTSize)});
    {self.rPtr_3}[{i}].x = 0; {self.rPtr_3}[{i}].y = 0;
    {self.rPtr_4}[{i}].x = 0; {self.rPtr_4}[{i}].y = 0;
    '''
            head += f'''
    __syncthreads();
    '''
        head += f'''
    int bid = 0;
    for(bid = (blockIdx.x / tb_gap) * tb_gap * {self.thread_bs} + blockIdx.x % tb_gap;
                bid_cnt < {self.thread_bs} && bid < ({N} * BS + {Ni * threadblock_bs} - 1) / {Ni * threadblock_bs}; bid += delta_bid)
    '''
        head += '''{
    bid_cnt += 1;
    '''
        return head
    
    def epilogue(self, ):
        epilogue = '''
    }
    
}
'''
        return epilogue

    def globalAccess(self, dim, global_tensor_shape, threadblock_bs_dim, threadblock_bs, 
                    WorkerFFTSize, blockorder, if_output=False, dict_output=None, if_twiddle=False, if_to_shared=False, if_correction=False):
        global_tensor_shape = th.as_tensor(global_tensor_shape)
        threadblock_tensor_shape = th.ones_like(global_tensor_shape)
        threadblock_tensor_shape[dim] = global_tensor_shape[dim]
        threadblock_tensor_shape[threadblock_bs_dim] = threadblock_bs

        T = int(global_tensor_shape[dim] / WorkerFFTSize)
        num_thread = (global_tensor_shape[dim] / WorkerFFTSize * threadblock_bs)
        globalAccess_code = f'''        
    bx = bid;
    tx = threadIdx.x;
    ''' 
        access_stride_bs_dim = 1
        access_stride_data_dim =  1
        access_stride_bs = 1
        access_stride_data =  1
        if if_output is False:
            globalAccess_code += f'''
            tmp_1.x = 0;
            tmp_1.y = 0;
            tmp_2.x = 0;
            tmp_2.y = 0;
            {self.gPtr} = {self.local_variable[self.gPtr][1]};
    '''
        else:
            globalAccess_code += f'''{self.gPtr} = outputs;
    '''
        if if_twiddle:
            globalAccess_code += '''global_j = 0;
    global_k = 0;
    '''

        access_stride = 1
        
        for i in blockorder[:-1]:
            stride = max(1, th.prod(global_tensor_shape[:i]))
            if i < dim and if_twiddle:
                globalAccess_code += f'''
    global_j += (bx % {global_tensor_shape[i] // threadblock_tensor_shape[i]}) * {threadblock_tensor_shape[i]} * {stride};
    '''
                if i == threadblock_bs_dim:
                    globalAccess_code += f'''
    global_j += (tx % {threadblock_bs}) * {stride};
    '''    
            if i == dim:
                if not if_to_shared:
                    
                    globalAccess_code += f'''
    {self.gPtr} += tx / {threadblock_bs} * {stride};
    '''
                    access_stride = global_tensor_shape[i] // WorkerFFTSize * stride
                    
                else:
                    globalAccess_code += f'''
    {self.gPtr} += tx % {global_tensor_shape[i]} * {stride};
    shared_offset_data = tx % {global_tensor_shape[i]}; 
    '''         
                    access_stride_data_dim = max(1, global_tensor_shape[dim] // (global_tensor_shape[dim] // WorkerFFTSize * threadblock_bs))
                    access_stride_data = (global_tensor_shape[dim] // WorkerFFTSize * threadblock_bs) * stride
            
            globalAccess_code += f'''
    {self.gPtr} += (bx % {global_tensor_shape[i] // threadblock_tensor_shape[i]}) * {threadblock_tensor_shape[i]} * {stride};
    bx = bx / {global_tensor_shape[i] // threadblock_tensor_shape[i]};
    '''
            if i == threadblock_bs_dim:
                if not if_to_shared:
                    globalAccess_code += f'''
    {self.gPtr} += tx % {threadblock_bs} * {stride};
    '''
                else:
                    globalAccess_code += f'''
    shared_offset_bs = tx / {global_tensor_shape[dim]}; 
    {self.gPtr} += tx / {global_tensor_shape[dim]} * {stride};
    '''         
                    access_stride_bs_dim = threadblock_bs if (global_tensor_shape[dim] // WorkerFFTSize * threadblock_bs) >= 1 else global_tensor_shape[dim] // WorkerFFTSize 
                    access_stride_bs = stride if (global_tensor_shape[dim] // WorkerFFTSize * threadblock_bs) >= 1 else ((WorkerFFTSize  * threadblock_bs) / global_tensor_shape[dim] * stride)

        globalAccess_code += f'''
    {self.gPtr} += (bx % BS * {th.prod(global_tensor_shape[:-1])});
    '''
        if if_twiddle:
            globalAccess_code += f'''
    global_k += tx / {threadblock_bs};
    '''
        
        if not if_to_shared:
            if self.ft == 1 and if_output and not if_correction:
                globalAccess_code += f'''
        tmp_2 = tmp_1;    
        '''
            for i in range(WorkerFFTSize):
                if not if_output:
                    if not if_correction:
                        globalAccess_code += f'''
        {self.rPtr}[{i}] = *({self.gPtr} + {i * access_stride});
        {self.rPtr_3}[{i}].x += {self.rPtr}[{i}].x;
        {self.rPtr_3}[{i}].y += {self.rPtr}[{i}].y;
        '''
                        if self.ft == 1 and not if_correction:
                            globalAccess_code += f'''
        // tmp = checksum_DFT[tx / {threadblock_bs} + {i * (global_tensor_shape[self.dim_] // WorkerFFTSize)}];
        // turboFFT_ZMUL_ACC(tmp_1, {self.rPtr}[{i}], tmp);
        turboFFT_ZMUL_ACC(tmp_1, {self.rPtr}[{i}], {self.rPtr_2}[{i}])
        '''               
                    else:
                        globalAccess_code += f'''
        {self.rPtr}[{i}] = {self.rPtr_3}[{i}];
        // printf("bx=%d tx=%d, ptr_3_x=%f ptr_3_y=%f, ptr_4_x=%f, ptr_4_y=%f\\n", blockIdx.x, 
        // threadIdx.x, {self.rPtr_3}[{i}].x, {self.rPtr_3}[{i}].y, {self.rPtr_4}[{i}].x, {self.rPtr_4}[{i}].y);
        '''
                else:
                    if self.ft == 1 and not if_correction:
                        globalAccess_code += f'''
        turboFFT_ZMUL_NACC(tmp_1,  {self.rPtr}[{dict_output[i]}], r[({i * (global_tensor_shape[dim] // WorkerFFTSize)} + tx / {threadblock_bs}) % 3])
        '''
                    if if_twiddle:
                        N = th.prod(global_tensor_shape[:(dim + 1)])
                        if i == 0:
                            globalAccess_code += f'''
        delta_angle = twiddle[{N - 1} + global_j * ({global_tensor_shape[dim] // WorkerFFTSize})];
        angle = twiddle[{N - 1} + global_j * global_k];
        '''                    
                        else:
                            globalAccess_code += f'''
        tmp = angle;
        turboFFT_ZMUL(angle, tmp, delta_angle);
        '''                              
                        globalAccess_code += f'''
            tmp = {self.rPtr}[{dict_output[i]}];
            turboFFT_ZMUL({self.rPtr}[{dict_output[i]}], tmp, angle);
            '''
                    if if_correction:
                        globalAccess_code += f'''
            // if(threadIdx.x == 0 && blockIdx.x == 0)printf("i=%d, rx=%f ry=%f r4x=%f r4y=%f\\n", {i}, {self.rPtr}[{dict_output[i]}].x
            //  ,{self.rPtr}[{dict_output[i]}].y, {self.rPtr_4}[{dict_output[i]}].x, {self.rPtr_4}[{dict_output[i]}].y );
            turboFFT_ZSUB({self.rPtr}[{dict_output[i]}], {self.rPtr}[{dict_output[i]}], {self.rPtr_4}[{i}]);
            // turboFFT_ZSUB({self.rPtr}[{dict_output[i]}], {self.rPtr}[{dict_output[i]}], {self.rPtr_2}[{i}]);
            '''
            if self.ft == 1 and self.if_err_injection and not if_output and len(self.fft_code) == len(self.shared_mem_size) - 1 and not if_correction:
                globalAccess_code += f'''
            // if(threadIdx.x == 0 && blockIdx.x == 0 && bid_cnt == 1) printf("asdadas\\n");
            // if(threadIdx.x == 0 && blockIdx.x == 0 && bid_cnt == 1) {self.rPtr}[0].x += 100;
            {self.rPtr}[0].x += ((threadIdx.x + blockIdx.x) == 0 && bid_cnt == 1) ? {self.err_inj}: 0;
            '''    
                
        else:
            for j in range(access_stride_bs_dim):
                for i in range(access_stride_data_dim):
                    globalAccess_code += f'''
        bs_id = shared_offset_bs + {int(j * max(1, (T * threadblock_bs) // global_tensor_shape[dim]))};
        data_id = shared_offset_data + {int(i * ((T * threadblock_bs)))};
        *({self.shPtr} + (bs_id + data_id) % {global_tensor_shape[dim]} + bs_id * {global_tensor_shape[dim]}) = *({self.gPtr} + {i * access_stride_data + j * access_stride_bs});
        '''
            globalAccess_code += f'''
        __syncthreads();
        '''
            for i in range(WorkerFFTSize):
                globalAccess_code += f'''
        bs_id = threadIdx.x % {threadblock_bs};
        data_id = threadIdx.x / {threadblock_bs} + {i * T};
        {self.rPtr}[{i}] = *({self.shPtr} + (bs_id + data_id) % {global_tensor_shape[dim]} + bs_id * {global_tensor_shape[dim]});
        '''

        if self.ft == 1 and if_output and not if_correction:
            globalAccess_code += f'''
        tmp_1.y = tmp_1.x + tmp_1.y;
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 16, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 8, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 4, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 2, 32);
        tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, 1, 32);

         tmp_1.x = tmp_2.x + tmp_2.y;
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 16, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 8, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 4, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 2, 32);
         tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, 1, 32);
        __syncthreads();
        // if(tx % 32 == 0) {self.shPtr}[tx / 32] = tmp_1;
        {self.shPtr}[tx / 32] = tmp_1;
        __syncthreads();
        // if(tx < 32)
        '''
            # globalAccess_code += '''{'''
            globalAccess_code += f'''
            tmp_1 = {self.shPtr}[tx % {num_thread // 32}];
        '''
            i = num_thread // 32
            while( i > 1):
                i //= 2
                globalAccess_code += f'''
                tmp_1.y += __shfl_xor_sync(0xffffffff, tmp_1.y, {i}, 32);
                tmp_1.x += __shfl_xor_sync(0xffffffff, tmp_1.x, {i}, 32);
        '''
            globalAccess_code  += f'''
            // if(tx == 0 && abs(tmp_1.y) / ({self.err_smoothing} + abs(tmp_1.x)) > 1e-3)printf("{len(self.fft_code)}, bid=%d bx=%d, by=%d, tx=%d: %f, %f, %f\\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            // if(abs(tmp_1.y / tmp_1.x) > {self.err_threshold})printf("{len(self.fft_code)}, bid=%d bx=%d, by=%d, tx=%d: %f, %f, %f\\n", bid, blockIdx.x, blockIdx.y, threadIdx.x, tmp_1.x, tmp_1.y, tmp_1.y / tmp_1.x);
            k = abs(tmp_1.y) / ({self.err_smoothing} + abs(tmp_1.x)) > {self.err_threshold} ? bid : k;
            // k = abs(tmp_1.y) > 10 ? bid : k;
            // if(tx == 0) *({self.gPtr}) = tmp_1;
            // if(tx == 0 && abs(tmp_1.y / tmp_1.x) > 1e-3)
            '''
            globalAccess_code += '''
            // }
            // }            
            '''
        for i in range(WorkerFFTSize):
            if if_output:
                if not if_correction:
                    globalAccess_code += f'''
            *({self.gPtr} + {i * access_stride}) = {self.rPtr}[{dict_output[i]}];
            {self.rPtr_4}[{i}].x += {self.rPtr}[{dict_output[i]}].x;
            {self.rPtr_4}[{i}].y += {self.rPtr}[{dict_output[i]}].y;
            // {self.rPtr_2}[{i}].x += {self.rPtr}[{dict_output[i]}].x;
            // {self.rPtr_2}[{i}].y += {self.rPtr}[{dict_output[i]}].y;
            '''     
                else:
                    globalAccess_code += f'''
            {self.rPtr_3}[{i}] = *({self.gPtr} + {i * access_stride});
            turboFFT_ZADD({self.rPtr_3}[{i}], {self.rPtr_3}[{i}], {self.rPtr}[{dict_output[i]}] );
            // printf("bx=%d tx=%d d_x=%f d_y=%f\\n", blockIdx.x, threadIdx.x, {self.rPtr}[{dict_output[i]}].x, {self.rPtr}[{dict_output[i]}].y);
            *({self.gPtr} + {i * access_stride}) = {self.rPtr_3}[{i}];
            '''                    
        return globalAccess_code

    def list_reverse(self, list_, st, end):
        if isinstance(list_, list):
            target = list_[st:end]
            target.reverse()
            target = list_[:st] + target + list_[end:]
        else:
            target = th.cat((list_[:st], list_[st:end].flip(0), list_[end:]), dim=0)
        return target
            
    def shared2reg(self, threadblock_bs, threadblock_tensor_shape, WorkerFFTSize, dim=None):
        shared2reg_code  = ''''''
        access_stride = int(threadblock_bs * th.prod(th.as_tensor(threadblock_tensor_shape))
                         / WorkerFFTSize)
        dim_0 = threadblock_tensor_shape[0]
        dim_1 = threadblock_tensor_shape[1]
        
        if len(threadblock_tensor_shape) == 2 and len(self.global_tensor_shape) == 2 :
            shared2reg_code += f'''
    offset = 0;
    '''
            shared2reg_code += f'''
    __syncthreads();
    '''

            
            for j in range(WorkerFFTSize):
                i = j % WorkerFFTSize
                shared2reg_code += f'''
        {self.rPtr}[{i}] = {self.shPtr}[{access_stride * i} + (tx / {dim_1}) * {dim_1} + (tx + {i}) % {dim_1}];
        '''
            return shared2reg_code

        if dim == 1 and len(self.global_tensor_shape) == 2 :
            shared2reg_code += f'''
    offset = 0;
    offset += (tx / {dim_0}) * {dim_0} + 
              ((tx % {dim_0}) / {dim_1}) * {dim_1} + (tx % {dim_1} + tx / {dim_0}) % {dim_1};
    '''
        else:
            shared2reg_code += f'''
    offset = 0;
    offset += tx;
    '''
        
        shared2reg_code += '''
    __syncthreads();
    '''
        for j in range(WorkerFFTSize):
            i = j % WorkerFFTSize
            shared2reg_code += f'''
    {self.rPtr}[{i}] = {self.shPtr}[offset + {access_stride * i}];
    '''
        return shared2reg_code

    def reg_output_remap(self, bs, reg_tensor_stride, WorkerFFTSize):
        logbs = int(log(bs, 2))
        # Keep the leading stride of batch size, flip the following tensor stride
        # reg_tensor_stride_reverse = th.cat((reg_tensor_stride[:logbs], reg_tensor_stride[logbs:].flip(0)), dim=0)
        
        reg_tensor_stride_reverse = self.list_reverse(reg_tensor_stride, logbs, len(reg_tensor_stride))
        
        # save output_id as dict so that we can visit it in an increment order
        dict_output = {}
        for i in range(WorkerFFTSize):
            output_id = int(th.dot(self.state_vec[i], reg_tensor_stride_reverse))
            dict_output[output_id] = i
        return dict_output

    def reg2shared(self, threadblock_bs, threadblock_tensor_shape, WorkerFFTSize, dim, dict_output):
        # Todo: swizzling not worked for dim1 < dim0
        reg2shared_code = '''
    j = 0;
    offset  = 0;
    '''
        if dim == len(threadblock_tensor_shape) - 1:
            access_stride = int(threadblock_bs * th.prod(th.as_tensor(threadblock_tensor_shape))
                         / WorkerFFTSize)
            reg2shared_code += f'''
    __syncthreads();
    '''
            for i in range(WorkerFFTSize):
                reg2shared_code += f'''
    {self.shPtr}[threadIdx.x + {access_stride * i}] = {self.rPtr}[{dict_output[i]}];
    '''
            reg2shared_code += f'''
    __syncthreads();
    '''
            return reg2shared_code

        # threadId to tensor coordinates
        tmp = 1
        stride = 1
        access_stride = 1
        bs_tensor_shape = [threadblock_bs] + threadblock_tensor_shape
        for i in range(len(bs_tensor_shape)):
            stride *= bs_tensor_shape[i]
            if i == dim + 1:
                access_stride = int(stride / bs_tensor_shape[i])
                reg2shared_code += f'''
    j = tx / {tmp};
    '''
                continue
            reg2shared_code += f'''
    offset += ((tx / {tmp}) % {bs_tensor_shape[i]}) * {int(stride / bs_tensor_shape[i])};
    '''
            tmp *= bs_tensor_shape[i]
            

        if dim == len(threadblock_tensor_shape) - 1:
            access_stride = int(threadblock_bs * th.prod(th.as_tensor(threadblock_tensor_shape)) / self.WorkerFFTSize)
        reg2shared_code += '''
    __syncthreads();
    '''
        N = th.prod(th.as_tensor(threadblock_tensor_shape[dim:]))
        for output_id in range(WorkerFFTSize): 
            # print(output_id, dict_output[output_id])
            if dim != len(threadblock_tensor_shape) - 1:
                if output_id == 0:
                    reg2shared_code += f'''
    delta_angle = twiddle[{N - 1} + j];
    angle.x = 1;
    angle.y = 0;
    '''       
                else:
                    reg2shared_code += f'''
    tmp = angle;
    turboFFT_ZMUL(angle, tmp, delta_angle);
    tmp = {self.rPtr}[{dict_output[output_id]}];
    turboFFT_ZMUL({self.rPtr}[{dict_output[output_id]}], tmp, angle);
    '''             

            if dim == 0 and len(self.global_tensor_shape) == 2 :
                reg2shared_code += f'''
    // {self.rPtr_2}[{output_id}] = {self.rPtr}[{dict_output[output_id]}];
    // {self.rPtr_2}[{output_id}] = {self.rPtr}[{dict_output[output_id]}];
    {self.shPtr}[offset + {access_stride} * ({output_id} + threadIdx.x % {threadblock_tensor_shape[1]}) % {threadblock_tensor_shape[1]} + ({output_id} / {threadblock_tensor_shape[1]}) * {threadblock_tensor_shape[1]}] = {self.rPtr}[{dict_output[output_id]}];
    '''             
            else:
                reg2shared_code += f'''
    {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr}[{dict_output[output_id]}];
    '''
        return reg2shared_code
    
    def fft_reg(self, threadblock_bs, threadblock_tensor_shape, WorkerFFTSize, dim, reg_tensor_stride):
        fft_reg_code = ''''''
        
        bs = WorkerFFTSize // threadblock_tensor_shape[-1]
        logbs = int(log(bs, 2))

        logWorkerFFTSize = int(log(WorkerFFTSize, 2))
        st = logWorkerFFTSize - 1
        
        for i in range(st, logbs - 1, -1):
            for j in range(WorkerFFTSize):
                if self.state_vec[j, i] == 1:
                    continue
                id_j1 = int(th.dot(self.state_vec[j], reg_tensor_stride))
                id_j2 = int(id_j1 + reg_tensor_stride[i])
                id_k = int(th.dot(self.state_vec[j, logbs:i], reg_tensor_stride[:i - logbs]))                
                tmp_angle = (-2 * id_k * 1 / (2 ** (i + 1 - logbs))) * pi
                rel_bounds = 1e-8
                abs_bounds = 1e-5
                fft_reg_code += f'''
    tmp = {self.rPtr}[{id_j1}];
    turboFFT_ZADD({self.rPtr}[{id_j1}], tmp, {self.rPtr}[{id_j2}]);
    turboFFT_ZSUB({self.rPtr}[{id_j2}], tmp, {self.rPtr}[{id_j2}]);
    tmp = {self.rPtr}[{id_j2}];
    '''
                if np.allclose(0, cos(tmp_angle), rel_bounds, abs_bounds):
                    if np.allclose(1, sin(tmp_angle), rel_bounds, abs_bounds):
                        fft_reg_code += f'''
    {self.rPtr}[{id_j2}].y = tmp.x;
    {self.rPtr}[{id_j2}].x = -tmp.y;
    '''
                    else:
                        fft_reg_code += f'''
    {self.rPtr}[{id_j2}].y = -tmp.x;
    {self.rPtr}[{id_j2}].x = tmp.y;
    '''
                elif np.allclose(0, sin(tmp_angle), rel_bounds, abs_bounds):
                    if np.allclose(1, cos(tmp_angle), rel_bounds, abs_bounds):
                        pass
                    else:
                        fft_reg_code += f'''
    {self.rPtr}[{id_j2}].x = -tmp.x;    
    {self.rPtr}[{id_j2}].y = -tmp.y;
    '''
                else:
                    if self.data_type == 'double2':
                        fft_reg_code += f'''
        angle.x = {cos(tmp_angle)};
        angle.y = {sin(tmp_angle)};
        turboFFT_ZMUL({self.rPtr}[{id_j2}], tmp, angle);
        '''
                    else:
                        fft_reg_code += f'''
        angle.x = {cos(tmp_angle)}f;
        angle.y = {sin(tmp_angle)}f;
        turboFFT_ZMUL({self.rPtr}[{id_j2}], tmp, angle);
        '''             
        
        return fft_reg_code

if __name__ == '__main__':
    params = []
    # datatype = 'float2'
    # datatype = 'double2'
    parser = argparse.ArgumentParser(description="turboFFT.")
    parser.add_argument('--if_ft', type=int, default=0, 
                        help='Flag to indicate feature transformation (0 for False, 1 for True)')
    parser.add_argument('--if_err_injection', type=int, default=0, 
                        help='Flag to indicate error injection (0 for False, 1 for True)')
    parser.add_argument('--err_smoothing', type=int, default=1000, 
                        help='Error smoothing parameter')
    parser.add_argument('--err_inj', type=int, default=100, 
                        help='Error injection parameter')
    parser.add_argument('--err_threshold', type=float, default=1e-3, 
                        help='Error threshold')
    parser.add_argument('--datatype', type=str, default="double2", 
                        help='Data type with a default value of "double2"')


    # Parse the arguments
    args = parser.parse_args()

    # Access arguments
    if_ft = args.if_ft
    if_err_injection = args.if_err_injection
    err_smoothing = args.err_smoothing
    err_inj = args.err_inj
    err_threshold = args.err_threshold
    datatype = args.datatype
    
    with open(f"../../param/param_A100_{datatype}.csv", 'r') as file:
        for line in file:
            # Splitting each line by comma
            split_elements = line.strip().split(',')
            row = [int(element) for element in split_elements]
            params.append(row)
    
    for row in params:
        global_tensor_shape = [2 ** i for i in row[2:(2 + row[1])]]
        threadblock_bs = row[5:(5 + row[1])]
        WorkerFFTSizes = row[8:(8 + row[1])]
        threadblock_bs.reverse()
        global_tensor_shape.reverse()
        WorkerFFTSizes.reverse()
        shared_mem_size = [global_tensor_shape[i] * threadblock_bs[i] for i in range(row[1])]
        # print(shared_mem_size, threadblock_bs, row[5:(5 + row[1])])
        
        st = row[1]
        if_special = False
        if row[1] == 1 and threadblock_bs[-1] != 1:
            st = 2
            if_special = True
            global_tensor_shape.append(threadblock_bs[-1])
            threadblock_bs.append(1)
            WorkerFFTSizes.append(2)
            
        global_tensor_shape.append(1)
        threadblock_bs_dim = [[1], [1, 0], [2, 0, 0]]
        
        fft = TurboFFT(global_tensor_shape=global_tensor_shape, WorkerFFTSizes=WorkerFFTSizes,
                    threadblock_bs=threadblock_bs, threadblock_bs_dim=threadblock_bs_dim[st - 1], shared_mem_size=shared_mem_size, data_type=datatype, if_special=if_special,
                    if_ft=if_ft, if_err_injection=if_err_injection, err_inj=err_inj, 
                    err_smoothing=err_smoothing, err_threshold=err_threshold)
        fft.codegen()
        fft.save_generated_code()
        main_code = main_codegen(datatype, fft.thread_bs)
        file_name = "../../../main.cu"
        with open(file_name, 'w') as f:
            f.write(main_code)
