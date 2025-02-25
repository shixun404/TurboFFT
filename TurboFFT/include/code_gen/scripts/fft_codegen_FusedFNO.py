import torch as th
from math import *
import argparse
import numpy as np
from main_codegen import main_codegen
import sys
M_PI = 3.141592653589793
class TurboFFT:
    def __init__(self, global_tensor_shape=[256, 1], radix=2, WorkerFFTSizes = [8],
                        threadblock_bs=[1], threadblock_bs_dim=[0], shared_mem_size=[0], 
                        data_type='double2', if_special=False, if_thread_ft=0, if_ft=0, 
                        if_err_injection=0,  err_smoothing=1000, err_inj=100,
                        err_threshold=1e-3, if_write=True):
        self.fft_code = []
        self.data_type = data_type
        self.gPtr = "gPtr"
        self.rPtr = "rPtr"
        self.rPtr_2 = "rPtr_2"
        self.rPtr_3 = "rPtr_3"
        self.rPtr_4 = "rPtr_4"
        self.shPtr = "shPtr"
        self.if_thread_ft = if_thread_ft
        self.if_err_injection = if_err_injection
        self.err_inj = err_inj
        self.err_smoothing = err_smoothing
        self.err_threshold = err_threshold
        self.shared_mem_size=shared_mem_size
        self.ft = if_ft
        self.WorkerFFTSizes = WorkerFFTSizes
        self.threadblock_bs = threadblock_bs
        self.threadblock_bs_dim = threadblock_bs_dim
        self.global_tensor_shape = global_tensor_shape
        self.radix = radix
        self.state_vec = th.zeros(64, 6)
        self.if_special = if_special
        self.if_write = if_write
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
            self.gPtr: (f"{self.data_type}*", "gPtr_1"),
            self.shPtr: (f"{self.data_type}*", "(float2*) shared_mem"),
            f"{self.rPtr}[{self.WorkerFFTSizes[dim]}]": (self.data_type, None),
            f"{self.rPtr_3}[{self.WorkerFFTSizes[dim] }]": (self.data_type, None),
            "tmp": (self.data_type, None),
            "angle": (self.data_type, None),
            "delta_angle": (self.data_type, None),
        }


    def save_generated_code(self, ):
        if not self.if_special:
            N = th.prod(th.as_tensor(self.global_tensor_shape[:-1]))
            for i in range(3):
                file_name = f"../generated/{self.data_type}/fft_radix_{self.radix}_logN_{int(log(N, 2))}_upload_{i}.cuh"
                if i >= len(self.global_tensor_shape) - 1:
                    with open(file_name, 'w') as f:
                        f.write("\n")
                else:
                    # if self.ft == 0:
                    if ( self.ft == 0 and self.if_thread_ft == 0) or self.if_write is True:
                        with open(file_name, 'w') as f:
                            f.write(self.fft_code[i])
                    else:
                        with open(file_name, 'a') as f:
                            f.write(self.fft_code[i])
                
        else:
            N = th.prod(th.as_tensor(self.global_tensor_shape[:-2]))
            for i in range(3):
                if i != 0:
                    file_name = f"../generated/{self.data_type}/fft_radix_{self.radix}_logN_{int(log(N, 2))}_upload_{i}.cuh"
                    with open(file_name, 'w') as f:
                        f.write("\n")
                else:
                    file_name = f"../generated/{self.data_type}/fft_radix_{self.radix}_logN_{int(log(N, 2))}_upload_{0}.cuh"
                    # if self.ft == 0:
                    if ( self.ft == 0 and self.if_thread_ft == 0) or self.if_write is True:
                        with open(file_name, 'w') as f:
                            f.write(self.fft_code[1])
                    else:
                        with open(file_name, 'a') as f:
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
                                            WorkerFFTSize, dim=threadblock_dim, global_dim=dim)
                fft_code += self.fft_reg(threadblock_bs, threadblock_tensor_shape, 
                                    WorkerFFTSize, threadblock_dim, reg_tensor_stride[:logWorkerFFTSize])
                dict_output = self.reg_output_remap(WorkerFFTSize // threadblock_tensor_shape[-1], 
                                                    reg_tensor_stride[:logWorkerFFTSize], WorkerFFTSize)
                threadblock_tensor_shape = threadblock_tensor_shape[:threadblock_dim] + \
                                            [threadblock_tensor_shape[-1]] + \
                                            threadblock_tensor_shape[threadblock_dim:-1]
                if threadblock_dim != len(threadblock_tensor_shape) - 1:
                
                    fft_code += self.reg2shared(threadblock_bs, threadblock_tensor_shape, 
                                        WorkerFFTSize, threadblock_dim, dict_output, dim)
            if dim == 0 and not if_special:
                blockorder = self.list_reverse(blockorder, 0, -1)
                global_tensor_shape = self.list_reverse(global_tensor_shape, 0, -1)
            fft_code += self.globalAccess(blockorder[dim], global_tensor_shape, 
                        blockorder[threadblock_bs_dim], threadblock_bs, WorkerFFTSize,
                            blockorder, if_output=True, dict_output=dict_output, if_twiddle=(dim!=0))
            fft_code += self.epilogue()
            self.fft_code.append(fft_code)

    def head(self,  dim, global_tensor_shape, threadblock_bs, WorkerFFTSize, dim_, smem_size):
        N = th.prod(th.as_tensor(self.global_tensor_shape[:-1]))
        global_tensor_shape = th.as_tensor(global_tensor_shape)
        num_thread = (global_tensor_shape[dim_] // WorkerFFTSize * threadblock_bs)
        Ni = self.global_tensor_shape[len(self.global_tensor_shape) - 2 -dim]
        threadblock_bs = self.threadblock_bs[len(self.global_tensor_shape) - 2 -dim]
        if self.if_special:
            N = th.prod(th.as_tensor(self.global_tensor_shape[:-2]))
            dim = 0
        head = f'''
extern __shared__ float shared_mem[];
__global__ void fft_{int(log(N, self.radix))}''' \
        + f'''(float2* gPtr_1, float2* outputs, int threadblock_bs)''' + ''' {
    int bid_cnt = 0;
    '''
    #     head += f'''
    # {self.data_type}* shared = ({self.data_type}*) ext_shared;
    # '''
        for key in self.local_variable.keys():
            head += f'''{self.local_variable[key][0]} {key};
    '''
        for key in self.local_variable.keys():
            if self.local_variable[key][1] is not None:
                head += f'''{key} = {self.local_variable[key][1]};
    '''
    
        head += f'''
    int bid = 0;
    '''
        return head
    
    def epilogue(self, ):
        epilogue = '''
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
        num_thread = (global_tensor_shape[dim] // WorkerFFTSize * threadblock_bs)
        globalAccess_code = f'''        
    bx = blockIdx.x;
    tx = threadIdx.x;
    ''' 
        if if_output is False:
            globalAccess_code += f'''
        {self.gPtr} += threadIdx.x % {global_tensor_shape[dim] // WorkerFFTSize};
    
        {self.gPtr} += (threadIdx.x / {global_tensor_shape[dim] // WorkerFFTSize}) * stride;
        __syncthreads();    
'''
        else:
            globalAccess_code += f'''{self.gPtr} = {self.shPtr};
        {self.gPtr} += threadIdx.x % {global_tensor_shape[dim] // WorkerFFTSize};
        {self.gPtr} += (threadIdx.x / {global_tensor_shape[dim] // WorkerFFTSize}) * {global_tensor_shape[dim]};
        __syncthreads();
    '''
#         globalAccess_code += f'''
#         {self.gPtr} += threadIdx.x % {global_tensor_shape[dim] // WorkerFFTSize};
    
#         {self.gPtr} += (threadIdx.x / {global_tensor_shape[dim] // WorkerFFTSize}) * {global_tensor_shape[dim]};
#         __syncthreads();
# '''
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
                    access_stride = global_tensor_shape[i] // WorkerFFTSize * stride

        if if_twiddle:
            globalAccess_code += f'''
    global_k += tx / {threadblock_bs};
    '''
        
        if not if_to_shared:
            for i in range(WorkerFFTSize):
                if not if_output:
                    if not if_correction:
                        globalAccess_code += f'''
        {self.rPtr}[{i}] = *({self.gPtr} + {i * access_stride});
        '''
                else:
                    if if_twiddle:
                        N = th.prod(global_tensor_shape[:(dim + 1)])
                        if i == 0:
                            globalAccess_code += f'''
        delta_angle = twiddle[{N - 1} + global_j * ({global_tensor_shape[dim] // WorkerFFTSize})];
        angle = twiddle[{N - 1} + global_j * global_k];
        ''' if self.data_type == 'double2' else f'''
        delta_angle.x = __cosf(global_j *  {-2.0 * M_PI * (global_tensor_shape[dim] // WorkerFFTSize) / float(N)}f);
        delta_angle.y = __sinf(global_j *  {-2.0 * M_PI * (global_tensor_shape[dim] // WorkerFFTSize) / float(N)}f);
        angle.x = __cosf( global_j * global_k * {-2.0 * M_PI / float(N)}f);
        angle.y = __sinf( global_j * global_k * {-2.0 * M_PI / float(N)}f);
        '''

                        else:
                            globalAccess_code += f'''
        tmp = angle;
        turboFFT_ZMUL{'_THREAD_FT' if self.if_thread_ft else ''}(angle, tmp, delta_angle);
        '''                              
                        globalAccess_code += f'''
            tmp = {self.rPtr}[{dict_output[i]}];
            turboFFT_ZMUL{'_THREAD_FT' if self.if_thread_ft else ''}({self.rPtr}[{dict_output[i]}], tmp, angle);
            '''

        for i in range(WorkerFFTSize):
            if if_output:
                globalAccess_code += f'''
            *({self.gPtr} + {i * access_stride}) = {self.rPtr}[{dict_output[i]}];
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
            
    def shared2reg(self, threadblock_bs, threadblock_tensor_shape, WorkerFFTSize, dim=None, global_dim=None):
        shared2reg_code  = ''''''
        access_stride = int(threadblock_bs * th.prod(th.as_tensor(threadblock_tensor_shape))
                         / WorkerFFTSize)
        dim_0 = threadblock_tensor_shape[0]
        dim_1 = threadblock_tensor_shape[1]
        
        # print("shared2reg", threadblock_tensor_shape, dim)
        # if dim == 1 and len(self.global_tensor_shape) == 2 :
        shared2reg_code += f'''
    offset = 0;
    offset += tx % {th.prod(th.as_tensor(threadblock_tensor_shape)) // WorkerFFTSize} + tx / {th.prod(th.as_tensor(threadblock_tensor_shape)) // WorkerFFTSize} * {th.prod(th.as_tensor(threadblock_tensor_shape))};
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

    def reg2shared(self, threadblock_bs, threadblock_tensor_shape, WorkerFFTSize, dim, dict_output, global_dim):
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
    j = (threadIdx.x % {self.global_tensor_shape[global_dim] // WorkerFFTSize}) / {tmp};
    '''
                continue
            reg2shared_code += f'''
    offset += ((threadIdx.x / {tmp}) % {bs_tensor_shape[i]}) * {int(stride / bs_tensor_shape[i])};
    '''
            tmp *= bs_tensor_shape[i]
        reg2shared_code += f'''
    offset += (threadIdx.x / {self.global_tensor_shape[global_dim] // WorkerFFTSize}) * {self.global_tensor_shape[global_dim]};
    '''    

        if dim == len(threadblock_tensor_shape) - 1:
            access_stride = int(threadblock_bs * th.prod(th.as_tensor(threadblock_tensor_shape)) / self.WorkerFFTSize)
        reg2shared_code += '''
    __syncthreads();
    '''
        N = th.prod(th.as_tensor(threadblock_tensor_shape[dim:]))
        print("reg2shared", self.global_tensor_shape, threadblock_tensor_shape, dim)
        
        
        for output_id in range(WorkerFFTSize): 
            # print(output_id, dict_output[output_id])
            if dim != len(threadblock_tensor_shape) - 1:
                if output_id == 0:
                    reg2shared_code += f'''
    delta_angle = twiddle[{N - 1} + j];
    ''' if self.data_type == 'double2' else f'''
    delta_angle.x = __cosf(j * {-2.0 * M_PI / N}f);
    delta_angle.y = __sinf(j * {-2.0 * M_PI / N}f);
    '''
                    reg2shared_code += f''' 
    angle.x = 1;
    angle.y = 0;
    '''       
                else:
                    reg2shared_code += f'''
    tmp = angle;
    turboFFT_ZMUL{'_THREAD_FT' if self.if_thread_ft else ''}(angle, tmp, delta_angle);
    tmp = {self.rPtr}[{dict_output[output_id]}];
    turboFFT_ZMUL{'_THREAD_FT' if self.if_thread_ft else ''}({self.rPtr}[{dict_output[output_id]}], tmp, angle);
    '''   
    #         if dim == 0 and len(self.global_tensor_shape) == 2 :
    #             reg2shared_code += f'''
    # // {self.shPtr}[offset + {access_stride} * (({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})] = {self.rPtr_3}[(({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})];
    # // {self.shPtr}[offset + {access_stride} * (({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})] = {self.rPtr}[{dict_output[output_id]}];
    # // {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr}[{dict_output[output_id]}];
    # //  {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr_3}[{output_id}];
    # '''
    #         elif dim == 1 and len(self.global_tensor_shape) == 2 :
    #             reg2shared_code += f'''
    # // {self.shPtr}[offset + {access_stride} * (({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})] = {self.rPtr_3}[(({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})];
    # // {self.shPtr}[offset + {access_stride} * (({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})] = {self.rPtr}[{dict_output[output_id]}];
    # // {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr}[{dict_output[output_id]}];
    # // {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr_3}[{output_id}];
    # '''
    #         else:
    #             reg2shared_code += f'''
    # {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr}[{dict_output[output_id]}];
    # '''
        for output_id in range(WorkerFFTSize): 
            reg2shared_code += f'''
            {self.rPtr_3}[{output_id}] = {self.rPtr}[{dict_output[output_id]}];
    '''
        for output_id in range(WorkerFFTSize):        
            if dim == 0 and len(self.global_tensor_shape) == 2 :
                reg2shared_code += f'''
    {self.shPtr}[offset + {access_stride} * (({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})] = {self.rPtr_3}[(({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})];
    // {self.shPtr}[offset + {access_stride} * (({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})] = {self.rPtr_3}[{output_id}];
     // {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr}[{dict_output[output_id]}];
    //  {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr_3}[{output_id}];
    '''
            elif dim == 1 and len(self.global_tensor_shape) == 2 :
                reg2shared_code += f'''
    {self.shPtr}[offset + {access_stride} * (({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})] = {self.rPtr_3}[(({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})];
    // {self.shPtr}[offset + {access_stride} * (({output_id} + (threadIdx.x / {(16 + WorkerFFTSize - 1) // WorkerFFTSize})) % {WorkerFFTSize})] = {self.rPtr_3}[{output_id}];
     // {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr}[{dict_output[output_id]}];
     // {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr_3}[{output_id}];
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
    turboFFT_ZADD{'_THREAD_FT' if self.if_thread_ft else ''}({self.rPtr}[{id_j1}], tmp, {self.rPtr}[{id_j2}]);
    turboFFT_ZSUB{'_THREAD_FT' if self.if_thread_ft else ''}({self.rPtr}[{id_j2}], tmp, {self.rPtr}[{id_j2}]);
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
        turboFFT_ZMUL{'_THREAD_FT' if self.if_thread_ft else ''}({self.rPtr}[{id_j2}], tmp, angle);
        '''
                    else:
                        fft_reg_code += f'''
        angle.x = {cos(tmp_angle)}f;
        angle.y = {sin(tmp_angle)}f;
        turboFFT_ZMUL{'_THREAD_FT' if self.if_thread_ft else ''}({self.rPtr}[{id_j2}], tmp, angle);
        '''             
        
        return fft_reg_code

if __name__ == '__main__':
    params = []
    parser = argparse.ArgumentParser(description="turboFFT.")
    parser.add_argument('--if_thread_ft', type=int, default=0, 
                        help='Flag to indicate thread level ft (0 for False, 1 for True)')
    parser.add_argument('--if_ft', type=int, default=0, 
                        help='Flag to indicate threadblock level ft (0 for False, 1 for True)')
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
    parser.add_argument('--gpu', type=str, default="A100", 
                        help='GPU spec a default value of "A100"')


    # Parse the arguments
    args = parser.parse_args()

    # Access arguments
    if_ft = args.if_ft
    if_thread_ft = args.if_thread_ft
    if_err_injection = args.if_err_injection
    err_smoothing = args.err_smoothing
    err_inj = args.err_inj
    err_threshold = args.err_threshold
    datatype = args.datatype
    gpu = args.gpu
    
    with open(f"../../param/{gpu}/param_{datatype}.csv", 'r') as file:
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
                    if_thread_ft=if_thread_ft, if_ft=if_ft, if_err_injection=if_err_injection, err_inj=err_inj, 
                    err_smoothing=err_smoothing, err_threshold=err_threshold, if_write=False)
        fft.codegen()
        fft.save_generated_code()
