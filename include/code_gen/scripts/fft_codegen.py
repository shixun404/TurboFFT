import torch as th
from math import *
import numpy as np
class TurboFFT:
    def __init__(self, global_tensor_shape=[256, 1], radix=2, WorkerFFTSizes = [8],
                        threadblock_bs=[1], threadblock_bs_dim=[0], data_type='double2'):
        self.fft_code = []
        self.data_type = data_type
        self.gPtr = "gPtr"
        self.rPtr = "rPtr"
        self.rPtr_2 = "rPtr_2"
        self.shPtr = "shPtr"
        self.WorkerFFTSizes = WorkerFFTSizes
        self.threadblock_bs = threadblock_bs
        self.threadblock_bs_dim = threadblock_bs_dim
        self.global_tensor_shape = global_tensor_shape
        self.radix = radix
        self.state_vec = th.zeros(64, 6)
        self.data = th.rand(self.WorkerFFTSizes, dtype=th.cfloat)
        self.output = th.zeros(self.WorkerFFTSizes, dtype=th.cfloat)
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
            "k" : ("int", "0"),
            "bx": ("int", "blockIdx.x"),
            "tx": ("int", "threadIdx.x"),
            "offset": ("int", "0"),
            self.gPtr: (f"{self.data_type}*", "inputs"),
            self.shPtr: (f"{self.data_type}*", "shared"),
            f"{self.rPtr}[{self.WorkerFFTSizes[dim]}]": (self.data_type, None),
            f"{self.rPtr_2}[{self.WorkerFFTSizes[dim]}]": (self.data_type, None),
            "tmp": (self.data_type, None),
            "angle": (self.data_type, None),
        }


    def save_generated_code(self, ):
        N = th.prod(th.as_tensor(self.global_tensor_shape[:-1]))
        for i in range(len(self.global_tensor_shape) - 1):
            file_name = f"../generated/fft_radix_{self.radix}_logN_{int(log(N, 2))}_upload_{i}.cuh"
            with open(file_name, 'w') as f:
                f.write(self.fft_code[i])

    def codegen(self,):

        reg_tensor_stride = th.as_tensor([1, 2, 4, 8, 16, 32, 64], dtype=th.float)
        state_vec = self.state_vec.clone()
        for dim in range(len(self.global_tensor_shape) - 2, -1, -1):
            self.init(dim)
            fft_code = self.head(len(self.global_tensor_shape) - 2 - dim)

            threadblock_tensor_shape =  self.threadblock_tensor_shape[dim]
            threadblock_bs = self.threadblock_bs[dim]
            threadblock_bs_dim = self.threadblock_bs_dim[dim]
            WorkerFFTSize = self.WorkerFFTSizes[dim]
            logWorkerFFTSize = int(log(WorkerFFTSize, 2))
            global_tensor_shape = self.global_tensor_shape
            blockorder = [i for i in range(len(global_tensor_shape))]
            fft_code += self.globalAccess(dim, global_tensor_shape, 
                                    threadblock_bs_dim, threadblock_bs, WorkerFFTSize,
                                     blockorder)
            for threadblock_dim in range(len(threadblock_tensor_shape)):
                self.state_vec = state_vec[:, :logWorkerFFTSize]
                if threadblock_dim != 0:
                    fft_code += self.shared2reg(threadblock_bs, threadblock_tensor_shape,
                                            WorkerFFTSize)
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
            if dim == 0:
                blockorder = self.list_reverse(blockorder, 0, -1)
                global_tensor_shape = self.list_reverse(global_tensor_shape, 0, -1)

            fft_code += self.globalAccess(blockorder[dim], global_tensor_shape, 
                        blockorder[threadblock_bs_dim], threadblock_bs, WorkerFFTSize,
                            blockorder, if_output=True, dict_output=dict_output, if_twiddle=(dim!=0))

            fft_code += self.epilogue()
            self.fft_code.append(fft_code)

    def head(self, dim):
        N = th.prod(th.as_tensor(self.global_tensor_shape[:-1]))
        head = f'''extern __shared__ {self.data_type} shared[];
__global__ void fft_radix_{self.radix}_logN_{int(log(N, self.radix))}_dim_{dim}''' \
        + f'''({self.data_type}* inputs, {self.data_type}* outputs, {self.data_type}* twiddle, int BS)''' + ''' {
    '''
        for key in self.local_variable.keys():
            head += f'''{self.local_variable[key][0]} {key};
    '''
        for key in self.local_variable.keys():
            if self.local_variable[key][1] is not None:
                head += f'''{key} = {self.local_variable[key][1]};
    '''
        return head
    
    def epilogue(self, ):
        epilogue = '''
}
'''
        return epilogue

    def globalAccess(self, dim, global_tensor_shape, threadblock_bs_dim, threadblock_bs, 
                    WorkerFFTSize, blockorder, if_output=False, dict_output=None, if_twiddle=False):
        global_tensor_shape = th.as_tensor(global_tensor_shape)
        threadblock_tensor_shape = th.ones_like(global_tensor_shape)
        threadblock_tensor_shape[dim] = global_tensor_shape[dim]
        threadblock_tensor_shape[threadblock_bs_dim] = threadblock_bs

        globalAccess_code = f'''bx = blockIdx.x;
    tx = threadIdx.x;
    ''' 
        if if_output is False:
            globalAccess_code += f'''{self.gPtr} = {self.local_variable[self.gPtr][1]};
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
                globalAccess_code += f'''
    {self.gPtr} += tx / {threadblock_bs} * {stride};
    '''
                access_stride = global_tensor_shape[i] // WorkerFFTSize * stride
            globalAccess_code += f'''
    {self.gPtr} += (bx % {global_tensor_shape[i] // threadblock_tensor_shape[i]}) * {threadblock_tensor_shape[i]} * {stride};
    bx = bx / {global_tensor_shape[i] // threadblock_tensor_shape[i]};
    '''
            if i == threadblock_bs_dim:
                globalAccess_code += f'''
    {self.gPtr} += tx % {threadblock_bs} * {stride};
    '''
            stride *= global_tensor_shape[i]
        if if_twiddle:
            globalAccess_code += f'''
    global_k += tx / {threadblock_bs};
    '''
        globalAccess_code += f'''
    {self.gPtr} += (bx % BS * {th.prod(global_tensor_shape[:-1])});
    '''
        
        for i in range(WorkerFFTSize):
            if not if_output:
                globalAccess_code += f'''
    {self.rPtr}[{i}] = *({self.gPtr} + {i * access_stride});
    '''
            else:
                if if_twiddle:
                    N = th.prod(global_tensor_shape[:(dim + 1)])
                    globalAccess_code += f'''
    // angle.x = cos(-2 * M_PI * global_j * (global_k + {i * global_tensor_shape[i] // WorkerFFTSize}) / {N});
    // angle.y = sin(-2 * M_PI * global_j * (global_k + {i * global_tensor_shape[i] // WorkerFFTSize}) / {N});
    angle = twiddle[{N - 1} + global_j * (global_k + {i * global_tensor_shape[i] // WorkerFFTSize})];
    tmp = {self.rPtr}[{dict_output[output_id]}];
    turboFFT_ZMUL({self.rPtr}[{dict_output[output_id]}], tmp, angle);
    '''
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
            
    def shared2reg(self, threadblock_bs, threadblock_tensor_shape, WorkerFFTSize):
        shared2reg_code  = ''''''
        access_stride = int(threadblock_bs * th.prod(th.as_tensor(threadblock_tensor_shape))
                         / WorkerFFTSize)
        shared2reg_code += f'''
    offset = 0;
    offset += tx;
    '''
        
        shared2reg_code += '''
    __syncthreads();
    '''
        for i in range(WorkerFFTSize):
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
        reg2shared_code = '''
    j = 0;
    offset  = 0;
    '''

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
            print(output_id, dict_output[output_id])
            if dim != len(threadblock_tensor_shape) - 1:
                reg2shared_code += f'''
    // angle.x = cos(M_PI * {-2 * output_id / N} * j);
    // angle.y = sin(M_PI * {-2 * output_id / N} * j);
    angle = twiddle[{N - 1} + {output_id} * j];
    tmp = {self.rPtr}[{dict_output[output_id]}];
    turboFFT_ZMUL({self.rPtr}[{dict_output[output_id]}], tmp, angle);
    '''
            if dim == 0:
                reg2shared_code += f'''
    {self.rPtr_2}[{output_id}] = {self.rPtr}[{dict_output[output_id]}];
    '''
            else:
                reg2shared_code += f'''
    {self.shPtr}[offset + {access_stride * output_id}] = {self.rPtr}[{dict_output[output_id]}];
    '''
        if dim == 0:
            for output_id in range(WorkerFFTSize):
                reg2shared_code += f'''
        k = (({output_id} + (threadIdx.x / {threadblock_bs}) % {WorkerFFTSize}) % {WorkerFFTSize});
        {self.shPtr}[offset + {access_stride} * k] = {self.rPtr_2}[k];
        '''
        return reg2shared_code
    
    def fft_reg(self, threadblock_bs, threadblock_tensor_shape, WorkerFFTSize, dim, reg_tensor_stride):
        fft_reg_code = ''''''
        
        bs = WorkerFFTSize // threadblock_tensor_shape[-1]
        logbs = int(log(bs, 2))

        logWorkerFFTSize = int(log(WorkerFFTSize, 2))
        st = logWorkerFFTSize - 1
    
        
        # Skip the leading batch size
        for i in range(st, logbs - 1, -1):
            print("*************************")
            print(reg_tensor_stride[i])
            for j in range(WorkerFFTSize):
                if self.state_vec[j, i] == 1:
                    continue
                id_j1 = int(th.dot(self.state_vec[j], reg_tensor_stride))
                id_j2 = int(id_j1 + reg_tensor_stride[i])
                id_k = int(th.dot(self.state_vec[j, logbs:i], reg_tensor_stride[:i - logbs]))                
                # print(id_k, self.state_vec[j, logbs:i], reg_tensor_stride[:i - logbs], reg_tensor_stride, logbs, i)
                print(id_j1, id_j2, id_k, i,  (2 ** (i + 1 - logbs)), st, logbs)
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
                    fft_reg_code += f'''
    angle.x = {cos(tmp_angle)};
    angle.y = {sin(tmp_angle)};
    turboFFT_ZMUL({self.rPtr}[{id_j2}], tmp, angle);
    '''

                tmp = self.data[id_j1].item()
                self.data[id_j1] = tmp + self.data[id_j2]
                self.data[id_j2] = tmp - self.data[id_j2]
                angle = cos(-2 * pi * id_k * 1 / (2 ** (i + 1 - logbs))) + sin(-2 * pi * id_k * 1 / (2 ** (i + 1 - logbs))) * 1.j
                self.data[id_j2] = self.data[id_j2] * angle
        
        return fft_reg_code

if __name__ == '__main__':
    params = []
    with open("../../param/param.csv", 'r') as file:
        for line in file:
            # Splitting each line by comma
            split_elements = line.strip().split(',')
            row = [int(element) for element in split_elements]
            params.append(row)
    
    # global_tensor_shape = [256, 256, 128, 1]
    for row in params:
        global_tensor_shape = [2 ** i for i in row[2:(2 + row[1])]]
        threadblock_bs = row[5:(5 + row[1])]
        WorkerFFTSizes = row[8:(8 + row[1])]
        threadblock_bs.reverse()
        global_tensor_shape.reverse()
        WorkerFFTSizes.reverse()
        global_tensor_shape.append(1)
        threadblock_bs_dim = [1, 0, 0]
        print(global_tensor_shape)
        print(WorkerFFTSizes)
        print(threadblock_bs)
        print(threadblock_bs_dim)
        fft = TurboFFT(global_tensor_shape=global_tensor_shape, WorkerFFTSizes=WorkerFFTSizes,
                        threadblock_bs=threadblock_bs, threadblock_bs_dim=threadblock_bs_dim[:row[1]])
        fft.codegen()
        fft.save_generated_code()