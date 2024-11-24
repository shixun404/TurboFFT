import torch as th
def load_data_single(file_name):
    data = th.zeros(26, 29, 3)
    data_list = []
    with open(file_name, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split(',')
            data_point = [line[0], int(line[1]), int(line[2]), 
                        float(line[3]), float(line[4]), float(line[5])]
            data_list.append(data_point)
            data[int(line[1]), int(line[2]), 0] = data_point[3]
            data[int(line[1]), int(line[2]), 1] = data_point[4]
            data[int(line[1]), int(line[2]), 2] = data_point[5]
            
    return data_list, data

def load_data_single_sc(file_name):
    data = th.zeros(26, 1, 3)
    data_list = []
    with open(file_name, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split(',')
            data_point = [line[0], int(line[1]), int(line[2]), 
                        float(line[3]), float(line[4]), float(line[5])]
            data_list.append(data_point)
            data[int(line[1]), 0, 0] = data_point[3]
            data[int(line[1]), 0, 1] = data_point[4]
            data[int(line[1]), 0, 2] = data_point[5]
            
    return data_list, data


def load_data(file_name):
    data = th.zeros(2, 26, 29, 3)
    data_list = []
    with open(file_name, 'r') as f:
        i = 0
        for line in f:
            
            line = line.strip()
            line = line.split(',')
            
            id = 0 if line[0] == 'turboFFT' else 1
            data_point = [line[0], int(line[1]), int(line[2]), 
                        float(line[3]), float(line[4]), float(line[5])]
            data_list.append(data_point)
            
            
            data[id, int(line[1]), int(line[2]), 0] = data_point[3]
            data[id, int(line[1]), int(line[2]), 1] = data_point[4]
            data[id, int(line[1]), int(line[2]), 2] = data_point[5]
            i += 1


    return data_list, data