import os
import subprocess
import math
import torch as th
import re

def test():
    filename = 'vkfft_output'
    file = open(filename, 'w')
    cnt = 1
    print(f"The default VkFFT executable only provides functionality to test a single data point. Therefore, we benchmarked all test points from Fig. 1, with 400 data points each for float and double precision C2C FFTs, namely 800 datapoints in total.")
    for precision in range(2):
        for x in range(1, 26):
            X_value = 2**x
            
            # Calculate the upper bound for B
            max_B_exponent = 28 - x
            B_upper_bound = 2**max_B_exponent
            
            # Iterate over B from 2^0 to 2^(28-X)
            for b in range(max_B_exponent + 1):
                B_value = 2**b
                
                # Formulate the command
                command = f"./VkFFT_TestSuite -benchmark_vkfft -X {X_value} -B {B_value} -P {precision}"
                
                # Execute the command
                process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                file.write(f"Command: {command}\n")
                file.write("STDOUT:\n" + stdout.decode())
                if stderr:
                    file.write("STDERR:\n" + stderr.decode())
                file.write("\n\n")  # Adding some spacing between entries for readability
                # Print output and error if any
                print(f"Benchmark VkFFT {cnt}/800")
                cnt += 1
                print("STDOUT:", stdout.decode())
                if stderr:
                    print("STDERR:", stderr.decode())
    print("Benchmarking VkFFT Finished")

if __name__ == "__main__":
    test()
    filename = "vkfft_output"
    f = open(filename, 'r')
    lines = f.readlines()
    x = 0
    bs = 0
    p = 0
    exec_time = 0
    result = th.ones(2, 30, 30)
    for line in lines:
        if 'Command' in line:
            x = int(math.log2(float(line.split(' ')[-5])))
            bs = int(math.log2(float(line.split(' ')[-3])))
            p = int(line.split(' ')[-1])
        if not(x + bs <= 28 and x <= 25 and x >= 1):
            continue
        if 'VkFFT System:' in line:
            match = re.search(r"avg_time_per_step:\s*([0-9.]+)\s*ms", line)
            if match:
                exec_time =  float(match.group(1))  # Convert the extracted value to float
                # print(exec_time, type(exec_time))
                result[p, x, bs] = exec_time / 2
            else:
                print(x, bs, p, 'not find') # Return None if the pattern is not found

    th.save(result, 'vkFFT.pt')