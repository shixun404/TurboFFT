import torch as th
import math
import random
a = th.rand(2000, 1024, dtype=th.cfloat)
e = th.ones(1024, dtype=th.cfloat)
b = a.clone()
# b[:1000, 0] += 1

import struct

# The float we'll be working with
for i in range(1000):
    f = b[i, 0].imag
    bit_to_flip = random.randint(0, 100) % 32
    # Convert the float to its IEEE 754 binary representation (as an integer)
    bytes_of_f = struct.pack('>f', f)
    int_representation = int.from_bytes(bytes_of_f, byteorder='big')

    # Flip the least significant bit
    flipped_int_representation = int_representation ^ (1 << bit_to_flip)

    # Convert the modified binary back to a float
    bytes_flipped = flipped_int_representation.to_bytes(4, byteorder='big')
    float_flipped = struct.unpack('>f', bytes_flipped)[0]
    if abs(float_flipped - f) / abs(f) > 1e-7:
        b[i, 0] = float_flipped + b[i,0].imag
    # print(float_flipped)
    # assert 0



for i in range(1024):
    angle = -1 * (i % 3) * math.pi * 2.0 / 3.0
    # e[i] = math.cos(angle) + math.sin(angle) * 1.j
    e[i] = 1 + 1.j
eW = th.fft.fft(e)
Wa = th.fft.fft(a)
Wb = th.fft.fft(b)
sys_err = th.zeros(2000)
inj_err = th.zeros(2000)
for i in range(2000):
    eW_a = th.dot(eW, a[i])
    e_Wa = th.dot(e, Wa[i])
    e_Wb = th.dot(e, Wb[i])
    sys_err[i] = abs(e_Wa.item() - eW_a.item()) / abs(eW_a)
    inj_err[i] = abs(e_Wb.item() - eW_a.item()) / abs(eW_a)
# err_threshold = inj_err.max().item()
err_threshold = 1e5
n_inj = []
n_sys = []
err_list = []
while err_threshold > 0:
    # print(err_threshold)
    if err_threshold < 1e-8:
        err_threshold = 0
    elements_larger_than_half = inj_err[:1000] > err_threshold
    # Count the number of True values in the comparison result
    count = elements_larger_than_half.sum().item()  # .item() t
    n_inj.append(count)
    elements_larger_than_half = sys_err > err_threshold
    # Count the number of True values in the comparison result
    count = elements_larger_than_half.sum().item()  # .item() t
    n_sys.append(count)
    err_list.append(err_threshold)
    err_threshold /= 2
n_sys = th.as_tensor(n_sys) / 2000
n_inj = th.as_tensor(n_inj) / 1000
err_list = th.as_tensor(err_list)

import matplotlib.pyplot as plt
plt.rcParams['lines.linewidth'] = 3
# plt.rcParams["font.family"] = "Times New Roman"
plt.rc('font', size=30, weight='bold')
# plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["hatch.color"] = 'white'
plt.rcParams['hatch.linewidth'] = 2.0
fig, ax = plt.subplots(ncols=2, figsize =(16, 4))
ax[0].grid(linestyle='--', linewidth=1.5, zorder=0)
ax[1].grid(linestyle='--', linewidth=1.5, zorder=0)
ax[0].plot(n_sys, n_inj, marker='^', ms=6,label='ROC FP32')

ax[1].plot(err_list,n_sys, color='green', linestyle='-', label='False Alarm FP32')
ax[1].plot(err_list,n_inj, color='red', linestyle='-', label='Fault Detection FP32')
ax[1].set_xscale("log", base=10)

# ax[1].set_xlim([0, 10])
# ax[1].legend(loc = 'upper right')
# ax[0].legend(loc = 'upper right')


a = th.rand(2000, 1024, dtype=th.cdouble)
e = th.ones(1024, dtype=th.cdouble)
b = a.clone()
# b[:1000, 0] += 1

import struct

# The float we'll be working with
for i in range(1000):
    f = b[i, 0].real
    bit_to_flip = random.randint(0, 100) % 64
    # Convert the float to its IEEE 754 binary representation (as an integer)
    bytes_of_f = struct.pack('>d', f)
    int_representation = int.from_bytes(bytes_of_f, byteorder='big')

    # Flip the least significant bit
    flipped_int_representation = int_representation ^ (1 << bit_to_flip)

    # Convert the modified binary back to a float
    bytes_flipped = flipped_int_representation.to_bytes(8, byteorder='big')
    float_flipped = struct.unpack('>d', bytes_flipped)[0]
    if abs(float_flipped - f) / abs(f) > 1e-17:
        # b[i, 0] = float_flipped + b[i,0].imag
        b[i, 0] = float_flipped *1.j + b[i,0].real
# print("adsadasd")
for i in range(1024):
    angle = -1 * (i % 3) * math.pi * 2.0 / 3.0
    # e[i] = math.cos(angle) + math.sin(angle) * 1.j
    e[i] = 1 + 1.j
eW = th.fft.fft(e)
Wa = th.fft.fft(a)
Wb = th.fft.fft(b)
sys_err = th.zeros(2000)
inj_err = th.zeros(2000)
for i in range(2000):
    # print(i)
    eW_a = th.dot(eW, a[i])
    e_Wa = th.dot(e, Wa[i])
    e_Wb = th.dot(e, Wb[i])
    sys_err[i] = abs(e_Wa.item() - eW_a.item()) / abs(eW_a)
    inj_err[i] = abs(e_Wb.item() - eW_a.item()) / abs(eW_a)
err_threshold = 1e5
n_inj = []
n_sys = []
err_list = []
while err_threshold > 0:
    # print(err_threshold)
    if err_threshold < 1e-18:
        err_threshold = 0
    elements_larger_than_half = inj_err[:1000] > err_threshold
    # Count the number of True values in the comparison result
    count = elements_larger_than_half.sum().item()  # .item() t
    n_inj.append(count)
    elements_larger_than_half = sys_err > err_threshold
    # Count the number of True values in the comparison result
    count = elements_larger_than_half.sum().item()  # .item() t
    n_sys.append(count)
    err_list.append(err_threshold)
    err_threshold /= 2
n_sys = th.as_tensor(n_sys) / 2000
n_inj = th.as_tensor(n_inj) / 1000
err_list = th.as_tensor(err_list)

# ax[0].grid(linestyle='--', linewidth=1.5, zorder=0)
# ax[1].grid(linestyle='--', linewidth=1.5, zorder=0)
ax[0].plot(n_sys, n_inj, marker='^', ms=6, label='ROC FP64')
ax[0].set_xlabel('(a) False Alarms%')
ax[0].set_ylabel('Fault Detection%')
ax[1].plot(err_list,n_sys, color='green', linestyle='--', label='False Alarm FP64')
ax[1].plot(err_list,n_inj, color='red', linestyle='--', label='Fault Detection FP64')
ax[1].set_xscale("log", base=10)
ax[1].set_ylabel('Rate%')
ax[1].set_xlabel('(b) Relative Error Threshold')
ax[1].set_xlim([0, 10])

ax[0].legend(loc = 'lower right')
ax[1].legend(loc = 'lower left', fontsize=20)
# ax[0].set_yticks(rotation=90)
ax[0].yaxis.set_tick_params(rotation=90)
ax[1].yaxis.set_tick_params(rotation=90)
# ax[1].set_yticks(rotation=90)
fig.subplots_adjust(hspace=0.1, wspace = 0.16)
plt.savefig("../artifact_figures/figure19.pdf", bbox_inches='tight')
print("./artifact_figures/figure19.pdf")

