# **TurboFFT: Co-Designed High-Performance and Fault-Tolerant Fast Fourier Transform on GPUs**

## Table of Contents
- [I. Overview of Contributions and Artifacts](#i-overview-of-contributions-and-artifacts)
  - [A. Paper's Main Contributions](#a-papers-main-contributions)
  - [B. Computational Artifacts](#b-computational-artifacts)
- [II. Introduction](#ii-introduction)
- [III. Getting Started](#iii-getting-started)
  - [A. Extract Code Repositories](#a-extract-code-repositories)
  - [B. Install Host Machine Compilation Prerequisites](#b-install-host-machine-compilation-prerequisites)
  - [C. Install Host Machine Codegen & Plot Prerequisites](#c-install-host-machine-codegen--plot-prerequisites)
- [IV. Reproducing Results](#iv-reproducing-results)
  - [A. How to Run](#a-how-to-run)
  - [B. Workflow Overview](#b-workflow-overview)
  - [C. Runtime Details](#c-runtime-details)

---

## I. Overview of Contributions and Artifacts

### A. Paper's Main Contributions

1. **TurboFFT without fault tolerance** outperforms the popular open-source library VkFFT, and is comparable to the state-of-the-art closed-source library, cuFFT.
2. **TurboFFT with two-side fault tolerance** efficiently fuses the checksum computation into the FFT kernel, minimizing the fault tolerance overhead compared to existing offline fault-tolerant FFT (FT-FFT).
3. **TurboFFT’s online error correction** protects FFT computation on-the-fly, obtaining lower error correction overhead compared to the time-redundant recomputation in offline FT-FFT under error injections.

### B. Computational Artifacts

| **Artifact ID** | **Contributions** | **Related Supported Paper Elements** |
|-----------------|-------------------|--------------------------------------|
| **A1**          | C1               | Figure 1, 10–14, 21                 |
|                 | C2               | Figure 16–18, 19–20                 |
|                 | C3               | Figure 19–20, 22                    |

---

## II. Introduction

This document demonstrates how to reproduce the results in the paper:

**TurboFFT: Co-Designed High-Performance and Fault-Tolerant Fast Fourier Transform on GPUs.**

All supplementary files are available on Zenodo. The repository `PPoPP25_Artifact_TurboFFT.zip` consists of all code, including two one-command scripts `run_A100.sh` and `run_T4.sh` to reproduce all result figures listed in [Table I](#b-computational-artifacts).

We executed all benchmarks in the paper using the hardware detailed in [Table II](#hardware-environment) and the software detailed in [Table III](#software-environment).

---

## III. Getting Started

This section guides you through the necessary steps to set up your machine. Please follow these steps before starting to reproduce the results.

### A. Extract Code Repositories

1. Start by downloading `PPoPP25_Artifact_TurboFFT.zip`.
2. Extract the archive into an empty directory and change into this directory using the commands below. Make sure that the absolute path to this directory contains no spaces.

```bash
# Create a reproduce directory
mkdir reproduce
cd reproduce
cp <path-to>/PPoPP25_Artifact_TurboFFT.zip ./

# Extract the artifact
unzip <path-to>/PPoPP25_Artifact_TurboFFT.zip
cd PPoPP25_Artifact_TurboFFT
```

Now, the folder `reproduce` contains all the necessary code to produce the results shown in the paper. The code consists of:
- **TurboFFT**: Source code of the high-performance FFT library.
- **Common**: Contains CUDA helper functions from `NVIDIA/cuda-samples`.

### B. Install Host Machine Compilation Prerequisites

1. **GCC**: Install a recent GCC version (≥ 11.2.0).
2. **CMake**: Install a recent CMake version (≥ 3.24.3).
3. **CUDA Toolkit**: 
   - Use CUDA Toolkit 12.0 for an A100 machine.
   - Use CUDA Toolkit 11.6 for a T4 machine.

```bash
# Check the version after installing
gcc --version
cmake --version
nvcc --version
```

### C. Install Host Machine Codegen & Plot Prerequisites

1. **Python**: Use a recent Python version (≥ 3.9).
2. **PyTorch**: Use a recent version (≥ 2.4).
3. **NumPy**: Use a recent version (≥ 2.0.2).
4. **Matplotlib**: Use a recent version (≥ 3.8.4).
5. **Seaborn**: Use a recent version (≥ 0.13.2).

Below is an example of setting up a Python environment:

```bash
python -m venv .venv --prompt turbofft
source .venv/bin/activate
pip install --upgrade pip
pip install torch
pip install numpy
pip install matplotlib
pip install seaborn
```

---

## IV. Reproducing Results

The one-command scripts `run_A100.sh` and `run_T4.sh` allow you to regenerate:
- **11 experimental result figures** on NVIDIA A100 GPUs (Figures 1, 10–14, and 16–20)
- **2 experimental result figures** on NVIDIA T4 GPUs (Figures 21–22)

### A. How to Run

1. Ensure all dependencies are installed (see [Section III](#iii-getting-started)).
2. Run the script:
   - On an NVIDIA A100 machine:
     ```bash
     ./run_A100.sh
     ```
   - On an NVIDIA T4 machine:
     ```bash
     ./run_T4.sh
     ```
3. View results:
   - Experimental data will be available in the `artifact_data` directory.
   - Figures will be saved in the `artifact_figures` directory.

### B. Workflow Overview

The scripts `run_A100.sh` and `run_T4.sh` execute the following steps:

1. **Environment Setup**: Configures environment variables.
2. **Code Generation**: Generates required CUDA kernels.
3. **Compilation**: Builds TurboFFT and related binaries.
4. **Benchmarking**: Runs benchmarks for TurboFFT.
5. **Plotting**: Produces figures matching the paper’s results.

### C. Runtime Details

[Table IV](#estimated-execution-time) shows the estimated execution time of `run_A100.sh` and `run_T4.sh`.

- The script `run_A100.sh` takes approximately **2 hours** on a machine with an AMD EPYC 7763 64-Core Processor and an NVIDIA A100 40GB GPU.
- The script `run_T4.sh` takes approximately **30 minutes** on a machine with an Intel(R) Xeon(R) Silver 4216 CPU and an NVIDIA T4 GPU.

#### Hardware Environment

| System Type | Description                                                                     |
|-------------|---------------------------------------------------------------------------------|
| **System A** | **GPU**: 1× NVIDIA A100-SXM4-40GB<br>**GPU Power**: 400 W<br>**CPU**: AMD EPYC 7713 64-Core Processor<br>**Cores per socket**: 64<br>**Threads per core**: 2<br>**Memory**: 256 GB |
| **System B** | **GPU**: 1× NVIDIA Tesla-T4<br>**GPU Power**: 70 W<br>**CPU**: Intel(R) Xeon(R) Silver 4216 CPU<br>**Cores per socket**: 16<br>**Threads per core**: 2<br>**Memory**: 192 GB |

#### Software Environment

| System    | Software      | Version  |
|-----------|-------------- |--------- |
| **System A** | gcc          | 12.3.0   |
|           | cmake         | 3.24.3   |
|           | cudatoolkit   | 12.0     |
|           | python        | 3.10.14  |
|           | torch         | 2.5.1    |
|           | numpy         | 2.1.3    |
|           | matplotlib    | 3.8.4    |
|           | seaborn       | 0.13.2   |
| **System B** | gcc          | 11.2.0   |
|           | cmake         | 3.26.4   |
|           | cudatoolkit   | 11.6     |
|           | python        | 3.9.18   |
|           | torch         | 2.5.1    |
|           | numpy         | 2.0.2    |
|           | matplotlib    | 3.9.2    |
|           | seaborn       | 0.13.2   |

#### Estimated Execution Time

| Script        | CodeGen | TurboFFT  | Baseline (cuFFT + VkFFT) | Plot | Total  |
|---------------|--------|-----------|---------------------------|------|--------|
| **run_A100.sh** | 20 s   | 15 min    | 90 min                    | 3 min| 2 hr    |
| **run_T4.sh**   | 20 s   | 10 min    | 10 min                    | 3 min| 30 min  |