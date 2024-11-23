#!/bin/bash
#SBATCH --job-name=check_gpu_info
#SBATCH --output=gpu_info_%j.out
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1

# Load the necessary module if required
# module load cuda/11.0  # Example for CUDA

# Get the GPU information
nvidia-smi --query-gpu=gpu_name --format=csv,noheader > gpu_name.txt

# Print the GPU information
cat gpu_name.txt
