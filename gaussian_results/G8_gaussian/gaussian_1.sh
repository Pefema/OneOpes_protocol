#!/bin/env bash
#SBATCH -J GausG8
#SBATCH -e gaus_%j.e
#SBATCH -o gaus_%j.o
#SBATCH --nodes=1
#SBATCH --constraint="V4|V5|V6|V9"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH -p shared-cpu
#SBATCH --mem=10G

ml gaussian

g16 < G8.com > G8.log