#! /bin/bash
#SBATCH -J A2_chose6
#SBATCH -p gpu_l40,gpu_l48
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --gres=gpu:1

sys_name="A2"
temp=300
mkdir $sys_name
mkdir $sys_name/$temp
python patched_mpipi.py --name $sys_name \
	--temp $temp \
	--time 10 \
	--fasta_file proteins.fasta \
	--phos 32,42,45,60,67,74,81,85,88,93,98,104,111,129,134,145,151 \
	--pH 5.5 \
	--ionic 0.025 \

