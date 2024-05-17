#!/bin/bash
#SBATCH -p LM
#SBATCH -N 1
#SBATCH --mem 200GB
#SBATCH -t 10:00:00

#srun ./extractSequences.sh 5476 /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/CA/17120D-10-03_S3_L001_1.fastq /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/CA/17120D-10-03_S3_L001_1.fastq.csv /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Filtered/CA/filtered.17120D-10-03_S3_L001.R1 0.10 0.75

#srun ./extractSequences.sh 5476 /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/CA/17120D-10-03_S3_L001_2.fastq /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/CA/17120D-10-03_S3_L001_1.fastq.csv /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Filtered/CA/filtered.17120D-10-03_S3_L001.R2 0.10 0.75


#exit

./extractSequences.sh 1280 /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/SA1/17120FL-10-05_S4_L001_1.fastq /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/SA1/17120FL-10-05_S4_L001_1.fastq.csv /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Filtered/SA1/17120FL-10-05_S4_L001_1.R1 0.10 0.75

./extractSequences.sh 1280 /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/SA1/17120FL-10-05_S4_L001_2.fastq /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/SA1/17120FL-10-05_S4_L001_1.fastq.csv /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Filtered/SA1/17120FL-10-05_S4_L001_1.R2 0.10 0.75

#exit

./extractSequences.sh 1280 /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/SA2/17120FL-10-06_S5_L001_1.fastq /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/SA2/17120FL-10-06_S5_L001_1.fastq.csv /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Filtered/SA2/17120FL-10-06_S5_L001_1.R1 0.10 0.75

./extractSequences.sh 1280 /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/SA2/17120FL-10-06_S5_L001_2.fastq /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Unfiltered/SA2/17120FL-10-06_S5_L001_1.fastq.csv /pylon5/ms5phsp/rachid/Data/AssemblyProjects/RawData/Illumina/Filtered/SA2/17120FL-10-06_S5_L001_1.R2 0.10 0.75
