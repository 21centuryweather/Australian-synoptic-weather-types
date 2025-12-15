#PBS -N precip_data_preparation
#PBS -P if69
#PBS -l ncpus=4
#PBS -l mem=256GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/if69+gdata/zv2+gdata/fs38+gdata/xp65
#PBS -q normalsr
#PBS -o precip_data_preparation.out
#PBS -e precip_data_preparation.err

# Load modules
module use /g/data/xp65/public/modules
module load singularity
module load conda/analysis3

# Move to submission directory
cd $PBS_O_WORKDIR

# Run your python script
python precip_data_preparation.py
