#!/bin/csh -f 

#$ -o log.co.main.txt					##$ output file name

#$ -N s1jul						##$ name of my job
#$ -S /bin/csh						##$ specify the shell
#$ -cwd							##$ job is submitted from here
#$ -V  							##$ uses current env variables / preserves your environment
#$ -q Students.q 
#$ -pe FillUp 1 					##$ set paralel environment to run parallel 
#$ -j y							##$ joint output and error into one file
#$ -M ehsan.mosadegh@dri.edu 
#$ -m abe


#source /etc/csh.cshrc					##$ to load the module for compiling
source /scratch/ehsanm/virtualEnv/virenvConda/bin/activate.csh virenvConda

module purge
module load intel/2015
module list
unlimit
limit

set mode='main'

echo "CO main job =====> submitted"
date 			

python -u spatial_CMAQ_analysis_original.py > log.cmaq.${mode}.txt 

echo "CO main job =====> ended"
date
