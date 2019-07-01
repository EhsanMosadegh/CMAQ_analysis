#!/bin/csh -f 

#$ -o ../logs/log.co.scen5.sep.txt	 			 ##$ output file name

#$ -N s5sep						##$ name of my job
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
 
echo "CO job =====> submitted"
date 			

python -u ../run_scripts/spatial_CMAQ_analysis_scen5_sep.py

echo "CO job =====> ended"
date
