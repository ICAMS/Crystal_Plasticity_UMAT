# Workflow Kanapy2RVE for fatigue simulations with Abaqus

# 
# On jupyter.icams.rub.de:
module load anaconda/2019
conda activate knpy
generate_rve.ipynb -> {geometry}.inp {material}.inp

#
# On cip-sim{n}:
module load anaconda/latest
module load abaqus/latest
cp -r abq_template {job_directory}
cp ../../source/* {job_directory}
cp {geometry}.inp {job_directory}
cp {material}.inp {job_directory}/Material.inp
mv {job_directory} /scratch/{user}
cd /scratch/{user}/{job_directory}
python Create_PeriodicBC_EDGE_3D.py {geometry}.inp
abaqus job={jobname} inp=Fatigue_BC_Amp_05.inp user=umat.f &

#
# After Abaqus completed:
abaqus python analyze_odb_structured_CP.py {jobname} -> {jobname}_sig_eps.csv
