import os
print("Be careful - running the entire ensemble submits 3900 jobs to the cluster")



file = "start_ensemble/latin_sh_save_runfile.txt"

num_lines = sum(1 for line in open(file))
#num_lines = 11


print("Reschedule file: latin_sh_save_runfile.txt")
die


with open(file) as fp:
    for cnt in range(0, num_lines):
        line = fp.readline()
        print(line)

        #iniate job script
        with open("job_submit.sh", "w+") as fh:
            fh.writelines("#!/bin/bash\n\n")

            #specifications of the job that should be submitted
            fh.writelines("#SBATCH --qos=short\n")
            fh.writelines("#SBATCH --job-name=socio_climate\n")
            fh.writelines("#SBATCH --account=dominoes\n\n")

            fh.writelines("#SBATCH --workdir=/p/projects/dominoes/nicowun/conceptual_tipping/uniform_distribution/socio_climate/out_err_files\n")
            fh.writelines("#SBATCH --output=outfile-%j.txt\n")
            fh.writelines("#SBATCH --error=error-%j.txt\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --ntasks-per-node=1\n")
            fh.writelines("#SBATCH --time=0-23:50:00\n\n")

            fh.writelines("module load anaconda/5.0.0_py3\n")
            fh.writelines("source activate earth_cascades\n\n")

            #job to be submitted
            fh.writelines("{}".format(line))
            fh.close()

        os.system("sbatch {}".format("job_submit.sh"))

