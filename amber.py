from multiprocessing import Process, Pool
import subprocess


def run_amber(conjoined_list):
    cd_file = 'm1.inpcrd'
    root_name = conjoined_list[0]
    coordinate_file = conjoined_list[1]
    if coordinate_file is not None:
        cd_file = coordinate_file
    subprocess.run(['sander', '-O', '-i', 'md_qmmm_amb.in', '-o', root_name+'.out', '-c', cd_file, '-p',
                    'm1.prmtop', '-r', root_name+'.rst', '-x', root_name+'.nc'])


def create_slurm_script(script_file_name, n_arrays, root_name, crd_file='m1.inpcrd'):
    job_script = '#!/bin/bash\n' \
                 '#SBATCH --job-name=NAESMD\n' \
                 '#SBATCH --output my_job\n' \
                 '#SBATCH --error my_job-#j.err #Error File\n' \
                 '#SBATCH --mail-type=FAIL,END #What emails to send\n' \
                 '#SBATCH --nodes=1 #No. computers requested\n' \
                 '#SBATCH --mem-per-cpu=2000mb #Per processor memory requested\n' \
                 '#SBATCH --array=1-' + str(n_arrays) + '\n' \
                 '#SBATCH --time=01:00:00 #Walltime\n\n' \
                 '$AMBERHOME/bin/sander -O -i md_qmmm_amb.in -o ' + root_name + '$SLURM_ARRAY_TASK_ID.out -c ' \
                 + crd_file + '.$SLURM_ARRAY_TASK_ID' + ' -p m1.prmtop -r ' + root_name + \
                 '$SLURM_ARRAY_TASK_ID.rst -x ' \
                 + root_name + '$SLURM_ARRAY_TASK_ID.nc\n'
    open(script_file_name, 'w').write(job_script)

# Returns a list of type gaussian output
def run_amber_parallel(root_names, coordinate_files, number_processors=4):
    p = Pool(number_processors)
    conjoined_list = []
    for i in range(len(root_names)):
        conjoined_list.append([root_names[i], coordinate_files[i]])
    p.map(run_amber, conjoined_list)
    p.close()
    p.join()

