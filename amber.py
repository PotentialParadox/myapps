from multiprocessing import Process, Pool
import subprocess


def run_amber(conjoined_list):
    cd_file = 'm1.inpcrd'
    root_name = conjoined_list[0]
    coordinate_file = conjoined_list[1]
    amber = conjoined_list[2]
    if coordinate_file is not None:
        cd_file = coordinate_file
    subprocess.run([amber, '-O', '-i', 'md_qmmm_amb.in', '-o', root_name+'.out', '-c', cd_file, '-p',
                    'm1.prmtop', '-r', root_name+'.rst', '-x', root_name+'.nc'])


def run_amber_parallel(pmemd_available, root_names, coordinate_files, number_processors=4):
    if pmemd_available:
        amber = 'pmemd'
    else:
        amber = 'sander'
    p = Pool(number_processors)
    conjoined_list = []
    for i in range(len(root_names)):
        conjoined_list.append([root_names[i], coordinate_files[i], amber])
    p.map(run_amber, conjoined_list)
    p.close()
    p.join()


def create_hpc_python_file(begin_index, end_index, n_processors_per_node, root_name):
    file_string = 'from multiprocessing import Process, Pool\n' \
                  'import subprocess\n\n\n' \
                  'def run_amber(conjoined_list):\n' \
                  "    cd_file = 'm1.inpcrd'\n" \
                  '    root_name = conjoined_list[0]\n' \
                  '    coordinate_file = conjoined_list[1]\n' \
                  '    amber = conjoined_list[2]\n' \
                  '    if coordinate_file is not None:\n' \
                  '        cd_file = coordinate_file\n' \
                  "    subprocess.run([amber, '-O', '-i', 'md_qmmm_amb.in', '-o', root_name+'.out', '-c'," \
                  " cd_file, '-p',\n" \
                  "                    'm1.prmtop', '-r', root_name+'.rst', '-x', root_name+'.nc'])\n\n\n" \
                  'def run_amber_parallel(pmemd_available, root_names, coordinate_files, number_processors=4):\n' \
                  '    if pmemd_available:\n' \
                  "        amber = 'pmemd'\n" \
                  "    else:\n" \
                  "        amber = 'sander'\n" \
                  "    p = Pool(number_processors)\n" \
                  "    conjoined_list = []\n" \
                  "    for i in range(len(root_names)):\n" \
                  "        conjoined_list.append([root_names[i], coordinate_files[i], amber])\n" \
                  "    p.map(run_amber, conjoined_list)\n" \
                  "    p.close()\n" \
                  "    p.join()\n\n\n" \
                  'root_names = []\n' \
                  'coordinate_files = []\n' \
                  'for i in range(' + str(begin_index) + ', ' + str(end_index) + '):\n' \
                  "    coordinate_files.append('ground_snap.'+str(i))\n" \
                  "    root_names.append('nasqm_abs_'+str(i))\n" \
                  "run_amber_parallel(False, root_names, coordinate_files, number_processors=" + \
                  str(n_processors_per_node)

    open('hpc_traj_'+str(begin_index)+'.py', 'w').write(file_string)


def run_hpc_trajectories(n_trajectories, n_processor_per_node, root_name):
    # root_name is the root name of the output file restart file etc. not the input
    # the input coordinates will be coming from ground_snap
    n_full_submissions = int(n_trajectories / n_processor_per_node)
    extra_submissions = int(n_trajectories % n_processor_per_node)
    for i in range(n_full_submissions):
        begin_index = i * n_processor_per_node + 1
        end_index = begin_index + n_processor_per_node
        create_hpc_python_file(begin_index, end_index, n_processor_per_node, root_name)
    if extra_submissions != 0:
        end_index = n_trajectories + 1
        begin_index = end_index - extra_submissions
        create_hpc_python_file(begin_index, end_index, n_processor_per_node, root_name)

def create_snapshot_slurm_script(script_file_name, n_trajectories, n_frames, root_name, crd_file='m1.inpcrd'):
    n_arrays = n_trajectories
    job_script = '#!/bin/bash\n' \
                 '#SBATCH --job-name=NAESMD\n' \
                 '#SBATCH --output my_job\n' \
                 '#SBATCH --error my_job-#j.err #Error File\n' \
                 '#SBATCH --mail-type=FAIL,END #What emails to send\n' \
                 '#SBATCH --nodes=1 #No. computers requested\n' \
                 '#SBATCH --mem-per-cpu=2000mb #Per processor memory requested\n' \
                 '#SBATCH --array=1-' + str(n_arrays) + '\n' \
                 '#SBATCH --time=01:00:00 #Walltime\n' \
                 '\n' \
                 'for i in `seq 1 ' + str(n_frames) + '`;\n' \
                 'do\n' \
                 '    $AMBERHOME/bin/sander -O -i md_qmmm_amb.in -o ' + root_name + '$SLURM_ARRAY_TASK_ID_$i.out -c ' \
                 + crd_file + '$SLURM_ARRAY_TASK_ID.' + '$i -p m1.prmtop -r ' + root_name + \
                 '$SLURM_ARRAY_TASK_ID_$i.rst -x ' + root_name + '$SLURM_ARRAY_TASK_ID_$i.nc\n' \
                 'done'
    open(script_file_name, 'w').write(job_script)



