from multiprocessing import Process, Pool
import os
import subprocess

am_p_string = 'from multiprocessing import Process, Pool\n' \
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
              "    p.join()\n\n\n"


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
    file_string = am_p_string
    file_string += 'root_names = []\n' \
                   'coordinate_files = []\n' \
                   'for i in range(' + str(begin_index) + ', ' + str(end_index) + '):\n' \
                   "    coordinate_files.append('ground_snap.'+str(i))\n" \
                   "    root_names.append('"+root_name+"'+str(i))\n" \
                   "run_amber_parallel(False, root_names, coordinate_files, number_processors=" + \
                   str(n_processors_per_node)+')\n\n'

    open('hpc_traj_'+str(begin_index)+'.py', 'w').write(file_string)


def moab_header(id, walltime):
    # Provide walltime in seconds
    hrs = int(walltime / 3600)
    mins = int((walltime % 3600) / 60)
    secs = int(walltime % 60)
    walltime = '{:02d}:{:02d}:{:02d}'.format(hrs, mins, secs)
    working_directory = os.getcwd()
    job_script = '#!/bin/bash\n' \
                 '#MSUB -N traj_sub_' + str(id) + '.out\n' \
                 '#MSUB -j oe\n' \
                 '#MSUB -V\n' \
                 '#MSUB -o traj_sub_' + str(id) + '.out.stdout\n' \
                 '#MSUB -l nodes=1:ppn=1\n' \
                 '#MSUB -l walltime='+walltime+'\n\n' \
                 'module load intel/16.0.3 mkl/11.3.3 python\n' \
                 'source /users/dtracy/.bashrc\n' \
                 'cd ' + working_directory + '\n\n'
    return job_script


def submit_job_script(id, begin_index, end_index, root_name):
    job_script = moab_header(id, 7200)
    job_script += 'for index in {' + str(begin_index) + '..' + str(end_index) + "}\n" \
                  'do\n' \
                  '  $AMBERHOME/bin/sander -O  -i md_qmmm_amb.in -o '+root_name+'$index.out -r ' \
                  + root_name + '$index.rst -p m1.prmtop -x '+root_name+'$index.nc -c ground_snap.$index &\n' \
                  'done\n' \
                  'wait\n'
    script_file = 'hpc_traj_' + str(id) + '.sh'
    open(script_file, 'w').write(job_script)
    subprocess.run(['msub', script_file])


def run_hpc_trajectories(n_trajectories, n_processor_per_node, root_name):
    # root_name is the root name of the output file restart file etc. not the input
    # the input coordinates will be coming from ground_snap
    n_full_submissions = int(n_trajectories / n_processor_per_node)
    extra_submissions = int(n_trajectories % n_processor_per_node)
    job_id = 0
    for i in range(1, n_full_submissions+1):
        job_id += 1
        end_index = n_trajectories - extra_submissions
        begin_index = end_index - i * n_processor_per_node + 1
        submit_job_script(job_id, begin_index, end_index, root_name)
    if extra_submissions != 0:
        job_id += 1
        end_index = n_trajectories
        begin_index = end_index - extra_submissions
        submit_job_script(job_id, begin_index, end_index, root_name)


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



