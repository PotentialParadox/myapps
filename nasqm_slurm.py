'''
Functions used to interface the NASQM automation routines with
the Slurm wrapper
'''
import math
from slurm import Slurm

def create_slurm_header(user_input):
    '''
    Returns the slurm header dictionary
    '''
    return {'email': user_input.email, 'email_options': user_input.email_options,
            'n_nodes': user_input.number_nodes, 'ppn': user_input.processors_per_node,
            'memory': user_input.memory_per_node, 'walltime': user_input.walltime,
            'max_jobs': user_input.max_jobs}


def run_hpc_trajectories(user_input, restart_root, output_root, title, n_trajectories):
    '''
    Run multiple sander trajectories over hpc
    '''
    n_arrays = math.ceil(n_trajectories/user_input.processors_per_node)
    slurm_header = create_slurm_header(user_input)
    command = "module load intel/2016.0.109\n\n"
    command += "for i in {1.."+n_arrays+"}\n" \
               "do\n" \
               '    multiplier="$((${SLURM_ARRAY_TASK_ID}))' \
    command += "$AMBERHOME/bin/sander - -i md_qmmm_amb.in -o "+output_root\
               +"${SLURM_ARRAY_TASK_ID}.out -c " \
              +restart_root+".${SLURM_ARRAY_TASK_ID} -p m1.prmtop -r " \
              +output_root+"${SLURM_ARRAY_TASK_ID}.rst -x "+output_root+"${SLURM_ARRAY_TASK_ID}.nc"
    slurm_script = Slurm(slurm_header)
    slurm_file = slurm_script.create_slurm_script(command, title, n_arrays)
    print(slurm_file)
