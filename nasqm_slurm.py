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
    command += "for i in {1.."+str(n_arrays)+"}\n" \
               "do\n" \
               '    MULTIPLIER="$((${SLURM_ARRAY_TASK_ID} - 1))\n' \
               '    FIRST_COUNT="$((${SLURM_CPUS_ON_NODE} * ${MULTIPLIER}))\n' \
               '    ID="$((${FIRST_COUNT} + ${i}))\n'
    command += "    $AMBERHOME/bin/sander -i md_qmmm_amb.in -o " + output_root \
              + "${ID}.out -c " \
              + restart_root+".${ID} -p m1.prmtop -r " \
              + output_root+ "${ID}.rst -x "+output_root \
              +"${ID}.nc\n" \
              +"done"
    slurm_script = Slurm(slurm_header)
    slurm_file = slurm_script.create_slurm_script(command, title, n_arrays)
    print(slurm_file)
