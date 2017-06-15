'''
Functions used to interface the NASQM automation routines with
the Slurm wrapper
'''
import math
from slurm import Slurm, run_slurm

def create_slurm_header(user_input):
    '''
    Returns the slurm header dictionary
    '''
    return {'email': user_input.email, 'email_options': user_input.email_options,
            'n_nodes': user_input.number_nodes, 'ppn': user_input.processors_per_node,
            'memory': user_input.memory_per_node, 'walltime': user_input.walltime,
            'max_jobs': user_input.max_jobs}

def build_command(amber, n_trajectories):
    '''
    Returns the command for the slurm script
    '''
    command = "module load intel/2016.0.109\n\n"
    command += "for i in $(seq 1 {})\n".format(n_trajectories)
    command += "do\n" \
               '    MULTIPLIER="$((${SLURM_ARRAY_TASK_ID} - 1))"\n' \
               '    FIRST_COUNT="$((${SLURM_CPUS_ON_NODE} * ${MULTIPLIER}))"\n' \
               '    ID="$((${FIRST_COUNT} + ${i}))"\n'
    command += "    $AMBERHOME/bin/sander -i {}${{ID}}.in -o {}".format(amber.input_roots[0],
                                                                        amber.output_roots[0])
    command += "${ID}.out -c "
    command += "{}.${{ID}} -p m1.prmtop -r ".format(amber.restart_roots[0])
    command += "{}${{ID}}.rst -x {}".format(amber.output_roots[0], amber.output_roots[0]) \
               +"${ID}.nc\n" \
               +"done\n" \
               +"wait\n"
    return command

def slurm_trajectory_files(user_input, amber, title, n_trajectories):
    '''
    Run multiple sander trajectories over hpc
    '''
    n_arrays_max = math.floor(n_trajectories/user_input.processors_per_node)
    n_trajectories_remaining = n_trajectories - n_arrays_max * user_input.processors_per_node
    slurm_header = create_slurm_header(user_input)
    slurm_script = Slurm(slurm_header)
    slurm_script_max = None
    slurm_script_nmax = None
    if n_arrays_max != 0:
        command = build_command(amber, user_input.processors_per_node)
        slurm_script_max = slurm_script.create_slurm_script(command, title, n_arrays_max)
    if n_trajectories_remaining != 0:
        command = build_command(amber, n_trajectories_remaining)
        slurm_script_nmax = slurm_script.create_slurm_script(command, title, 1)
    return slurm_script_max, slurm_script_nmax

def run_nasqm_slurm_files(slurm_files):
    '''
    Run the files produced by slurm_trajectory_files
    '''
    run_slurm(slurm_files[0])
    run_slurm(slurm_files[1])
