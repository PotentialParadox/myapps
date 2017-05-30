'''
Allow the NASQM automation script to interface with slurm
'''
def create_slurm_header(user_input):
    '''
    Returns the slurm header dictionary
    '''
    return {'email': user_input.email, 'email_options': user_input.email_options,
            'n_nodes': user_input.number_nodes, 'ppn': user_input.processors_per_node,
            'memory': user_input.memory_per_node, 'walltime': user_input.walltime,
            'max_jobs': user_input.max_jobs}
