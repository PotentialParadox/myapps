'''
A unit tester for slurm.py
'''
from slurm import Slurm

SLURM_HEADER = {'email': "dtracy.uf@gmail.com", 'email_options': 3,
                'n_nodes': 1, 'ppn': 2, 'memory': "3000mb", 'walltime': "00:01:00",
                'max_jobs': 8}

def test_email_options():
    '''
    test the email option
    '''
    slurm_object = Slurm(SLURM_HEADER)
    assert slurm_object.email_preferences == "BEGIN,END"

def test_create_slurm_script():
    '''
    Compare to a known working version of slurm
    '''
    command = 'echo "In this case SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID" >>' \
              ' result_${SLURM_ARRAY_TASK_ID}.out'
    slurm_object = Slurm(SLURM_HEADER)
    job_script = slurm_object.create_slurm_script(command, "MyJob", 3)
    # open('slurm_test.txt', 'w').write(job_script)
    test_script = str(open('slurm_test.txt', 'r').read())
    assert job_script == test_script