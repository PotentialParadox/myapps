'''
This the User Input file for the NASQM
Automation Script
'''
class UserInput:
    '''
    Contains the data for the users to change at will
    '''
    def __init__(self):
        # Change here whether your simulation is qmmm
        self.is_qmmm = False
        # Change here whether you are working on your personal computer
        # or an HPC with SLUM
        self.is_hpc = False
        # Are you performing tully surface hopping
        self.is_tully = False
        # How many processors will be on a node?
        self.processors_per_node = 8
        # Do you want to run ground state dynamics
        self.run_ground_state_dynamics = True
        # Do you want to run the trajectories used for the abjorption specta
        self.run_absorption_trajectories = False
        # Do you want to run the single point snapshots from these
        # absorption spectra trajectories?
        self.run_absorption_snapshots = False
        # Do you want to collect the data from the absorption calculations?
        self.run_absorption_collection = False
        # Do you want to run the exctied state trajectories?
        self.run_excited_state_trajectories = True
        # Do you want to collect the data from the exctied state trajectory
        # calculations?
        self.run_fluorescence_collection = False

        # Change here the number of snapshots you wish to take
        # from the initial ground state trajectory to run the
        # further ground state dynamics
        self.n_snapshots_gs = 1000

        # Change here the number of states you wish to
        # calculate in the absorption singlpoint calculations
        self.n_abs_exc = 1

        # Change here the number of snapshots you wish to take
        # from the initial ground state trajectory to run the
        # new excited state dynamics
        self.n_snapshots_ex = 8

        # Change here the time step that will be shared by
        # each trajectory
        self.time_step = 0.5  # fs

        # Change here the runtime of the initial ground state MD
        self.ground_state_run_time = 1 # ps

        # Change here the runtime for the the trajectories
        # used to create calculated the absorption
        self.abs_run_time = 1 # ps

        # Change here the runtime for the the trajectories
        # used to create calculated the fluorescence
        self.exc_run_time = 1 # ps

        # Change here the number of excited states you
        # with to have in the CIS calculation
        self.n_exc_states_propagate_ex_param = 1

        # Change here the initial state
        self.exc_state_init_ex_param = 1

        # Some time will be needed for the molecule to equilibrate
        # from jumping from the ground state to the excited state.
        # We don't want to include this data in the calculation
        # of the fluorescence. We therefore set a time delay.
        self.fluorescene_time_delay = 1 # PS


