'''
A class to contain the data use to compare ouputs from nexmd
'''
import io
import pynasqm.amberout as amber_out
import numpy as np

class NexmdOutput:
    def __init__(self, file_name):
        self.file_name = file_name
        self.n_excited_states = None
        self.gs_energies = None
        self.es_energies = None
        self.omegas = None
        self.orbitals = None
        self.read_file()

    def find_gs_energies(self):
        file_stream = open(self.file_name, 'r')
        gs_energies = amber_out.find_ground_energies(file_stream, None)
        self.gs_energies = np.fromstring(gs_energies, dtype=float, sep="\n")
        return self.gs_energies

    def find_es_energies(self):
        file_stream = open(self.file_name, 'r')
        es_range = range(1, self.n_excited_states+1)
        excited_states = amber_out.find_excited_energies(file_stream, None, es_range)
        self.es_energies = np.fromstring(excited_states, dtype=float, sep="\n")

    def find_omegas(self):
        file_stream = open(self.file_name, 'r')
        omegas = amber_out.find_nasqm_excited_state(file_stream, None,
                                                    range(1, self.n_excited_states+1))
        self.omegas = np.fromstring(omegas, dtype=float, sep=" ")
        return self.omegas

    def find_orbitals(self):
        file_stream = open(self.file_name, 'r')
        orbitals = amber_out.find_molecular_orbitals(file_stream, None)
        self.orbitals = np.fromstring(orbitals, dtype=float, sep=" ")
        return self.orbitals

    def read_file(self):
        self.n_excited_states = amber_out.find_number_excited_states(open(self.file_name, 'r'))
        self.find_gs_energies()
        self.find_es_energies()
        self.find_omegas()
        self.find_orbitals()
