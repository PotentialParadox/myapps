'''
Interface to Amber
'''
from multiprocessing import Pool
import subprocess


class Amber:
    '''
    A class to control calls to amber
    '''
    def __init__(self, input_files=None, output_files=None, coordinate_files=None,
                 prmtop_files=None, restart_files=None, export_files=None):
        self.input_files = input_files
        self.output_files = output_files
        self.coordinate_files = coordinate_files
        self.prmtop_files = prmtop_files
        self.restart_files = restart_files
        self.export_files = export_files

    def run_amber(self, conjoined_list):
        '''
        Single thread amber run
        '''
        input_file = conjoined_list[0]
        output_file = conjoined_list[1]
        coordinate_file = conjoined_list[2]
        prmtop_file = conjoined_list[3]
        restart_file = conjoined_list[4]
        export_file = conjoined_list[5]
        subprocess.run(['sander', '-O', '-i', input_file, '-o', output_file, '-c',
                        coordinate_file, '-p', prmtop_file, '-r', restart_file, '-x', export_file])

    def run_amber_parallel(self, number_processors=4):
        '''
        Multithreaded amber run
        '''
        pool = Pool(number_processors)
        conjoined_list = []
        for i in range(len(self.input_files)):
            conjoined_list.append([self.input_files[i], self.output_files[i],
                                   self.coordinate_files[i], self.prmtop_files[i],
                                   self.restart_files[i], self.export_files[i]])
        pool.map(self.run_amber, conjoined_list)
        pool.close()
        pool.join()
