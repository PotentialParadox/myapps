from multiprocessing import Process, Pool
import os
import subprocess

def run_amber(conjoined_list):
    cd_file = 'm1.inpcrd'
    input_name = conjoined_list[0]
    root_name = conjoined_list[1]
    coordinate_file = conjoined_list[2]
    amber = conjoined_list[3]
    if coordinate_file is not None:
        cd_file = coordinate_file
    subprocess.run([amber, '-O', '-i', input_name+'.in', '-o', root_name+'.out', '-c', cd_file, '-p',
                    'm1.prmtop', '-r', root_name+'.rst', '-x', root_name+'.nc'])


def run_amber_parallel(pmemd_available, input_names, root_names, coordinate_files, number_processors=4):
    if pmemd_available:
        amber = 'pmemd'
    else:
        amber = 'sander'
    p = Pool(number_processors)
    conjoined_list = []
    for i in range(len(root_names)):
        conjoined_list.append([input_names[i], root_names[i], coordinate_files[i], amber])
    p.map(run_amber, conjoined_list)
    p.close()
    p.join()
