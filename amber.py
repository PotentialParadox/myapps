from multiprocessing import Process, Pool
import os
import subprocess


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





