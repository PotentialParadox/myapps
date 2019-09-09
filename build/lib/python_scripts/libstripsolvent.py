from amber import run_amber_parallel
import subprocess

def stripper(n_ground_snaps, n_atoms, ppn):
    last_line = int(n_atoms / 2)
    snap_restarts=[]
    snap_trajectories=[]
    pmemd_available = False
    for i in range(n_ground_snaps):
        file_in = "ground_snap." + str(i+1)
        file_out = "ground_stripped." + str(i+1) + ".rst"
        file_string = ""
        cpptraj_string = "parm m1.prmtop\n"\
                       "trajin " + file_in + "\n"\
                       "strip :CL3\n"\
                       "trajout " + file_out + " restart\n"\
                       "run\n"\
                       "quit"
        open("stripper.traj", 'w').write(cpptraj_string)
        subprocess.run(['cpptraj','-i','stripper.traj'])
        snap_restarts.append("ground_stripped."+str(i+1)+".rst")
        snap_trajectories.append("nasqm_abs_" + str(i + 1))
    parm_stripper = "parm m1.prmtop\n"\
                    "parmstrip :CL3\n"\
                    "parmwrite out stripped.prmtop\n"\
                    "run\n"\
                    "quit\n"
    open("parm_stripper.traj", 'w').write(parm_stripper)
    subprocess.run(['cpptraj','-i','parm_stripper.traj'])
    subprocess.run(['cp','m1.prmtop','solvent.prmtop'])
    subprocess.run(['cp','stripped.prmtop','m1.prmtop'])
    run_amber_parallel(pmemd_available, snap_trajectories, snap_restarts, number_processors=ppn)
    subprocess.run(['cp','solvent.prmtop','m1.prmtop'])

