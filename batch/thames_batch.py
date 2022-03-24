#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from io import StringIO
from smt.sampling_methods import LHS
from tqdm import tqdm
import time
from datetime import datetime
import multiprocessing

def thames_input(filepath, tstart, tend, points, scid, atol, rtol, statetype, degree):
    with open(filepath, "w") as f:
        f.write(f"{tstart}, {tend}, {scid}, {statetype}, {degree}, {atol}, {rtol}\n")
        for point in points:
            f.write(f"{point[0]:+8e}, {point[1]:+8e}, {point[2]:+8e}, {point[3]:+8e}, {point[4]:+8e}, {point[5]:+8e}\n")

def thames_output(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    states = " ".join(lines[:-1])
    states = np.genfromtxt(StringIO(states), delimiter=",")

    param = lines[-1].split(",")
    tstart = float(param[0])
    tend = float(param[1])
    scid = int(param[2])
    statetype = int(param[3])
    degree = int(param[4])
    atol = float(param[5])
    rtol = float(param[6])

    return tstart, tend, scid, states, statetype, degree, atol, rtol 

def thames_run(command, filepathin, filepathout, tstart, tend, states, degree, atol, rtol, statetype=0, scid=0):
    # Generate input files
    thames_input(filepathin, tstart, tend, states, scid, atol, rtol, statetype, degree)

    # Run THAMES
    start = time.perf_counter_ns()
    os.system(f"../bin/{command} {filepathin} {filepathout}")
    end = time.perf_counter_ns()
    delta = end - start

    # Load output
    _, _, _, states_propagated, _, _, _, _ = thames_output(filepathout)
    states_propagated = pd.DataFrame(states_propagated, columns=["X", "Y", "Z", "VX", "VY", "VZ"])
    states_propagated["proptype"] = command
    states_propagated["degree"] = degree
    states_propagated["atol"] = atol
    states_propagated["rtol"] = rtol
    states_propagated["time_total"] = delta*10**-9

    # Return output
    return states_propagated


filepathin = "./test_in.txt"
filepathout = "./test_out.txt"

tstart = 0.0
tend = 114000 #5*5400.0

# RV = np.array([-4102, 5093, -1866, -2.9, -4.3, -5.6])
# RVunc = np.array([0.1, 0.1, 0.1, 0.01, 0.01, 0.01])
# RVupper = RV + RVunc
# RVlower = RV - RVunc
RVupper = np.array([-4416782.92246406, -5263218.34932045, 2773.68818819729, 7656.46288470595, -6412.88169544488, 0.882777270214449])/1000.0
RVlower = np.array([-4422432.26516413, -5270764.01693497, -2693.42920650658, 7641.65569328887, -6424.59052100732, -0.819487850875413])/1000.0
RVbound = np.vstack([RVlower, RVupper]).T
sampling = LHS(xlimits=RVbound)
states = sampling(10**3)

command = ["thames_cowell_discrete_point"] #, "thames_geqoe_discrete_point"]
command_poly = ["thames_cowell_polynomial_point", "thames_geqoe_polynomial_point"]
degree = range(2, 11, 2)
tstep = np.array([15, 30, 45, 60, 90, 120, 180, 240, 360, 480, 1200, 1920, 2880, 3840, 5760, 7680])
# atol = 10.0**np.arange(-14, -9, 1)
atol = tstep

perms = [(commandi, degreei, atoli) for commandi in command for degreei in [0] for atoli in [10**(-13)]]
perms_poly = [(commandi, degreei, atoli) for commandi in command_poly for degreei in degree for atoli in atol]
perms.extend(perms_poly)
perms.reverse()
print(perms)

states_propagated = []

def worker(perm):
    commandi, degreei, atoli = perm
    name = multiprocessing.current_process().name
    statei = thames_run(commandi, f"{filepathin}{name}", f"{filepathout}{name}", tstart, tend, states, degreei, atoli, atoli)
    return statei

# for perm in tqdm(perms):
#     statei = worker(perm)
#     states_propagated.append(statei)

with multiprocessing.Pool() as pool:
    with tqdm(total=len(perms)) as pbar:
        for statei in pool.imap_unordered(worker, perms):
            states_propagated.append(statei)
            pbar.update(1)

now = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
pkl_path = f"sweep_{now}.pkl"

states_propagated = pd.concat(states_propagated)
states_propagated.to_pickle(pkl_path)
print(states_propagated.head())
print(states_propagated.tail())