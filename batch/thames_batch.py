#!/usr/bin/env python3

# Built-in packages
from datetime import datetime
from io import StringIO
import multiprocessing
import os
import random
import time

# Third-party packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from smt.sampling_methods import LHS
from tqdm import tqdm

def thames_input(filepath, tstart, tend, points, scid, atol, rtol, statetype, degree):
    # Open filepath for input file
    with open(filepath, "w") as f:
        # Write header
        f.write(f"{tstart}, {tend}, {scid}, {statetype}, {degree}, {atol}, {rtol}\n")

        # Iterate through points, and write to input file
        for point in points:
            f.write(f"{point[0]:+8e}, {point[1]:+8e}, {point[2]:+8e}, {point[3]:+8e}, {point[4]:+8e}, {point[5]:+8e}\n")

def thames_output(filepath):
    # Open filepath for output file
    with open(filepath, "r") as f:
        # Read file
        lines = f.readlines()

    # Read output states, and convert to a Numpy array
    states = " ".join(lines[:-1])
    states = np.genfromtxt(StringIO(states), delimiter=",")

    # Split footer by comma
    param = lines[-1].split(",")

    # Read footer parameters
    tstart = float(param[0])
    tend = float(param[1])
    scid = int(param[2])
    statetype = int(param[3])
    degree = int(param[4])
    atol = float(param[5])
    rtol = float(param[6])

    # Return parameters and states
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

# Set I/O filepaths
filepathin = "./batch_in"
filepathout = "./batch_out"

# Set propagation start and end time
tstart = 0.0
tend = 5400.0

# Set initial states and uncertainties
RV = np.array([-4102, 5093, -1866, -2.9, -4.3, -5.6])
RVunc = np.array([0.1, 0.1, 0.1, 0.01, 0.01, 0.01])

# Generate initial state boundaries
RVupper = RV + RVunc
RVlower = RV - RVunc
RVbound = np.vstack([RVlower, RVupper]).T

# Sample initial boundaries
sampling = LHS(xlimits=RVbound)
states = sampling(10**3)

# Set desired point, and polynomial commands
command = ["thames_cowell_discrete_point", "thames_geqoe_discrete_point"]
command_poly = ["thames_cowell_polynomial_point", "thames_geqoe_polynomial_point"]

# Set desired polynomial degrees
degree = range(2, 11, 2)

# Set desired solver tolerances
atol = 10.0**np.arange(-14, -9, 1)

# Generate permutations for point, and polynomial commands
perms = [(commandi, degreei, atoli) for commandi in command for degreei in [0] for atoli in atol]
perms_poly = [(commandi, degreei, atoli) for commandi in command_poly for degreei in degree for atoli in atol]
perms.extend(perms_poly)

# Sort to start with most computationally expensive permutations
perms.sort(key= lambda x: (-x[1], x[2]))

# Declare empty list to store propagated states
states_propagated = []

def worker(perm):
    # Unpack parameters
    commandi, degreei, atoli = perm

    # Retrieve worker name
    name = multiprocessing.current_process().name

    # Execute THAMES, appending the worker name to the filepaths to avoid interaction between the workers
    statei = thames_run(commandi, f"{filepathin}_{name}.txt", f"{filepathout}_{name}.txt", tstart, tend, states, degreei, atoli, atoli)

    # Return propagated state
    return statei

# Create worker pool
with multiprocessing.Pool() as pool:
    # Create progress bar
    with tqdm(total=len(perms)) as pbar:
        # Iterate through permutations
        for statei in pool.imap_unordered(worker, perms):
            # Store propagated state
            states_propagated.append(statei)

            # Update progress bar
            pbar.update(1)

# Calculate current time for output file name
now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
csv_path = f"batch_{now}.csv"

# Generate pandas table of propagated states
states_propagated = pd.concat(states_propagated)

# Save propagated states
states_propagated.to_csv(csv_path)
