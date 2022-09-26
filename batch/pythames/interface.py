#!/usr/bin/env python3

# MIT License
#
# Copyright (c) 2021-2022 Max Hallgarten La Casta
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import datetime
import multiprocessing
import os
from typing import List, Tuple
import uuid

import tqdm

from .dataclasses import Parameters

def save(filepath: str, parameters: Parameters) -> None:
    # Update times
    now = datetime.datetime.utcnow().isoformat(sep='T', timespec='milliseconds') + "Z"
    parameters.metadata.datetimeCreated = now
    parameters.metadata.datetimeModified = now

    # Open file
    with open(filepath, "w") as fid:
        # Output parameters in JSON format
        fid.write(parameters.to_json(indent=4))

def load(filepath: str) -> Parameters:
    # Read file
    with open(filepath, "r") as fid:
        raw = fid.read()

    # Create Parameters
    parameters = Parameters.from_json(raw)

    # Return parameters
    return parameters

def run(command: str, parametersin: Parameters, filepathin: str, filepathout: str) -> Parameters:
    # Save input file
    save(filepathin, parametersin)

    # Run command
    os.system(f"{command} {filepathin} {filepathout}")

    # Load output file
    parametersout = load(filepathout)

    # Return parameters
    return parametersout

def worker_run(input: Tuple[str, Parameters, str, str]) -> Parameters:
    # Run command
    return run(*input)

def batch_run(command: str, parametersin: List[Parameters], batchpath=None, parallel=False) -> List[Parameters]:
    # Declare output list
    parametersout = []

    # Define (and create) output directory
    folderpath = os.getcwd()
    if batchpath is None: batchpath = os.path.join(folderpath, "output", datetime.datetime.utcnow().isoformat(sep='T', timespec="seconds"))
    if not os.path.exists(batchpath): os.makedirs(batchpath)

    # Length
    nlen = len(parametersin)

    # Generate list of IDs
    ids = [str(uuid.uuid4().hex) for _ in range(nlen)]

    # Generate filepaths
    filepaths = [os.path.join(batchpath, f"{id}.json") for id in ids]

    # Generate inputs
    inputs = zip([command]*nlen, parametersin, filepaths, filepaths)

    # Iterate through input parameters
    if parallel:
        with multiprocessing.Pool() as pool:
             # Create progress bar
            with tqdm.tqdm(total=nlen) as pbar:
                # Iterate through permutations
                for iparamout in pool.imap_unordered(worker_run, inputs):
                    # Store propagated state
                    parametersout.append(iparamout)

                    # Update progress bar
                    pbar.update(1)
    else:
        for input in tqdm.tqdm(inputs, total=nlen):
            # Run propagation
            iparamout = run(*input)

            # Append to output list
            parametersout.append(iparamout)

    # Return output
    return parametersout