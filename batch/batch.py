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

import dataclasses
import datetime
import os

import numpy as np
import pandas as pd
import smt.sampling_methods

import pythames.dataclasses
import pythames.interface
import pythames.permutations
import pythames.statistics

def main():
    # Set filepaths
    BATCH_PARENT_DIR = os.path.dirname(os.path.realpath(__file__))
    COMMAND = os.path.join(BATCH_PARENT_DIR, "..", "bin", "thames_main")
    DATETIME = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    BATCH_OUTPUT_DIR = os.path.join(BATCH_PARENT_DIR, "output", DATETIME)
    CSV_OUTPUT_PATH = os.path.join(BATCH_OUTPUT_DIR, f"{DATETIME}_bulk.csv")

    ## Define batch parameters
    # Generate samples
    RV = np.array([6800.0, 0.0, 0.0, 0.0, 5.4137654, 5.4137654])
    RVunc = np.array([0.1, 0.1, 0.1, 0.01, 0.01, 0.01])
    T = 5580.515898

    # Generate initial state boundaries
    RVupper = RV + RVunc
    RVlower = RV - RVunc
    RVbound = np.vstack([RVlower, RVupper]).T

    # Sample initial boundaries
    lhs = smt.sampling_methods.LHS(xlimits=RVbound)
    RVsamples = lhs(int(1e4)).tolist()

    ## Generate parameter set
    # Generate propagator sets
    propagator = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.PropagatorParameters,
        pythames.permutations.PROPAGATORPARAMETERS_DEFAULT,
        endTime=[5*T],
        isFixedStep=[False],
        absoluteTolerance=[1e-13],
        relativeTolerance=[1e-13]
    )
    propagator_polynomial = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.PropagatorParameters,
        pythames.permutations.PROPAGATORPARAMETERS_DEFAULT,
        endTime=[5*T],
        isFixedStep=[True],
        timeStep=[T/400, T/200, T/150, T/125, T/100, T/80, T/60, T/50, T/40, T/30, T/20, T/10, T/5, T/2, T, 2.5*T, 5*T],
        equations=["Cowell", "GEqOE"]
    )

    # Generate perturbation sets
    geopotential = pythames.dataclasses.GeopotentialPerturbationParameters(True, "J2", 0, 0)
    atmosphere = pythames.dataclasses.AtmospherePerturbationParameters(True, "Wertz")
    perturbation = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.PerturbationParameters,
        pythames.permutations.PERTURBATIONPARAMETERS_DEFAULT,
        geopotential=[geopotential],
        atmosphere=[atmosphere]
    )

    # Generate state sets
    states = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.StateParameters,
        pythames.permutations.STATEPARAMETERS_DEFAULT,
        states=[RVsamples]
    )

    # Generate polynomial sets
    polynomial = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.PolynomialParameters,
        pythames.permutations.POLYNOMIALPARAMETERS_DEFAULT,
        isEnabled=[True],
        type=["Taylor"],
        maxDegree=[1, 2, 3, 4, 5, 6, 7, 8]
    )

    # Generate parameter sets
    parameters_point = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.Parameters,
        pythames.permutations.PARAMETERS_DEFAULT,
        propagator=propagator,
        perturbation=perturbation,
        states=[states]
    )
    parameters_polynomial = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.Parameters,
        pythames.permutations.PARAMETERS_DEFAULT,
        propagator=propagator_polynomial,
        perturbation=perturbation,
        states=[states],
        polynomial=polynomial
    )
    parameters = parameters_point + parameters_polynomial

    # Sort permutations
    parameters = sorted(parameters, key=lambda x: (x.propagator.timeStep, x.propagator.absoluteTolerance, -x.polynomial.maxDegree))

    ## Execute batch propagations 
    # Execute propagations
    param_prop = pythames.interface.batch_run(COMMAND, parameters, batchpath=BATCH_OUTPUT_DIR, parallel=True)

    ## Calculate statistics
    # Calculate bulk statistics
    param_prop_df = pythames.statistics.bulk_statistics(param_prop)

    ## Output statistics
    # Define output columns
    output_columns = [
        "propagator.timeStep",
        "propagator.equations",
        "polynomial.isEnabled",
        "polynomial.maxDegree",
        "datetime",
        "dr_rms","dv_rms", "dr_max", "dv_max",
        "rsw_r_rms","rsw_s_rms","rsw_w_rms","rsw_dom",
        "statistics.propagationTime"
    ]
    
    # Save statistics
    param_prop_df[output_columns].to_csv(CSV_OUTPUT_PATH, index=False)

if __name__ == "__main__":
    main()