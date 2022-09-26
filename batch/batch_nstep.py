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
    numrange = np.logspace(2, 6, 11).tolist()
    RVsamples = [lhs(int(num)).tolist() for num in numrange]

    ## Generate parameter set
    # Generate propagator sets
    propagator = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.PropagatorParameters,
        pythames.permutations.PROPAGATORPARAMETERS_DEFAULT,
        endTime=[5*T],
        isFixedStep=[True],
        timeStep=[30],
        equations=["Cowell"]
    )
    propagator_polynomial = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.PropagatorParameters,
        pythames.permutations.PROPAGATORPARAMETERS_DEFAULT,
        endTime=[5*T],
        isFixedStep=[True],
        timeStep=[15],
        equations=["Cowell"]
    )

    # Generate perturbation sets
    geopotential = pythames.dataclasses.GeopotentialPerturbationParameters(True, "J2", 0, 0)
    atmosphere = pythames.dataclasses.AtmospherePerturbationParameters(True, "Wertz-P5")
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
        states=RVsamples
    )

    # Generate polynomial sets
    polynomial = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.PolynomialParameters,
        pythames.permutations.POLYNOMIALPARAMETERS_DEFAULT,
        isEnabled=[True],
        type=["Taylor"],
        maxDegree=[6]
    )

    # Generate parameter sets
    parameters_point = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.Parameters,
        pythames.permutations.PARAMETERS_DEFAULT,
        propagator=propagator,
        perturbation=perturbation,
        states=[[state] for state in states]
    )
    parameters_polynomial = pythames.permutations.dataclass_permutations(
        pythames.dataclasses.Parameters,
        pythames.permutations.PARAMETERS_DEFAULT,
        propagator=propagator_polynomial,
        perturbation=perturbation,
        states=[[state] for state in states],
        polynomial=polynomial
    )
    parameters = parameters_point + parameters_polynomial

    # Sort permutations
    # parameters = sorted(parameters, key=lambda x: (x.propagator.timeStep, x.propagator.absoluteTolerance, -x.polynomial.maxDegree))
    parameters = sorted(parameters, key=lambda x: -len(x.states[0].states))

    ## Execute batch propagations 
    # Execute propagations
    param_prop = pythames.interface.batch_run(COMMAND, parameters, batchpath=BATCH_OUTPUT_DIR, parallel=True)

    ## Output statistics
    isPoly = [param.polynomial.isEnabled for param in param_prop]
    nSamp = [len(param.states[0].states) for param in param_prop]
    compTime = [param.statistics.propagationTime for param in param_prop]

    df = pd.DataFrame(
        {
            "isPolynomial": isPoly,
            "n_samples": nSamp,
            "tcpu": compTime
        }
    )
    df.to_csv(CSV_OUTPUT_PATH, index=False)

if __name__ == "__main__":
    main()