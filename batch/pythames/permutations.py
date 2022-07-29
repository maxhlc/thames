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

import itertools
from typing import List

from .dataclasses import *

def dataclass_permutations(cls: type, default: dict, **kwargs) -> List[object]:
    # Copy default dictionary
    default_copy = default.copy()

    # Update using keyword arguments
    for key, value in kwargs.items():
        default_copy[key] = value

    # Generate permutations
    keys, values = zip(*default_copy.items())
    permutations = [dict(zip(keys, v)) for v in itertools.product(*values)]

    # Generate metadata objects
    cls_permutations = [cls(**iperm) for iperm in permutations]

    # Return permutations
    return cls_permutations

METADATA_DEFAULT = {
    "name" : ["Batch run"],
    "description": ["Automatically generated with Python"],
    "datetimeCreated": [""],
    "datetimeModified": [""],
    "isInputFile": [True]
}

SPACECRAFTPARAMETERS_DEFAULT = {
    "mass": [0],
    "dragArea": [0],
    "Cd": [0]
}

GEOPOTENTIALPERTURBATIONPARAMETERS_DEFAULT = {
    "isEnabled": [False],
    "model": [""],
    "maxOrder": [0],
    "maxDegree": [0]
}

ATMOSPHEREPERTURBATIONPARAMETERS_DEFAULT = {
    "isEnabled": [False],
    "model": [""]
}

PERTURBATIONPARAMETERS_DEFAULT = {
    "geopotential": dataclass_permutations(GeopotentialPerturbationParameters, GEOPOTENTIALPERTURBATIONPARAMETERS_DEFAULT),
    "atmosphere": dataclass_permutations(AtmospherePerturbationParameters, ATMOSPHEREPERTURBATIONPARAMETERS_DEFAULT)
}

PROPAGATORPARAMETERS_DEFAULT = {
    "startTime": [0],
    "endTime": [86400.0],
    "equations": ["Cowell"],
    "isNonDimensional": [True],
    "isFixedStep": [True],
    "intermediateOutput": [False],
    "timeStepIntermediate": [30],
    "timeStep": [30],
    "absoluteTolerance": [1e-14],
    "relativeTolerance": [1e-14]
}

POLYNOMIALPARAMETERS_DEFAULT = {
    "isEnabled": [False],
    "type": [""],
    "maxDegree": [0]
}

STATEPARAMETERS_DEFAULT = {
    "datetime": [0.0],
    "states": [
        [[7000.0, 0.0, 0.0, 0.0, 8.0, 0.0]]
    ],
    "statetype": ["Cartesian"]
}

EXECUTIONSTATISTICS_DEFAULT = {
    "propagationTime": [0.0]
}

PARAMETERS_DEFAULT = {
    "metadata": dataclass_permutations(Metadata, METADATA_DEFAULT),
    "spacecraft": dataclass_permutations(SpacecraftParameters, SPACECRAFTPARAMETERS_DEFAULT),
    "perturbation": dataclass_permutations(PerturbationParameters, PERTURBATIONPARAMETERS_DEFAULT),
    "propagator": dataclass_permutations(PropagatorParameters, PROPAGATORPARAMETERS_DEFAULT),
    "polynomial": dataclass_permutations(PolynomialParameters, POLYNOMIALPARAMETERS_DEFAULT),
    "states": [dataclass_permutations(StateParameters, STATEPARAMETERS_DEFAULT)],
    "statistics": dataclass_permutations(ExecutionStatistics, EXECUTIONSTATISTICS_DEFAULT)
}