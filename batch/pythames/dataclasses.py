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
from typing import List

import dataclasses_json

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class Metadata:
    name: str
    description: str
    datetimeCreated: str
    datetimeModified: str
    isInputFile: bool

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class SpacecraftParameters:
    mass: float
    dragArea: float
    Cd: float

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class GeopotentialPerturbationParameters:
    isEnabled: bool
    model: str
    maxOrder: int
    maxDegree: int

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class AtmospherePerturbationParameters:
    isEnabled: bool
    model: str

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class PerturbationParameters:
    geopotential: GeopotentialPerturbationParameters
    atmosphere: AtmospherePerturbationParameters

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class PropagatorParameters:
    startTime: float
    endTime: float
    equations: str
    isNonDimensional: bool
    isFixedStep: bool
    intermediateOutput: bool
    timeStepIntermediate: float
    timeStep: float
    absoluteTolerance: float
    relativeTolerance: float

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class PolynomialParameters:
    isEnabled: bool
    type: str
    maxDegree: int

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class StateParameters:
    datetime: float
    states: List[List[float]]
    statetype: str

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class ExecutionStatistics:
    propagationTime: float

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class Parameters:
    metadata: Metadata
    spacecraft: SpacecraftParameters
    perturbation: PerturbationParameters
    propagator: PropagatorParameters
    polynomial: PolynomialParameters
    states: List[StateParameters]
    statistics: ExecutionStatistics

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class ParametersSet:
    parameters: List[Parameters]