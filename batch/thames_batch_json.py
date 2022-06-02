#!/usr/bin/env python3

import dataclasses
import datetime
import os
from typing import List
import itertools
import multiprocessing

import dataclasses_json
import numpy as np
import smt.sampling_methods
import tqdm

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
class Parameters:
    metadata: Metadata
    spacecraft: SpacecraftParameters
    perturbation: PerturbationParameters
    propagator: PropagatorParameters
    polynomial: PolynomialParameters
    states: List[StateParameters]

@dataclasses_json.dataclass_json
@dataclasses.dataclass
class ParametersSet:
    parameters: List[Parameters]

def parameter_permutations(parameter_lists: dict, perturbation: PerturbationParameters, states: StateParameters) -> List[Parameters]:
    # Declare empty dictionary for parameters
    parameters = dict()

    # Iterate through parameter fields
    for key, value in parameter_lists.items():
        # Generate permutations
        field_permutations = itertools.product(*(value[ikey] for ikey in value))
        
        # Append to overall dictionary
        parameters[key] = list(field_permutations)

    # Generate permutations
    permutations = itertools.product(*(parameters[ikey] for ikey in parameters))

    # Construct parameter objects
    parameters = [
        Parameters(
            Metadata(*metadata),
            SpacecraftParameters(*spacecraft),
            perturbation,
            PropagatorParameters(*propagator), 
            PolynomialParameters(*polynomial),
            [states]
            )
        for metadata, spacecraft, propagator, polynomial in permutations
    ]

    # Prune parameters where atol != rtol
    parameters = [
        parameter
        for parameter in parameters
        if parameter.propagator.absoluteTolerance == parameter.propagator.relativeTolerance
    ]

    # Return list of parameters
    return parameters

def save(filepath: str, parameters: Parameters) -> None:
    # Update times
    now = datetime.datetime.utcnow().isoformat(sep='T', timespec='milliseconds') + "Z"
    parameters.metadata.datetimeCreated = now
    parameters.metadata.datetimeModified = now

    # Open file
    with open(filepath, "w") as fid:
        # Output parameters in JSON format
        fid.write(parameters.to_json())

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

# Define worker function
def worker(parameter: Parameters) -> Parameters:
    # Find process name
    name = multiprocessing.current_process().name

    # Define I/O file directory, and executable path
    filedir = "batch"
    command = "../bin/thames_main"

    # Generate I/O filepaths
    filepathin = f"{filedir}_{name}_in.json"
    filepathout = f"{filedir}_{name}_out.json"

    # Execute THAMES
    parameters_propagated = run(command, parameter, filepathin, filepathout)

    # Return propagated parameters
    return parameters_propagated

def main():
    # Set initial states and uncertainties
    RV = np.array([-4102, 5093, -1866, -2.9, -4.3, -5.6])
    RVunc = np.array([0.1, 0.1, 0.1, 0.01, 0.01, 0.01])

    # Generate initial state boundaries
    RVupper = RV + RVunc
    RVlower = RV - RVunc
    RVbound = np.vstack([RVlower, RVupper]).T

    # Sample initial boundaries
    sampling = smt.sampling_methods.LHS(xlimits=RVbound)
    states = sampling(10**3)
    states = states.tolist()

    # Set lists for generating permutations
    parameter_lists = {
        "metadata" : {
            "name" : ["Batch run"],
            "description": ["Automatically generated with Python"],
            "datetimeCreated": [""],
            "datetimeModified": [""],
            "isInputFile": [True]
        },
        "spacecraft": {
            "mass": [200],
            "dragArea": [10e-6],
            "Cd": [2.2]
        },
        "propagator": {
            "startTime": [0],
            "endTime": [5400],
            "equations": ["Cowell", "GEqOE"],
            "isNonDimensional" : [True],
            "isFixedStep": [False],
            "intermediateOutput": [False],
            "timeStep": [30],
            "absoluteTolerance": [1e-14],
            "relativeTolerance": [1e-14]
        },
        "polynomial": {
            "isEnabled": [False],
            "type": ["Taylor"],
            "maxDegree": [2, 4, 6, 8, 10]
        },
    }

    # Set constant parameters (i.e. perturbation model, and states)
    perturbation = PerturbationParameters(
        GeopotentialPerturbationParameters(True, "J2", 0, 0),
        AtmospherePerturbationParameters(True, "USSA76")
    )
    states = StateParameters(
        0,
        states,
        "Cartesian"
    )

    # Generate parameter permutations
    parameters = parameter_permutations(parameter_lists, perturbation, states)

    # Iterate through permutations
    parameters_propagated = list()
    # Create worker pool
    with multiprocessing.Pool() as pool:
        # Create progress bar
        with tqdm.tqdm(total=len(parameters)) as pbar:
            # Iterate through permutations
            for parami in pool.imap_unordered(worker, parameters):
                # Store propagated state
                parameters_propagated.append(parami)

                # Update progress bar
                pbar.update(1)
    
    # Create set of parameters
    parameters_propagated_set = ParametersSet(parameters_propagated)

    # Calculate current time for output file name
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filepathbatch = f"batch_{now}.json"

    # Save parameter sets
    with open(filepathbatch, "w") as fid:
        # Output parameters in JSON format
        fid.write(parameters_propagated_set.to_json())

if __name__ == "__main__":
    main()