#!/usr/bin/env python3

import copy
import dataclasses
import datetime
import os
import statistics
from typing import List
import itertools
import multiprocessing

import dataclasses_json
import numpy as np
import smt.sampling_methods
import tqdm
import pandas as pd

BATCH_PATH = os.path.dirname(os.path.realpath(__file__))

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
            [states],
            ExecutionStatistics(0.0)
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

# Define worker function
def worker(parameter: Parameters) -> Parameters:
    # Find process name
    name = multiprocessing.current_process().name

    # Generate I/O filepaths
    command = os.path.join(BATCH_PATH, "..", "bin", "thames_main")
    filepathin = os.path.join(BATCH_PATH, f"batch_{name}_in.json")
    filepathout =  os.path.join(BATCH_PATH, f"batch_{name}_out.json")

    # Execute THAMES
    parameters_propagated = run(command, parameter, filepathin, filepathout)

    # Return propagated parameters
    return parameters_propagated

def main():
    # Set parameters
    rms_limits = [(0, 3), (3, 6)]

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
            "equations": ["Cowell"],
            "isNonDimensional" : [True],
            "isFixedStep": [False],
            "intermediateOutput": [False],
            "timeStepIntermediate": [30],
            "timeStep": [30],
            "absoluteTolerance": [1e-14],
            "relativeTolerance": [1e-14]
        },
        "polynomial": {
            "isEnabled": [False],
            "type": [""],
            "maxDegree": [0]
        },
    }

    # Generate parameter sets for polynomial propagations
    parameter_lists_polynomial = copy.deepcopy(parameter_lists)
    # Set propagator parameters
    parameter_lists_polynomial["propagator"]["equations"] = ["GEqOE"]
    parameter_lists_polynomial["propagator"]["isNonDimensional"] = [True, False] 
    parameter_lists_polynomial["propagator"]["isFixedStep"] = [False]
    parameter_lists_polynomial["propagator"]["absoluteTolerance"] = [1e-13]
    parameter_lists_polynomial["propagator"]["relativeTolerance"] = [1e-13]
    # Set polynomial parameters
    parameter_lists_polynomial["polynomial"]["isEnabled"] = [True]
    parameter_lists_polynomial["polynomial"]["type"] = ["Taylor"]
    parameter_lists_polynomial["polynomial"]["maxDegree"] = [6]

    # Set constant parameters (i.e. perturbation model, and states)
    perturbation = PerturbationParameters(
        GeopotentialPerturbationParameters(True, "J2", 0, 0),
        AtmospherePerturbationParameters(False, "Wertz")
    )
    states = StateParameters(
        0,
        states,
        "Cartesian"
    )

    # Generate parameter permutations
    parameters = [
        *parameter_permutations(parameter_lists, perturbation, states),
        *parameter_permutations(parameter_lists_polynomial, perturbation, states)
    ]

    # Remove permutations where the absolute tolerance does not match the relative tolerance
    parameters = [param for param in parameters if param.propagator.absoluteTolerance == param.propagator.relativeTolerance]

    # Sort permutations
    parameters = sorted(parameters, key=lambda x: (x.propagator.timeStep, x.propagator.absoluteTolerance, -x.polynomial.maxDegree))

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

    # Generate tables of results
    parameters_propagated_tables = list()
    for ip in parameters_propagated:
        # Declare table for each propagation
        ip_table = pd.DataFrame()

        # Iterate through states
        for state in ip.states:
            # Declare table for each set of states
            state_table = pd.DataFrame()

            # Generate labels
            state_labels = [f"state_{stri}" for stri in range(0, len(state.states[0]))]

            # Add states, datetime, and statetype to the table
            state_table[state_labels] = np.array(state.states)
            state_table["datetime"] = state.datetime
            state_table["statetype"] = state.statetype

            # Append data to overall propagation table
            ip_table = pd.concat((ip_table, state_table))

        # Add metadata
        propagator = ip.propagator.to_dict()
        polynomial = ip.polynomial.to_dict()
        statistics = ip.statistics.to_dict()
        ip_table[list(propagator.keys())] = list(propagator.values())
        ip_table[list(polynomial.keys())] = list(polynomial.values())
        ip_table[list(statistics.keys())] = list(statistics.values())

        # Store table
        parameters_propagated_tables.append(ip_table)

    # Pick reference solution
    ref_solution = sorted([iparam for iparam in parameters_propagated_tables if not iparam["isEnabled"].all()], key=lambda x: (x["absoluteTolerance"].iloc[0], x["timeStep"].iloc[0]))[0]

    # Set pivot parameters
    pivot_parameters = []
    pivot_parameters.append("datetime")
    pivot_parameters.extend([param for param in list(PropagatorParameters.__dataclass_fields__)])
    pivot_parameters.extend([param for param in list(PolynomialParameters.__dataclass_fields__)])
    pivot_parameters.extend([param for param in list(ExecutionStatistics.__dataclass_fields__)])

    # Calculate bulk statistics
    bulk_table = pd.DataFrame()
    for table in tqdm.tqdm(parameters_propagated_tables):
        # Calculate deltas in state variables
        state_cols = [col for col in table.columns if "state_" in col]
        dstate_cols = [f"delta_{col}" for col in state_cols]
        table[dstate_cols] = table[state_cols] - ref_solution[state_cols]

        # Iterate through slices
        for slice in rms_limits:
            # Define related columns
            state_cols = [f"delta_state_{col}" for col in range(*slice)]
            dstate_cols = f"delta_state_{slice[0]}_{slice[-1]-1}"

            # Iterate through times
            for t in np.unique(table["datetime"]):
                # Calculate slice magnitude
                table.loc[(table["datetime"] == t), dstate_cols] = np.linalg.norm(table.loc[(table["datetime"] == t), state_cols], axis=1)

        # Condense, and append to bulk table
        bulk_table = pd.concat((
            bulk_table,
            pd.pivot_table(
                table,
                values=[f"delta_state_{slice[0]}_{slice[-1]-1}" for slice in rms_limits],
                index=pivot_parameters,
                aggfunc=lambda x: np.sqrt(np.mean(x**2))
                )
            )
        )

    # Calculate current time for output file name
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    csv_path = os.path.join(BATCH_PATH, f"batch_{now}.csv") 

    # Save bulk table
    bulk_table.to_csv(csv_path)

if __name__ == "__main__":
    main()