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

import numpy as np
import pandas as pd

from .dataclasses import Parameters

def bulk_statistics(parameters: List[Parameters]) -> pd.DataFrame:
    # Condense parameters in dataframe
    # TODO: clean up
    param_df = pd.concat([pd.concat([pd.json_normalize(dataclasses.asdict(ip)).explode("states").drop(["states"], axis=1), pd.json_normalize(dataclasses.asdict(ip)).explode("states")["states"].apply(pd.Series)], axis=1).reset_index(drop=True) for ip in parameters])

    # Expand states, and ensure they are numpy arrays
    param_df["states"] = param_df["states"].apply(lambda x: np.array(x, dtype=np.float))

    # Extract position and velocity vectors
    param_df["r"] = param_df["states"].apply(lambda x: x[:, 0:3])
    param_df["v"] = param_df["states"].apply(lambda x: x[:, 3:6])

    # Select reference solution
    ref_sol = pd.DataFrame(param_df[~param_df["polynomial.isEnabled"]])
    # Calculate RSW frame
    ref_sol["rh"] = ref_sol.apply(lambda x: x.r/np.linalg.norm(x.r, axis=1).reshape(-1,1), axis=1)
    ref_sol["wh"] = ref_sol.apply(lambda x: np.cross(x.r, x.v, axis=1)/np.linalg.norm(np.cross(x.r, x.v, axis=1), axis=1).reshape(-1,1), axis=1)
    ref_sol["sh"] = ref_sol.apply(lambda x: np.cross(x.wh, x.r, axis=1)/np.linalg.norm(np.cross(x.wh, x.r, axis=1), axis=1).reshape(-1,1), axis=1)

    # Join reference solution
    param_df = param_df.join(ref_sol[["r", "v", "rh", "wh", "sh"]], rsuffix="_ref")

    # Calculate solution errors
    param_df["dr"] = param_df["r"] - param_df["r_ref"]
    param_df["dv"] = param_df["v"] - param_df["v_ref"]
    param_df["rsw"] = param_df.apply(lambda x: np.column_stack((np.sum(x.dr*x.rh, axis=1), np.sum(x.dr*x.sh, axis=1), np.sum(x.dr*x.wh, axis=1))), axis=1)

    # Calculate RMSEs
    param_df["dr_rms"] = param_df["dr"].apply(lambda x: np.sqrt(np.mean(np.linalg.norm(x, axis=1)**2)))
    param_df["dv_rms"] = param_df["dv"].apply(lambda x: np.sqrt(np.mean(np.linalg.norm(x, axis=1)**2)))
    param_df["rsw_r_rms"] = param_df["rsw"].apply(lambda x: np.sqrt(np.mean(x[:,0]**2)))
    param_df["rsw_s_rms"] = param_df["rsw"].apply(lambda x: np.sqrt(np.mean(x[:,1]**2)))
    param_df["rsw_w_rms"] = param_df["rsw"].apply(lambda x: np.sqrt(np.mean(x[:,2]**2)))

    # Calculate mean RSW vector, and dominant error direction
    param_df["rsw_mean"] = param_df["rsw"].apply(lambda x: np.mean(x, axis=0))
    param_df["rsw_dom"] = param_df["rsw_mean"].apply(lambda x: np.argmax(np.abs(x)))

    # Return bulk statistics
    return param_df