#!/usr/bin/env python3

import pandas as pd

batchname = "sweep_2022-03-23_23:15:31.pkl"
batch = pd.read_pickle(batchname)
batch.to_csv(f"{batchname}.csv")