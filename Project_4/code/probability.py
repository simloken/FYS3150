from functions import monteCarlo, sorter, Stabilize, normer, probability
import multiprocessing as mp
import numpy as np
import itertools as it

import matplotlib.pyplot as plt

def proba(T, L, cycles):
    state=True; proba = True; accept=False
    stead = int(1e4)
    cycles = cycles
    store = monteCarlo(T,L, cycles, False, state, accept, proba)
    probability(store, L, T, stead)