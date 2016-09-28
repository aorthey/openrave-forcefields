from TOPP import Utilities
import TOPP
import numpy as np
from topp_interface import *
from environment_force_blood_stream import *

Ndim = 4
discrtimestep= 0.001

env = EnvironmentBloodStream()
TOPPInterface.LoadFromFile('topp/clc', env)
