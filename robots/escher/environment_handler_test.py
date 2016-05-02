
import numpy as np
import openravepy as rave
import time
import random
import math

class environment_handler:
    def __init__(self,mode='box_world'):
        self.env = rave.Environment()  # create openrave environment
        self.env.SetViewer('qtcoin')  # attach viewer (optional)
        self.env.GetViewer().SetCamera([
            [0.31499128, 0.09759726, -0.94406317, 6.81987572],
            [0.94805905, 0.01409698, 0.31778187, -2.29564428],
            [0.04432308, -0.99512615, -0.08808754, 1.60788679],
            [0., 0., 0., 1.]
        ])

        # predefined backgound structures
        self.env.Load('ground_plane_2.xml')

        self.surface_list = None
        self.surface_mesh_list = None
        self.structures = None
        self.mode = mode
