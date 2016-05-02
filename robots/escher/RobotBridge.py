from __future__ import print_function
from __future__ import print_function
import openravepy as rave
import numpy as np
import math
import copy
import time

class AttributePassthrough(object):
    def __init__(self, getter, getAll):
        self.getter = getter
        self.getAll = getAll

    def __getattr__(self, item):
        return self.getter(item)

    def __getitem__(self, item):
        return self.getter(item)

    def __iter__(self):
        return iter(self.getAll())

class RobotBridge:
    def __init__(self, env, urdf_path, srdf_path):
        self.env = env

        # Load robot
        module = rave.RaveCreateModule(self.env, 'urdf')
        robot_name = module.SendCommand('load {} {}'.format(urdf_path, srdf_path))

        self.robot = self.env.GetRobot(robot_name)

        self.manip = AttributePassthrough(self.robot.GetManipulator, self.robot.GetManipulators)


    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return

    
