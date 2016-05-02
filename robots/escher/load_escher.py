from RobotBridge import RobotBridge
import numpy as np
import openravepy as rave

urdf = 'file://escher_v2.kinbody.urdf'
srdf = 'file://escher.robot.srdf'

#initialize link configuration
def escher(env):
    escher = RobotBridge(env, urdf_path=urdf, srdf_path=srdf)

    escher.manip.l_arm.SetLocalToolDirection(np.array([1, 0, 0]))
    escher.manip.l_arm.SetLocalToolTransform(np.array([
        [-1,  0, 0, -0.036],
        [ 0, -1, 0, 0],
        [ 0,  0, 1, 0],
        [ 0,  0, 0, 1]])
    )

    escher.manip.r_arm.SetLocalToolDirection(np.array([1, 0, 0]))
    escher.manip.r_arm.SetLocalToolTransform(np.array([
        [ 1,  0, 0, 0.036],
        [ 0,  1, 0, 0],
        [ 0,  0, 1, 0],
        [ 0,  0, 0, 1]])
    )

    escher.manip.l_leg.SetLocalToolDirection(np.array([0, 0, -1]))
    escher.manip.r_leg.SetLocalToolDirection(np.array([0, 0, -1]))

    return escher


