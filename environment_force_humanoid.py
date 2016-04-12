from environment_force import *
import numpy as np
import math
import openravepy as rave

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


class EnvironmentHumanoid(ForceEnvironment):
        def __init__(self):
                ForceEnvironment.__init__(self)
                xmlenv='environments/plane.env.xml'
                #xmlenv='robots/escher/ground_plane_2.xml'
                urdf = 'robots/escher/escher_v2.kinbody.urdf'
                srdf = 'robots/escher/escher.robot.srdf'
                with self.env:
                        self.env_xml = xmlenv
                        self.robot_urdf = urdf
                        self.robot_srdf = srdf
                        #self.env.Add(self.env.ReadRobotXMLFile(self.robot_xml))
                        module = rave.RaveCreateModule(self.env, 'urdf')
                        self.robot_name = module.SendCommand('load {} {}'.format(urdf, srdf))
                        self.robot = self.env.GetRobot(self.robot_name)
                        self.manip = AttributePassthrough(self.robot.GetManipulator, self.robot.GetManipulators)
                        self.manip.l_arm.SetLocalToolDirection(np.array([1, 0, 0]))
                        self.manip.l_arm.SetLocalToolTransform(np.array([
                            [-1,  0, 0, -0.036],
                            [ 0, -1, 0, 0],
                            [ 0,  0, 1, 0],
                            [ 0,  0, 0, 1]])
                        )

                        self.manip.r_arm.SetLocalToolDirection(np.array([1, 0, 0]))
                        self.manip.r_arm.SetLocalToolTransform(np.array([
                            [ 1,  0, 0, 0.036],
                            [ 0,  1, 0, 0],
                            [ 0,  0, 1, 0],
                            [ 0,  0, 0, 1]])
                        )

                        self.manip.l_leg.SetLocalToolDirection(np.array([0, 0, -1]))
                        self.manip.r_leg.SetLocalToolDirection(np.array([0, 0, -1]))
                        l_arm_indices = self.robot.GetManipulator('l_arm').GetArmIndices()
                        r_arm_indices = self.robot.GetManipulator('r_arm').GetArmIndices()
                        l_leg_indices = self.robot.GetManipulator('l_leg').GetArmIndices()
                        r_leg_indices = self.robot.GetManipulator('r_leg').GetArmIndices()
                        additional_active_DOFs = ['x_prismatic_joint','y_prismatic_joint',
                                                'z_prismatic_joint','roll_revolute_joint',
                                                'pitch_revolute_joint' ,'yaw_revolute_joint',
                                                'waist_yaw']
   
                        additional_active_DOF_indices = [None]*len(additional_active_DOFs)
                        for index,j in enumerate(additional_active_DOFs):
                                additional_active_DOF_indices[index] = self.robot.GetJoint(j).GetDOFIndex()
    
                        robot_z = 0.9
                        whole_body_indices = np.concatenate((l_arm_indices, r_arm_indices, l_leg_indices, r_leg_indices, additional_active_DOF_indices),axis=0)

                        self.robot.SetActiveDOFs(whole_body_indices)
                        self.robot.SetTransform(np.array([[1,0,0,0],[0,1,0,0],[0,0,1,robot_z],[0,0,0,1]]))


                        # surrender posture
                        DOFValues = self.robot.GetDOFValues()
                        DOFValues[6] = math.pi/2
                        DOFValues[19] = math.pi/2
                        DOFValues[24] = -math.pi/2
                        DOFValues[37] = -math.pi/2
                        self.robot.SetDOFValues(DOFValues)

                        OriginalDOFValues = self.robot.GetDOFValues()

                        self.env.Load(self.env_xml)

                        #print DOFValues
                        #print OriginalDOFValues
                        #sys.exit(0)


        def GetCells(self):
                C = self.GetCellsAll()
                self.cells = C[0:2]
                return self.cells

        def GetForces(self):
                ##
                self.forces = np.array((0.0,0.0,0.0))
                self.forces = np.vstack([self.forces,(0.0,1.0,0.0)])
                return self.forces

        def RobotGetInitialPosition(self):
                return [-2.5,0.0,0.1,-pi,0,0,0,0]

        def RobotGetGoalPosition(self):
                return [-4.5,-0.0,0.1,-pi,0,0,0,0]

        def GetRobot(self):
                with self.env:
                        self.env.GetViewer().SetCamera([
                            [0.31499128, 0.09759726, -0.94406317, 6.81987572],
                            [0.94805905, 0.01409698, 0.31778187, -2.29564428],
                            [0.04432308, -0.99512615, -0.08808754, 1.60788679],
                            [0., 0., 0., 1.]
                        ])
                        return self.robot

if __name__ == "__main__":
        env = EnvironmentHumanoid()
