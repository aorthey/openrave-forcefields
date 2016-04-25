import numpy as np
import os.path

def COM_from_path(rave_path, robot, env):
        active_dofs = robot.GetActiveConfigurationSpecification()
        N = len(robot.GetActiveDOFValues())
        M = rave_path.GetNumWaypoints()
        COM_original = np.zeros((3,M))
        COM_gik = np.zeros((3,M))
        q_original = np.zeros((N, M))
        q_gik = np.zeros((N, M))

        i = 0
        [qlimL,qlimU]=robot.GetActiveDOFLimits()

        print "Waypoints:",M," - Dimension:",N
        with env.env:
                while i < M:
                        q = rave_path.GetWaypoint(i,active_dofs)
                        q_original[:,i] = env.EnforceLimits(q,qlimL,qlimU)
                        robot.SetActiveDOFValues(q_original[:,i])
                        COM_original[:,i] = robot.GetCenterOfMass()
                        i = i+1

        return [q_original, COM_original]

def GIK_from_COM(COM_path, q_original, robot, env, recompute=False, DEBUG=False):
        q_gik_fname = 'tmp/q_gik.numpy'
        COM_gik_fname = 'tmp/COM_gik.numpy'

        if not os.path.isfile(q_gik_fname+'.npy') or recompute:
                i = 0
                while i < M:
                        if DEBUG:
                                print "------------------------------------------------------------------"
                                print "------------ WAYPOINT",i,"/",M,": recompute GIK ------------ "
                                print "------------------------------------------------------------------"
                        try:
                                with env.env:
                                        robot.SetActiveDOFValues(q_original[:,i])

                                        left_leg_tf = robot.GetManipulator('l_leg').GetTransform()
                                        right_leg_tf = robot.GetManipulator('r_leg').GetTransform()
                                        left_arm_tf = robot.GetManipulator('l_arm').GetTransform()
                                        right_arm_tf = robot.GetManipulator('r_arm').GetTransform()

                                        #maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf),('l_arm',left_arm_tf),('r_arm',right_arm_tf)]
                                        maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf)]
                                        friction_coefficient = 0.8
                                        support_list = [('l_leg',friction_coefficient),
                                                        ('r_leg',friction_coefficient),
                                                        ('l_arm',friction_coefficient),
                                                        ('r_arm',friction_coefficient)]

                                        cog = COM_offset[:,i]
                                        q_res = cbirrt.DoGeneralIK(
                                                        movecog=cog,
                                                        maniptm=maniptm_list,
                                                        support=support_list,
                                                        printcommand=True)

                                        if q_res is None:
                                                print "------------ WAYPOINT",i,"/",M,": recompute GIK ------------ "
                                                print "No solution found GIK"
                                                sys.exit(0)
                                        else:
                                                q_gik[:,i] = q_res
                                        #print "DIST q,q_gik:",np.linalg.norm(q_original[:,i]-q_gik[:,i])
                                        robot.SetActiveDOFValues(q_gik[:,i])
                                        COM_gik[:,i] = robot.GetCenterOfMass()

                        except Exception as e:
                                print "Exception in GIK, waypoint",i,"/",M
                                print e
                                sys.exit(0)

                        if DEBUG:
                                dcc = np.linalg.norm(cog[:,i]-COM_gik[:,i])
                                print "ERROR GIK      :",dcc
                                print "INPUT GIK  COM :",cog
                                print "OUTPUT GIK COM :",COM_gik[:,i]

                        i = i+1

                np.save(q_gik_fname,q_gik)
                np.save(COM_gik_fname,COM_gik)
        else:
                q_gik = np.load(q_gik_fname+'.npy')
                COM_gik = np.load(COM_gik_fname+'.npy')

        return [q_gik, COM_gik]
