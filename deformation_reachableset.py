from deformation import *
from util_mvc import *
import sys

class DeformationReachableSet(Deformation):

        ## change only traj_deformed here
        Nc_handle = []

        ###############################################################
        ### LAMBDA1: move against force
        ### LAMBDA2: project onto reachable set
        ### LAMBDA3: orientation
        ### LAMBDA4: stretching
        ###############################################################

        #lambda_1 = 0.001
        #lambda_1 = 0.0005
        #lambda_1 = 0.0005
        #lambda_2 = 1
        #lambda_3 = 0.5*1e-2


        #lambda_1 = 0.001
        #lambda_2 = 0.01
        #lambda_3 = 0.1*1e-2
        #lambda_4 = 0.01

        lambda_1 = 0.0
        lambda_2 = 1.0
        lambda_3 = 0.0
        lambda_4 = 0.0

        smoothing_factor = 20.0

        def deform_onestep(self, computeNewCriticalPoint = True):

                DeformInfo = self.extractInfoFromTrajectory(self.traj_deformed)

                Ndim = DeformInfo['Ndim']
                Nwaypoints = DeformInfo['Nwaypoints']
                traj = DeformInfo['traj']
                Wori = DeformInfo['Wori']
                dWori = DeformInfo['dWori']
                eta = DeformInfo['eta']

                dU = np.zeros((Ndim,Nwaypoints))

                ###############################################################
                ## check if path dynamically feasible => return if yes
                ###############################################################

                if self.IsDynamicallyFeasible(DeformInfo):
                        return DEFORM_SUCCESS

                #print "X2=np.array(",Wori[0:2,0:100].tolist(),")"

                from deformation_module_counterwrench import *
                from deformation_module_projection_reachable_set2 import *
                from deformation_module_stretch import *
                from deformation_module_orientation import *
                from deformation_module_endpoint_projection import *

                d1 = DeformationModuleCounterWrench( DeformInfo )
                d2 = DeformationModuleProjectionReachableSet2( DeformInfo )
                d3 = DeformationModuleOrientation( DeformInfo )
                #d4 = DeformationModuleStretch( DeformInfo )

                dU += d1.get_update( self.lambda_1 )
                dU += d2.get_update( self.lambda_2 )
                dU += d3.get_update( self.lambda_3 )
                #dU += d4.get_update( self.lambda_4 )

                DeformInfo['dU'] = dU
                dend = DeformationModuleEndPointProjection( DeformInfo )
                dU += dend.get_update(0)

                #################################################################
                ## update 
                #################################################################
                self.SafetyCheckUpdate(dU)
                Wnext = Wori + eta*dU

                if np.linalg.norm(Wori-Wnext)<1e-10:
                        print "no deformation achieved with current critical point"
                        return DEFORM_NOPROGRESS

                if self.COLLISION_ENABLED:
                        if self.traj_deformed.IsInCollision(self.env, Wnext):
                                print "no deformation possible -> collision"
                                return DEFORM_COLLISION

                #Wnext = traj.InsertMissingWaypoints(Wnext, traj.DISCRETIZATION_TIME_STEP/6)
                Wnext = traj.RepairTrajectory(Wnext,1e-3)
                self.traj_deformed.new_from_waypoints(Wnext)
                #[Wori,dWori,ddWori] = self.traj_deformed.get_waypoints_second_order()
                #print "X1=np.array(",Wori[0:2,0:100].tolist(),")"
                return DEFORM_OK

