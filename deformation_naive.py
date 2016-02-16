from deformation_factory import *

class DeformationNaive(Deformation):

        lambda_position_peak = 0.25
        lambda_position_sigma = 10.0

        lambda_endpoint_peak = 1.0
        lambda_endpoint_sigma = 0.5

        lambda_stretch_peak = 1.0
        lambda_stretch_sigma = 5.0

        Ninsert = 20
        lambda_insert = 0.2

        ## change only traj_deformed here
        def deform_onestep(self):

                [Wori,dWori] = self.traj_deformed.get_waypoints(N=100)#(self.traj_ori, N = 100)
                Ndim = Wori.shape[0]

                Nwaypoints = Wori.shape[1]
                F = self.GetForcesAtWaypoints(Wori)
                print F.shape

                idx_invalid = self.GetFirstInfeasibleWaypoint(Wori, dWori, F)
                dW_invalid = dWori[:,idx_invalid]/(np.linalg.norm(dWori[:,idx_invalid]))
                print "############################################################"
                print "INVALID POINT FOUND"
                print idx_invalid
                print Wori[:,idx_invalid], dW_invalid
                print "############################################################"

                Wdown = Wori[:,0:idx_invalid]
                Wup = Wori[:,idx_invalid:]

                W = np.zeros((Ndim,self.Ninsert+Nwaypoints))

                W[:,0:idx_invalid] = Wdown
                W[:,idx_invalid+self.Ninsert:] = Wup

                Wmiddle = np.zeros((Ndim,self.Ninsert))

                for i in range(0,self.Ninsert):
                        Wstretch = (self.Ninsert-(i))*self.lambda_insert*dW_invalid
                        Wmiddle[:,i] = Wori[:,idx_invalid] - Wstretch

                W[:,idx_invalid:idx_invalid+self.Ninsert] = Wmiddle

                Wupdate_stretch = self.GetForceStretch(idx_invalid, W, -(self.Ninsert+1)*self.lambda_insert*dW_invalid)
                Wupdate_position = self.GetForcePosition(idx_invalid+self.Ninsert, W, -F[:,idx_invalid])

                Wupdate = Wupdate_position+Wupdate_stretch

                Winit=self.env.RobotGetInitialPosition()
                Wgoal=self.env.RobotGetGoalPosition()

                Wupdate_correction = self.GetForceEndpointCorrection(Wupdate, W, Winit[:Ndim], Wgoal[:Ndim])

                W = W + Wupdate_correction

                self.traj_deformed.new_from_waypoints(W)


        def gaussian(self, a, b, c, t):
                return a*np.exp(-((t-b)*(t-b))/(2*c*c))


        def lambda_position(self, t0, tidx):
                a = self.lambda_position_peak
                c = self.lambda_position_sigma
                #a = 0.25
                #c = 15.0
                return self.gaussian(a,0,c,abs(t0-tidx))

        def lambda_endpoint(self, t0, endT):
                #a = 1.0
                #c = 0.5
                a = self.lambda_endpoint_peak
                c = self.lambda_endpoint_sigma
                return self.gaussian(a,0,c,abs(t0-endT))

        def lambda_stretch(self, t0, tidx):
                ## a has to be 1.0 to make sure that we stretch the last
                ##point the correct amount
                #a = 1.0 
                #c = 10.0
                a = self.lambda_stretch_peak
                c = self.lambda_stretch_sigma
                return self.gaussian(a,0,c,abs(t0-tidx))

        ## idx: index at which waypoint will be shifted
        ## W  : array of waypoints
        ## dF : direction of force applied at W[:,idx]
        ## 
        ## return: a force direction to each waypoint
        def GetForcePosition(self, idx, W, dF):
                Ndim = W.shape[0]
                N = W.shape[1]
                dW = np.zeros((Ndim,N))
                for i in range(0,N):
                        dW[:,i] = self.lambda_position(i, idx)*dF

                return dW

        def GetForceCellCorrection(self, W, Cells, dF):

                Ndim = W.shape[0]
                N = W.shape[1]
                dW = np.zeros((Ndim,N))
                for i in range(0,N):
                        dW[:,i] = self.lambda_position(i, idx)*dF
                return dW

        def GetForceStretch(self, idx, W, dF):

                if idx < 1:
                        print "\nERROR: initial point was not feasible, cannot \n stretch this without changing initial \n point. Solution: Change your trajectory \n to make sure that your initial point is feasible\n"
                        sys.exit(0)

                assert(idx>=1)

                Ndim = W.shape[0]
                N = W.shape[1]
                dW = np.zeros((Ndim,N))

                for i in range(0,idx):
                        dW[:,i] = self.lambda_stretch(i, idx)*dF
                return dW

        ### remove any updates on the endpoint, keep them static
        def GetForceEndpointCorrection(self, Wupdate, W, Winit, Wgoal):

                Ndim = W.shape[0]
                N = W.shape[1]
                dW = np.zeros((Ndim,N))

                assert(W.shape[1]==Wupdate.shape[1])

                for i in range(0,N):
                        dW[:,i] = dW[:,i] + self.lambda_endpoint(i, 0)*(-Wupdate[:,0])
                        dW[:,i] = dW[:,i] + self.lambda_endpoint(i, N-1)*(-Wupdate[:,-1])

                if np.linalg.norm(dW[:,-1])>0.001:
                        print "update correction not successful!"
                        sys.exit(0)

                dWinit = Winit - W[:,0]
                dWgoal = Wgoal - W[:,-1]

                for i in range(0,N):
                        dW[:,i] = dW[:,i] + self.lambda_endpoint(i, 0)*(dWinit)
                        dW[:,i] = dW[:,i] + self.lambda_endpoint(i, N-1)*(dWgoal)

                Wcorrection = Wupdate + dW
                return Wcorrection


