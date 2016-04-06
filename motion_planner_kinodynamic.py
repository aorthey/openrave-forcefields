from motion_planner_factory import MotionPlanner

class MotionPlannerKinodynamic(MotionPlanner):
        def GetPath(self):
                planner_name = 'kinodynamicrrt'
                return self.PlanPath(planner_name)
