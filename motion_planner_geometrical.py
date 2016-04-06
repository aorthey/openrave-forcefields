from motion_planner_factory import MotionPlanner

class MotionPlannerGeometrical(MotionPlanner):
        def GetPath(self):
                planner_name = 'birrt'
                return self.PlanPath(planner_name)
