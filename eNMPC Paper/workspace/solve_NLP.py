from pyomo.environ import SolverFactory, TerminationCondition

# def solve_problem(step, m):
#     solver = SolverFactory("ipopt")
#     print("iteration = ", step, ", Solve " + m.name)
#     results = solver.solve(m, tee = True)
#     #check results
#     if results["Solver"][0]["Message"][-22:] != 'Optimal Solution Found':
#         raise RuntimeError("Fail to solve " + m.name + " @ " + str(step))

def solve_problem(step, solver, m):
    print("")
    print("iteration = ", step, ", Solve " + m.name)
    results = solver.solve(m, tee = True)
    #check results
    assert results.solver.termination_condition == TerminationCondition.optimal
