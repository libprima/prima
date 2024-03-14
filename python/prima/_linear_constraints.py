import numpy as np
from scipy.optimize import LinearConstraint

def process_single_linear_constraint(constraint):
  # Convert lb and ub to vectors if they are scalars
  num_constraints = constraint.A.shape[0]
  try:
    len_lb = len(constraint.lb)
  except TypeError:
    constraint.lb = [constraint.lb]*num_constraints
    len_lb = num_constraints
  if len_lb != num_constraints and len_lb == 1:
    constraint.lb = [constraint.lb[0]]*num_constraints
  elif len_lb != num_constraints:
    raise ValueError('Length of lb must match number of rows in A')
  # Now ub
  try:
    len_ub = len(constraint.ub)
  except TypeError:
    constraint.ub = [constraint.ub]*num_constraints
    len_ub = num_constraints
  if len_ub != num_constraints and len_ub == 1:
    constraint.ub = [constraint.ub[0]]*num_constraints
  elif len_ub != num_constraints:
    raise ValueError('Length of ub must match number of rows in A')
  return constraint


def process_multiple_linear_constraints(constraints):
  # Need to combine A, and also upgrade lb, ub to vectors if they are scalars
  constraint = process_single_linear_constraint(constraints[0])
  full_A = constraint.A
  full_lb = constraint.lb
  full_ub = constraint.ub
  for constraint in constraints[1:]:
    constraint = process_single_linear_constraint(constraint)
    full_A = np.concatenate((full_A, constraint.A), axis=0)
    full_lb = np.concatenate((full_lb, constraint.lb), axis=0)
    full_ub = np.concatenate((full_ub, constraint.ub), axis=0)
  return LinearConstraint(full_A, full_lb, full_ub)


def separate_LC_into_eq_and_ineq(linear_constraint):
  # Two things:
    # 1. PRIMA prefers A <= b as opposed to lb <= A <= ub
    # 2. PRIMA has both A_eq and A_ineq (and b_eq and b_ineq)
    # As such, we must:
    # 1. Convert lb <= A <= ub to A <= b
    # 2. Split A <= b into A_eq == b_eq and A_ineq <= b_ineq
    # Fortunately we can do both at the same time
    A_eq = []
    b_eq = []
    A_ineq = []
    b_ineq = []
    for i in range(len(linear_constraint.lb)):
      if linear_constraint.lb[i] == linear_constraint.ub[i]:
        A_eq.append(linear_constraint.A[i])
        b_eq.append(linear_constraint.lb[i])
      else:
        A_ineq.append(linear_constraint.A[i])
        b_ineq.append(linear_constraint.ub[i])
        # Flip the lb to to take format preferred by PRIMA, as long as it's not -inf
        if linear_constraint.lb[i] != -np.inf:
          A_ineq.append( - linear_constraint.A[i])
          b_ineq.append( - linear_constraint.lb[i])
    # Convert to numpy arrays, or set to None if empty
    A_eq = np.array(A_eq, dtype=np.float64) if len(A_eq) > 0 else None
    b_eq = np.array(b_eq, dtype=np.float64) if len(b_eq) > 0 else None
    A_ineq = np.array(A_ineq, dtype=np.float64) if len(A_ineq) > 0 else None
    b_ineq = np.array(b_ineq, dtype=np.float64) if len(b_ineq) > 0 else None
    return A_eq, b_eq, A_ineq, b_ineq
