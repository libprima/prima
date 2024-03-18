Requirements
==============

These are high level requirements for the python package. Next to the requirement is the name of the test file and the name of the test
within that file which implements it. The testfile/testname is written in such a way that it can be passed as an argument to `pytest` in
order to run that test by itself.


## Requirements for basic functionality, algorithm autoselection, and explicit algorithm selection
- [x] providing nonlinear constraints alone calls COBYLA and the constraints are successfully applied - test_basic_functionality.py::test_provide_nonlinear_constraints_alone
- [x] providing nonlinear constraints alone and selecting COBYLA calls COBYLA and the constraints are successfully applied - test_basic_functionality.py::test_provide_nonlinear_constraints_alone_and_select_COBYLA
- [x] providing linear constraints alone calls LINCOA and the constraints are successfully applied - test_basic_functionality.py::test_provide_linear_constraints_alone
- [x] providing linear constraints alone and selecting LINCOA calls LINCOA and the constraints are successfully applied - test_basic_functionality.py::test_provide_linear_constraints_alone_and_select_LINCOA
- [x] providing bounds alone calls BOBYQA and the bounds are successfully applied - test_basic_functionality.py::test_provide_bounds_alone
- [x] providing bounds alone and selecting BOBYQA calls BOBYQA and the bounds are successfully applied - test_basic_functionality.py::test_provide_bounds_alone_and_select_BOBYQA
- [x] not providing any bounds, linear constraints, or nonlinear constraints calls NEWUOA and provides the optimal unbounded/unconstrained solution - test_basic_functionality.py::test_not_providing_bounds_linear_constraints_or_nonlinear_constraints
- [x] not providing any bounds, linear constraints, or nonlinear constraints and selecting NEWUOA calls NEWUOA and provides the optimal unbounded/unconstrained solution - test_basic_functionality.py::test_not_providing_bounds_linear_constraints_or_nonlinear_constraints_and_select_NEWUOA
- [x] not providing any bounds, linear constraints, or nonlinear constraints and selecting UOBYQA calls UOBYQA and provides the optimal unbounded/unconstrained solution - test_basic_functionality.py::test_not_providing_bounds_linear_constraints_or_nonlinear_constraints_and_select_UOBYQA


## Requirements for combining constraints and bounds
- [x] providing linear and nonlinear constraints together calls COBYLA and the constraints are successfully applied - test_combining_constraints.py::test_providing_linear_and_nonlinear_constraints
- [x] providing bounds and linear constraints together calls LINCOA and the constraints/bounds are successfully applied - test_combining_constraints.py::test_providing_bounds_and_linear_constraints
- [x] providing bounds and nonlinear constraints together calls COBYLA and the constraints/bounds are successfully applied - test_combining_constraints.py::test_providing_bounds_and_nonlinear_constraints
- [x] providing bounds and linear and nonlinear constraints together calls COBYLA and the constraints/bounds are successfully applied - test_combining_constraints.py::test_providing_bounds_and_linear_and_nonlinear_constraints


## Requirements for scalar/list/array
- [x] providing a nonlinear constraint function returning a numpy array works without warnings or errors - test_data_types.py::test_constraint_function_returns_numpy_array
- [x] providing a nonlinear constraint function returning a Python list works without warnings or errors - test_data_types.py::test_constraint_function_returns_list
- [x] providing a nonlinear constraint function returning a scalar works without warnings or errors - test_data_types.py::test_constraint_function_returns_scalar
- [x] providing x0 as numpy array works without warnings or errors - test_data_types.py::test_x0_as_list
- [x] providing x0 as Python list works without warnings or errors - test_data_types.py::test_x0_as_array
- [x] providing x0 as scalar works without warnings or errors - test_data_types.py::test_x0_as_scalar
- [x] providing A as either a scalar, list, or numpy array, providing lb as either a scalar, list, or numpy array, and providing ub as either a scalar, list or numpy array works without warning or errors in all combinations - test_data_types.py::test_linear_constraint_data_types
- [x] providing a nonlinear constraint with scalar lower bound and scalar upper bound works without warnings or errors - test_data_types.py::test_nonlinear_constraint_lb_scalar_ub_scalar
- [x] providing a nonlinear constraint with scalar lower bound and list/array upper bound works without warnings or errors - test_data_types.py::test_nonlinear_constraint_lb_scalar_ub_not_scalar
- [x] providing a nonlinear constraint with list/array lower bound and scalar upper bound works without warnings or errors - test_data_types.py::test_nonlinear_constraint_lb_not_scalar_ub_scalar
- [x] providing a nonlinear constraint with list/array lower bound and list/array upper bound works without warnings or errors - test_data_types.py::test_nonlinear_constraint_lb_not_scalar_ub_not_scalar


## Requirements for various options
- [x] various options can be successfully provided
  - [x] ftarget - test_options.py::test_ftarget
  - [x] iprint - test_options.py::test_iprint
  - [x] max function evaluations via maxfun (name used by SciPy) - test_options.py::test_maxfun
  - [x] max function evaluations via maxfev (name used by PRIMA) - test_options.py::test_maxfev
  - [x] npt - test_options.py::test_npt
  - [x] rhobeg - test_options.py::test_rhobeg
  - [x] rhoend - test_options.py::test_rhoend
  - [x] quiet - test_options.py::test_quiet
- [x] an objective function without args can be used successfully - test_options.py::test_normal_function
- [x] an objective function with args can be used successfully - test_options.py::test_function_with_regular_args
- [x] an objective function with *args can be used successfully - test_options.py::test_function_with_star_args
- [x] providing a callback leads to the callback being successfully called - test_options.py::test_callback
- [x] providing a callback that returns True leads to early termination - test_options.py::test_callback_early_termination


## Requirements for compatibility with existing APIs
- [x] compatible with scipy.optimize.minimize API - test_compatibility_scipy.py::test_scipy
- [x] compatible with PDFO API - test_compatibility_pdfo.py::test_pdfo


## Requirements for processing nonlinear constraints
- [x] providing a list of nonlinear constraints with either scalar or vector for lb/ub correctly constructs the constraint function - test_process_nonlinear_constraints.py::test_multiple_nl_constraints_various_data_types
- [x] providing a single nonlinear constraints with either scalar or vector for lb/ub correctly constructs the constraint function - test_process_nonlinear_constraints.py::test_single_nl_constraint
- [x] providing a nonlinear constraint with upper and/or lower bounds as vectors of different length than the length of the output of the constraint function raises an exception - test_process_nonlinear_constraints.py::test_length_nlc_fun_not_equal_to_length_lb_ub
- [x] providing a nonlinear constraint with upper and lower bound as vectors of where upper bounds has correct length and lower bounds does not as compared to the length of the output of the constraint function raises an exception - test_process_nonlinear_constraints.py::test_length_nlc_fun_ne_to_length_ub
- [x] providing a nonlinear constraint with scalar lb at -inf and scalar ub at inf raises an exception - test_process_nonlinear_constraints.py::test_lb_neg_inf_ub_inf_raises
- [x] providing a nonlinear constraint with scalar lb at -inf and vector ub containing at least 1 inf raises an exception - test_process_nonlinear_constraints.py::test_lb_neg_inf_ub_vector_w_inf_raises
- [x] providing a nonlinear constraint with vector lb containing at least 1 -inf and scalar ub at inf raises an exception - test_process_nonlinear_constraints.py::test_lb_vector_with_neg_inf_ub_inf_raises
- [x] providing a nonlinear constraint with vector lb and vector ub containing -inf and +inf at the same index raises an exception - test_process_nonlinear_constraints.py::test_lb_vector_with_neg_inf_ub_vector_w_inf_at_same_index_raises


## Requirements for processing linear constraints
- [x] providing a list of linear constraints leads to them being successfully combined - test_process_linear_constraints.py::test_multiple_linear_constraints_implementation
- [x] providing a list of linear constraints leads to them being successfully combined and applied - test_process_linear_constraints.py::test_multiple_linear_constraints_high_level
- [x] providing linear constraints that contain both equality and inequality constraints leads to a correct separation into A_eq/b_eq and A_ineq/b_ineq - test_process_linear_constraints.py::test_separate_LC_into_eq_and_ineq


## Requirements for constraint types
- [x] PRIMA's LinearConstraint is the same as SciPy's scipy.optimize.LinearConstraint - test_constraint_types.py::test_prima_lc_is_scipy_lc
- [x] providing a linear constraint as a type scipy.optimize.LinearConstraint works without warnings or errors - test_constraint_types.py::test_linear_constraint_object
- [x] providing a linear constraint as a dictionary with keys 'A', 'lb', and 'ub' works without warnings or errors - test_constraint_types.py::test_linear_constraint_dict
- [x] PRIMA's NonlinearConstraint is the same as SciPy's scipy.optimize.NonlinearConstraint - test_constraint_types.py::test_prima_lc_is_scipy_lc
- [x] providing a nonlinear constraint as a type scipy.optimize.NonlinearConstraint works without warnings or errors - test_constraint_types.py::test_nonlinear_constraint_object
- [x] providing a nonlinear constraint as a dictionary with keys 'fun', 'lb', and 'ub' works without warnings or errors - test_constraint_types.py::test_nonlinear_constraint_dict
- [x] providing an unsupported constraint type raises an exception - test_constraint_types.py::test_unsupported_type_raises_exception


## Requirements for warnings and exceptions related to method selection
- [x] providing a method that is not one of the 5 provided raises an exception - test_method_selection_warnings_and_exceptions.py::test_method_not_recognized
- [x] providing nonlinear constraints to a method that cannot handle them raises an exception - test_method_selection_warnings_and_exceptions.py::test_providing_nonlinear_constraints_to_non_cobyla_method
- [x] providing linear constraints to a method that cannot handle them raises an exception - test_method_selection_warnings_and_exceptions.py::test_providing_linear_constraints_to_non_cobyla_non_lincoa_method
- [x] providing bounds to a method that cannot handle them raises an exception - test_method_selection_warnings_and_exceptions.py::test_providing_bounds_to_non_cobyla_non_lincoa_non_bobyqa_method

## Requirements for regression tests
These tests are for behavior observed during testing that we want to make sure remains fixed.
- [x] providing anonymous lambda as: objective function, constraint function, callback works without warnings or errors - test_anonymous_lambda.py::test_anonymous_lambda
- [x] calling with a method that is not available throws an exception and exits cleanly (i.e. no other warnings or error or hanging of the interpreter) (relates to the test regarding anonymous lambda functions) - test_anonymous_lambda.py::test_anonymous_lambda_unclean_exit
