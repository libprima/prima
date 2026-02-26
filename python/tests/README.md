Requirements
==============

These are high level requirements for the python package. Next to the requirement is the name of the test file and the name of the test
within that file which implements it. The testfile/testname is written in such a way that it can be passed as an argument to `pytest` in
order to run that test by itself.

Not all tests apply to all backends. In particular, the Python backend only implements COBYLA and so cannot participate in any of the tests related to non-COBLYA algorithms. Each requirement will be prefaced with an array specifying the backends to which it applies. `[F, P]` indicates that it applies to both the Fortran and Python backend and should be implemented as such. `[F]` or `[P]` implies that it only applies to the Fortran or to the Python backend. `[N/A]` implies that the test is focused on the interface logic and the choice of backend is irrelevant.


## Requirements for basic functionality, algorithm autoselection, and explicit algorithm selection
- [x] [F, P] providing nonlinear constraints alone calls COBYLA and the constraints are successfully applied - test_basic_functionality.py::test_provide_nonlinear_constraints_alone
- [x] [F, P] providing nonlinear constraints alone and selecting COBYLA calls COBYLA and the constraints are successfully applied - test_basic_functionality.py::test_provide_nonlinear_constraints_alone_and_select_COBYLA
- [x] [F] providing linear constraints alone calls LINCOA and the constraints are successfully applied - test_basic_functionality.py::test_provide_linear_constraints_alone
- [x] [F] providing linear constraints alone and selecting LINCOA calls LINCOA and the constraints are successfully applied - test_basic_functionality.py::test_provide_linear_constraints_alone_and_select_LINCOA
- [x] [F] providing bounds alone calls BOBYQA and the bounds are successfully applied - test_basic_functionality.py::test_provide_bounds_alone
- [x] [F] providing bounds alone and selecting BOBYQA calls BOBYQA and the bounds are successfully applied - test_basic_functionality.py::test_provide_bounds_alone_and_select_BOBYQA
- [x] [F] not providing any bounds, linear constraints, or nonlinear constraints calls NEWUOA and provides the optimal unbounded/unconstrained solution - test_basic_functionality.py::test_not_providing_bounds_linear_constraints_or_nonlinear_constraints
- [x] [F] not providing any bounds, linear constraints, or nonlinear constraints and selecting any algorithm calls that algorithm and provides the optimal unbounded/unconstrained solution - test_basic_functionality.py::test_not_providing_bounds_linear_constraints_or_nonlinear_constraints_and_selecting_any_algorithm
- [x] [F, P] selecting a backend which is not recognized leads to an exception being raised - test_basic_functionality.py::test_invalid_backend
- [x] [P] Selecting an algorithm which is not implemented by the Python backend leads to a warning and auto-selection of the Fortran backend. - test_basic_functionality.py::test_select_algorithm_not_provided_by_python_implementation


## Requirements for combining constraints and bounds
- [x] [F, P] providing linear and nonlinear constraints together calls COBYLA and the constraints are successfully applied - test_combining_constraints.py::test_providing_linear_and_nonlinear_constraints
- [x] [F] providing bounds and linear constraints together calls LINCOA and the constraints/bounds are successfully applied - test_combining_constraints.py::test_providing_bounds_and_linear_constraints
- [x] [F, P] providing bounds and nonlinear constraints together calls COBYLA and the constraints/bounds are successfully applied - test_combining_constraints.py::test_providing_bounds_and_nonlinear_constraints
- [x] [N/A] providing bounds and linear and nonlinear constraints together calls COBYLA and the constraints/bounds are successfully applied - test_combining_constraints.py::test_providing_bounds_and_linear_and_nonlinear_constraints


## Requirements for scalar/list/array
- [x] [F, P] providing a nonlinear constraint function returning a numpy array works without warnings or errors - test_data_types.py::test_constraint_function_returns_numpy_array
- [x] [F, P] providing a nonlinear constraint function returning a Python list works without warnings or errors - test_data_types.py::test_constraint_function_returns_list
- [x] [F, P] providing a nonlinear constraint function returning a scalar works without warnings or errors - test_data_types.py::test_constraint_function_returns_scalar
- [x] [F, P] providing x0 as numpy array works without warnings or errors - test_data_types.py::test_x0_as_list
- [x] [F, P] providing x0 as Python list works without warnings or errors - test_data_types.py::test_x0_as_array
- [x] [F, P] providing x0 as scalar works without warnings or errors - test_data_types.py::test_x0_as_scalar
- [x] [F, P] providing A as either a scalar, list, or numpy array, providing lb as either a scalar, list, or numpy array, and providing ub as either a scalar, list or numpy array works without warning or errors in all combinations - test_data_types.py::test_linear_constraint_data_types
- [x] [F, P] providing a nonlinear constraint with scalar lower bound and scalar upper bound works without warnings or errors - test_data_types.py::test_nonlinear_constraint_lb_scalar_ub_scalar
- [x] [F, P] providing a nonlinear constraint with scalar lower bound and list/array upper bound works without warnings or errors - test_data_types.py::test_nonlinear_constraint_lb_scalar_ub_not_scalar
- [x] [F, P] providing a nonlinear constraint with list/array lower bound and scalar upper bound works without warnings or errors - test_data_types.py::test_nonlinear_constraint_lb_not_scalar_ub_scalar
- [x] [F, P] providing a nonlinear constraint with list/array lower bound and list/array upper bound works without warnings or errors - test_data_types.py::test_nonlinear_constraint_lb_not_scalar_ub_not_scalar


## Requirements for various options
- [x] [F, P] various options can be successfully provided
  - [x] [F, P] ftarget - test_options.py::test_ftarget
  - [x] [F, P] iprint - test_options.py::test_iprint
  - [x] [F, P] max function evaluations via maxfev (name used by PRIMA) - test_options.py::test_maxfev
  - [x] [F] npt - test_options.py::test_npt
  - [x] [F, P] rhobeg - test_options.py::test_rhobeg
  - [x] [F, P] rhoend - test_options.py::test_rhoend
  - [x] [N/A] quiet - test_options.py::test_quiet
- [x] [F, P] an objective function without args can be used successfully - test_options.py::test_normal_function
- [x] [F, P] an objective function with args can be used successfully - test_options.py::test_function_with_regular_args
- [x] [F, P] an objective function with *args can be used successfully - test_options.py::test_function_with_star_args
- [x] [F, P] providing a callback leads to the callback being successfully called - test_options.py::test_callback
- [x] [F, P] providing a callback that returns True leads to early termination - test_options.py::test_callback_early_termination

# Requirements for invalid options
These tests are only applied to the Python backend since I am unable to properly capture stdout from Fortran with these tests, despite being able to capture it in test_iprint above.
- [x] [P] Supplying invalid rhobeg should trigger a warning message - test_invalid_options.py::test_invalid_rhobeg
- [x] [P] Supplying eta2 without a corresponding eta1 and eta2 out of range should trigger a warning message - test_invalid_options.py::test_eta2_without_eta1_and_eta2_out_of_range
- [x] [P] Supplying eta1 without a corresponding eta2 and eta1 out of range should trigger a warning message - test_invalid_options.py::test_eta1_without_eta2_and_eta1_out_of_range
- [x] [P] Supplying an invalid iprint should trigger a warning message - test_invalid_options.py::test_invalid_iprint


## Requirements for compatibility with existing APIs
- [x] [N/A] compatible with scipy.optimize.minimize API - test_compatibility_scipy.py::test_scipy
- [x] [N/A] compatible with PDFO API - test_compatibility_pdfo.py::test_pdfo


## Requirements for processing nonlinear constraints
- [x] [N/A] providing a list of nonlinear constraints with either scalar or vector for lb/ub correctly constructs the constraint function - test_process_nonlinear_constraints.py::test_multiple_nl_constraints_various_data_types
- [x] [N/A] providing a single nonlinear constraints with either scalar or vector for lb/ub correctly constructs the constraint function - test_process_nonlinear_constraints.py::test_single_nl_constraint
- [x] [N/A] providing a nonlinear constraint with upper and/or lower bounds as vectors of different length than the length of the output of the constraint function raises an exception - test_process_nonlinear_constraints.py::test_length_nlc_fun_not_equal_to_length_lb_ub
- [x] [N/A] providing a nonlinear constraint with upper and lower bound as vectors of where upper bounds has correct length and lower bounds does not as compared to the length of the output of the constraint function raises an exception - test_process_nonlinear_constraints.py::test_length_nlc_fun_ne_to_length_ub
- [x] [N/A] providing a nonlinear constraint with scalar lb at -inf and scalar ub at inf raises an exception - test_process_nonlinear_constraints.py::test_lb_neg_inf_ub_inf_raises
- [x] [N/A] providing a nonlinear constraint with scalar lb at -inf and vector ub containing at least 1 inf raises an exception - test_process_nonlinear_constraints.py::test_lb_neg_inf_ub_vector_w_inf_raises
- [x] [N/A] providing a nonlinear constraint with vector lb containing at least 1 -inf and scalar ub at inf raises an exception - test_process_nonlinear_constraints.py::test_lb_vector_with_neg_inf_ub_inf_raises
- [x] [N/A] providing a nonlinear constraint with vector lb and vector ub containing -inf and +inf at the same index raises an exception - test_process_nonlinear_constraints.py::test_lb_vector_with_neg_inf_ub_vector_w_inf_at_same_index_raises


## Requirements for processing linear constraints
- [x] [F] providing a list of linear constraints leads to them being successfully combined - test_process_linear_constraints.py::test_multiple_linear_constraints_implementation
- [x] [F] providing a list of linear constraints leads to them being successfully combined and applied - test_process_linear_constraints.py::test_multiple_linear_constraints_high_level
- [x] [F] providing linear constraints that contain both equality and inequality constraints leads to a correct separation into A_eq/b_eq and A_ineq/b_ineq - test_process_linear_constraints.py::test_separate_LC_into_eq_and_ineq


## Requirements for constraint types
- [x] [N/A] PRIMA's LinearConstraint is the same as SciPy's scipy.optimize.LinearConstraint - test_constraint_types.py::test_prima_lc_is_scipy_lc
- [x] [N/A] providing a linear constraint as a type scipy.optimize.LinearConstraint works without warnings or errors - test_constraint_types.py::test_linear_constraint_object
- [x] [N/A] providing a linear constraint as a dictionary with keys 'A', 'lb', and 'ub' works without warnings or errors - test_constraint_types.py::test_linear_constraint_dict
- [x] [N/A] PRIMA's NonlinearConstraint is the same as SciPy's scipy.optimize.NonlinearConstraint - test_constraint_types.py::test_prima_lc_is_scipy_lc
- [x] [N/A] providing a nonlinear constraint as a type scipy.optimize.NonlinearConstraint works without warnings or errors - test_constraint_types.py::test_nonlinear_constraint_object
- [x] [N/A] providing a nonlinear constraint as a dictionary with keys 'fun', 'lb', and 'ub' works without warnings or errors - test_constraint_types.py::test_nonlinear_constraint_dict
- [x] [N/A] providing an unsupported constraint type raises an exception - test_constraint_types.py::test_unsupported_type_raises_exception


## Requirements for warnings and exceptions related to method selection
- [x] [N/A] providing a method that is not one of the 5 provided raises an exception - test_method_selection_warnings_and_exceptions.py::test_method_not_recognized
- [x] [N/A] providing nonlinear constraints to a method that cannot handle them raises an exception - test_method_selection_warnings_and_exceptions.py::test_providing_nonlinear_constraints_to_non_cobyla_method
- [x] [N/A] providing linear constraints to a method that cannot handle them raises an exception - test_method_selection_warnings_and_exceptions.py::test_providing_linear_constraints_to_non_cobyla_non_lincoa_method
- [x] [N/A] providing bounds to a method that cannot handle them raises an exception - test_method_selection_warnings_and_exceptions.py::test_providing_bounds_to_non_cobyla_non_lincoa_non_bobyqa_method

## Requirements for package functionality
- [x] [P] Producing the python package without the Fortran bindings should still allow for the use of the pure python backend - test_without_pybindings.py::test_without_pybindings

## Requirements for regression tests
These tests are for behavior observed during testing that we want to make sure remains fixed.
- [x] [F] providing anonymous lambda as: objective function, constraint function, callback works without warnings or errors - test_anonymous_lambda.py::test_anonymous_lambda
- [x] [F] calling with a method that is not available throws an exception and exits cleanly (i.e. no other warnings or error or hanging of the interpreter) (relates to the test regarding anonymous lambda functions) - test_anonymous_lambda.py::test_anonymous_lambda_unclean_exit
