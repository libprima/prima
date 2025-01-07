from optiprofiler.problems import load_cutest_problem
from primapy import minimize, Bounds, LinearConstraint, NonlinearConstraint
import numpy as np

'''
This module tests various problem from the CUTEST set in order to really stress test
the implementation and also cover some cases not covered by the naive tests in
test_end_to_end.py. The list is semi-arbitrary, some of these helped to fine bugs
when testing the Python implementation against the Fortran one.
'''

def get_constraints(problem):
    constraints = []
    if problem.m_linear_ub > 0:
        constraints.append(LinearConstraint(problem.a_ub, -np.inf, problem.b_ub))
    if problem.m_linear_eq > 0:
        constraints.append(LinearConstraint(problem.a_eq, problem.b_eq, problem.b_eq))
    if problem.m_nonlinear_ub > 0:
        constraints.append(NonlinearConstraint(problem.c_ub, -np.inf, np.zeros(problem.m_nonlinear_ub)))
    if problem.m_nonlinear_eq > 0:
        constraints.append(NonlinearConstraint(problem.c_eq, np.zeros(problem.m_nonlinear_eq), np.zeros(problem.m_nonlinear_eq)))
    return constraints


def run_problem(name, expected_x, expected_f, expected_constraints, expected_nf):
    problem = load_cutest_problem(name)
    constraints = get_constraints(problem)
    bounds = Bounds(problem.lb, problem.ub)
    result = minimize(problem.fun, problem.x0, method='cobyla', constraints=constraints, bounds=bounds)
    assert np.allclose(result.x, expected_x, atol=1e-8)
    assert np.isclose(result.f, expected_f, atol=1e-8)
    assert np.allclose(result.constr, expected_constraints, atol=1e-8)
    assert result.nf == expected_nf


def test_errinbar():
    # Expected values are just obtained from running the problem and collecting the results
    # If future changes improve the algorithm, these values may need to be updated.
    expected_x = np.array([
            3.23416763e+00, -5.38155691e+01,  1.19635776e+03, -5.38155691e+01,
           -1.89932469e+02,  1.78041499e+02, -1.07283385e+02, -1.76693103e+01,
            1.65553655e+02, -2.37056910e+00, -2.37056910e+00,  1.11111855e+00,
           -2.37056910e+00, -2.37056910e+00, -2.37056910e+00,  1.12930585e+00,
            9.19784743e-02, -2.37056910e+00])
    expected_f = 382.65742398099866
    expected_constraints = np.array([
              3.0155691 ,    3.0155691 , -228.84149931,  -33.13068968,
           -164.90865503,    3.0155691 ,    3.0155691 ,   -0.46611855,
              3.0155691 ,    3.0155691 ,    3.0155691 ,   -0.48430585,
              0.55302153,    3.0155691 ,  -22.86      ,   -3.0154404 ,
             -3.01556889,    3.01551024,    3.0157243 ,    3.01556564,
             -3.01557263,   -3.01535278,    3.01572578,    3.0154404 ,
              3.01556889,   -3.01551024,   -3.0157243 ,   -3.01556564,
              3.01557263,    3.01535278,   -3.01572578])
    expected_nf = 9000
    run_problem('ERRINBAR', expected_x, expected_f, expected_constraints, expected_nf)


def test_palmer2c():
    expected_x = np.array([0.92409968,  0.95594354,  0.95322078,  0.74881236,  0.53994945,
                        -0.01907517,  0.05613685, -0.01960721])
    expected_f = 1402.0983572109658
    expected_constraints = np.array([])
    expected_nf = 4000
    run_problem('PALMER2C', expected_x, expected_f, expected_constraints, expected_nf)


def test_palmer3b():
    expected_x = np.array([1.11911978, 7.56621466, 0.44987163, 0.03976417])
    expected_f = 319.67553113156816
    expected_constraints = np.array([-0.44986163, -0.03975417])
    expected_nf = 2000
    run_problem('PALMER3B', expected_x, expected_f, expected_constraints, expected_nf)


def test_tfi3():
    expected_x = np.array([1, 0.5, 0])
    expected_f = 5.367003099159174
    expected_constraints = np.array([
            0.        , -0.00509999, -0.01039984, -0.01589919, -0.02159744,
           -0.02749377, -0.03358709, -0.03987611, -0.0463593 , -0.05303492,
           -0.05990099, -0.06695534, -0.07419558, -0.08161914, -0.08922323,
           -0.09700489, -0.104961  , -0.11308825, -0.12138319, -0.1298422 ,
           -0.13846154, -0.14723733, -0.15616559, -0.16524219, -0.17446294,
           -0.18382353, -0.19331959, -0.20294669, -0.2127003 , -0.22257587,
           -0.23256881, -0.24267448, -0.25288824, -0.26320543, -0.27362137,
           -0.2841314 , -0.29473088, -0.30541516, -0.31617966, -0.32701979,
           -0.33793103, -0.34890891, -0.359949  , -0.37104692, -0.38219839,
           -0.39339917, -0.4046451 , -0.4159321 , -0.42725618, -0.43861342,
           -0.45      , -0.46141219, -0.47284635, -0.48429893, -0.49576649,
           -0.50724568, -0.51873325, -0.53022605, -0.54172104, -0.55321527,
           -0.56470588, -0.57619015, -0.58766542, -0.59912914, -0.61057889,
           -0.6220123 , -0.63342714, -0.64482124, -0.65619256, -0.66753912,
           -0.67885906, -0.69015059, -0.70141201, -0.71264173, -0.7238382 ,
           -0.735     , -0.74612576, -0.7572142 , -0.76826411, -0.77927437,
           -0.7902439 , -0.80117173, -0.81205692, -0.82289863, -0.83369606,
           -0.84444848, -0.85515521, -0.86581564, -0.87642922, -0.88699542,
           -0.89751381, -0.90798397, -0.91840555, -0.92877822, -0.93910172,
           -0.94937582, -0.95960033, -0.9697751 , -0.97990002, -0.989975  ,
           -1.        ])
    expected_nf = 15
    run_problem('TFI3', expected_x, expected_f, expected_constraints, expected_nf)


def test_hs103():
    # This one hits the section in trustregion.py which scales the problem if A
    # contains large values.
    expected_x = np.array([4.39403665, 0.85447094, 2.84319504, 3.40002056, 0.72292401,
            0.8704116 , 0.02463953])
    expected_f = 543.6679593101205
    expected_constraints = np.array([
           -4.29403665e+00, -7.54470945e-01, -2.74319504e+00, -3.30002056e+00,
           -6.22924011e-01, -7.70411597e-01, -1.46395250e-02, -5.60596335e+00,
           -9.14552906e+00, -7.15680496e+00, -6.59997944e+00, -9.27707599e+00,
           -9.12958840e+00, -9.97536047e+00, -5.16980275e-11, -1.57108661e-12,
            2.90542244e-14, -2.60019228e-11, -4.43667959e+02, -2.45633204e+03])
    expected_nf = 1380
    run_problem('HS103', expected_x, expected_f, expected_constraints, expected_nf)


def test_cresc4():
    expected_x = np.array([-3.07632335e+01,  9.87373538e-01,  2.09562218e+00,  1.40743208e+01,
                            5.92899122e-04,  3.90051552e-01])
    expected_f = 7.940863533363242
    expected_constraints = np.array([
        -2.09562217e+00, -1.30743208e+01, -5.92899122e-04, -5.15523260e-05,
        -6.28259230e+00, -2.14662504e-02, -1.16793954e+02, -4.56874771e+00,
        -5.32927399e+01, -5.49304601e-01, -5.72422340e+01, -2.04024620e+00,
        -8.52807205e+01])
    expected_nf = 3000
    run_problem('CRESC4', expected_x, expected_f, expected_constraints, expected_nf)


def test_mgh10ls():
    # This one also hits the section in trustregion.py which scales the problem if A
    # contains large values.
    expected_x = np.array([1.49546689e-03, 3.99999495e+05, 2.50010062e+04])
    expected_f = 1366860355.936367
    expected_constraints = np.array([])
    expected_nf = 47
    run_problem('MGH10LS', expected_x, expected_f, expected_constraints, expected_nf)


def test_tenbars1():
    expected_x = np.array([
        5.52815824e+02,  2.66960799e+02, -5.25681197e+02,  3.64853094e+02,
        1.06566613e+02,  4.10799791e+01, -1.09350022e+02,  1.51475777e+02,
        2.07490306e+00, -1.93398099e+00, -1.93398099e+00,  2.28004667e+00,
        -1.93388575e+00, -4.78440467e-01, -1.93398099e+00, -1.93398099e+00,
        -9.03505082e-01, -1.93398099e+00])
    expected_f = -29.95537645604265
    expected_constraints = np.array([
        -317.76079933, -415.65309429,  -91.87997913, -202.27577738,
          -1.42990306,    2.57898099,    2.57898099,   -1.63504667,
           2.57888575,    1.12344047,    2.57898099,    2.57898099,
           1.54850508,    2.57898099, -120.75229495,   -0.8009505 ,
           2.57901841,   -2.57930764,    2.57919261,    2.57931513,
           2.57900977,   -2.35582445,    2.57911697,    0.8009505 ,
          -2.57901841,    2.57930764,   -2.57919261,   -2.57931513,
          -2.57900977,    2.35582445,   -2.57911697])
    expected_nf = 9000
    run_problem('TENBARS1', expected_x, expected_f, expected_constraints, expected_nf)


def test_biggs3():
    expected_x = np.array([0.99996557, 9.99841456, 4.99911898])
    expected_f = 1.5888794025529085e-08
    expected_constraints = np.array([])
    expected_nf = 1500
    run_problem('BIGGS3', expected_x, expected_f, expected_constraints, expected_nf)


def test_biggs6():
    expected_x = np.array([1.11201132, 8.8273265 , 1.06178977, 4.35354572, 3.7127316, 2.82550501])
    expected_f = 0.003956138314365143
    expected_constraints = np.array([])
    expected_nf = 3000
    run_problem('BIGGS6', expected_x, expected_f, expected_constraints, expected_nf)


def test_degenlpb():
    expected_x = np.array([
         2.50481372e-01,  8.16528208e-04,  2.80868772e-02,  9.98363594e-02,
         3.97295802e-06,  1.36025420e-04,  4.83428272e-04,  5.66993565e-04,
         1.02478447e-03,  1.96029655e-01, -2.53908006e-13, -2.53908005e-13,
        -2.53908006e-13, -9.58150253e-14, -2.53534800e-13,  1.99620738e-03,
        -2.53908006e-13,  1.65864565e-06,  3.92756445e-03,  1.20111604e-03])
    expected_f = -30.731246817982942
    expected_constraints = np.array([
           -2.50481372e-01, -8.16528208e-04, -2.80868772e-02, -9.98363594e-02,
           -3.97295802e-06, -1.36025420e-04, -4.83428272e-04, -5.66993565e-04,
           -1.02478447e-03, -1.96029655e-01,  2.53908006e-13,  2.53908005e-13,
            2.53908006e-13,  9.58150253e-14,  2.53534800e-13, -1.99620738e-03,
            2.53908006e-13, -1.65864565e-06, -3.92756445e-03, -1.20111604e-03,
           -7.49518628e-01, -9.99183472e-01, -9.71913123e-01, -9.00163641e-01,
           -9.99996027e-01, -9.99863975e-01, -9.99516572e-01, -9.99433006e-01,
           -9.98975216e-01, -8.03970345e-01, -1.00000000e+00, -1.00000000e+00,
           -1.00000000e+00, -1.00000000e+00, -1.00000000e+00, -9.98003793e-01,
           -1.00000000e+00, -9.99998341e-01, -9.96072436e-01, -9.98798884e-01,
           -2.54352095e-13, -2.87186941e-13, -2.42722509e-13, -2.53897597e-13,
           -2.75710010e-13, -2.54016643e-13, -2.53913210e-13, -2.42673937e-13,
            2.53912559e-13, -2.54027810e-13,  2.53952756e-13,  2.53908388e-13,
           -2.53907836e-13, -2.54019028e-13,  2.53914945e-13,  2.54352095e-13,
            2.87186941e-13,  2.42722509e-13,  2.53897597e-13,  2.75710010e-13,
            2.54016643e-13,  2.53913210e-13,  2.42673937e-13, -2.53912559e-13,
            2.54027810e-13, -2.53952756e-13, -2.53908388e-13,  2.53907836e-13,
            2.54019028e-13, -2.53914945e-13])
    expected_nf = 97
    run_problem('DEGENLPB', expected_x, expected_f, expected_constraints, expected_nf)
