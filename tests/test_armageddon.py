import os
import pandas as pd
import numpy as np

from collections import OrderedDict
from pytest import fixture

from armageddon.solver import SolveODE
from armageddon.locator import PostcodeLocator

pl = PostcodeLocator()


# Use pytest fixtures to generate objects we know we'll reuse.
# This makes sure tests run quickly

@fixture(scope='module')
def armageddon():
    import armageddon
    return armageddon


@fixture(scope='module')
def planet(armageddon):
    return armageddon.Planet()


@fixture(scope='module')
def loc(armageddon):
    return armageddon.PostcodeLocator()


@fixture(scope='module')
def result(planet):
    input = {'radius': 1.,
             'velocity': 2.0e4,
             'density': 3000.,
             'strength': 1e5,
             'angle': 30.0,
             'init_altitude': 1e3
             }

    result = planet.solve_atmospheric_entry(**input)

    return result


@fixture(scope='module')
def energy(planet, result):
    outcome = planet.calculate_energy(result=result)
    return outcome


@fixture(scope='module')
def outcome(planet, energy):
    outcome = planet.analyse_outcome(result=energy)
    return outcome


def test_import(armageddon):
    assert armageddon


def test_planet_signature(armageddon):
    inputs = OrderedDict(atmos_func='exponential',
                         atmos_filename=None,
                         Cd=1., Ch=0.1, Q=1e7, Cl=1e-3,
                         alpha=0.3, Rp=6371e3,
                         g=9.81, H=8000., rho0=1.2)

    # call by keyword
    planet = armageddon.Planet(**inputs)  # noqa

    # call by position
    planet = armageddon.Planet(*inputs.values())  # noqa


def test_attributes(planet):
    for key in ('Cd', 'Ch', 'Q', 'Cl',
                'alpha', 'Rp', 'g', 'H', 'rho0'):
        assert hasattr(planet, key)


def test_atmos_filename(planet):

    assert os.path.isfile(planet.atmos_filename)


def test_solve_atmospheric_entry(result):

    assert type(result) is pd.DataFrame

    for key in ('velocity', 'mass', 'angle', 'altitude',
                'distance', 'radius', 'time'):
        assert key in result.columns


def test_simple_model():
    Cd = 1.
    H = 8000.
    rho0 = 1.2
    density = 3000
    velocity = 19000
    radius = 35
    altitude = 100000
    mass = (4/3) * np.pi * radius**3 * density
    theta = 45*np.pi/180
    distance = 0.0
    A = np.pi * radius**2

    # Simplified model RHS
    def simple_rhs(t, y):
        dy = np.zeros_like(y)
        dy[0] = -Cd*A*rho0*0.5*np.exp(-y[1]/H)*y[0]**2/mass
        dy[1] = -y[0]*np.sin(theta)
        dy[2] = y[0]*np.cos(theta)
        return dy

    # Solve the model as requested
    y0 = [velocity, altitude, distance]
    solver = SolveODE(simple_rhs, (0.0, 50.0), y0, tol=1e-8, maxstep=0.2)
    solver.solve()
    sol = solver.sol
    vel = sol[:, 0]
    z = sol[:, 1]

    lhs = np.log(vel / velocity)
    rhs = (-H * Cd * A * rho0 / (2 * mass * np.sin(theta))) * \
          (np.exp(-z / H) - np.exp(-altitude / H))

    assert np.isclose(lhs, rhs).all()


def test_rk4_solver():
    def func(x, y):
        return np.exp(x)
    solver = SolveODE(func, (0, 1), 1, maxstep=0.05, method='RK4')
    solver.solve()

    assert np.allclose(solver.sol, np.exp(solver.t))


def test_rk45_solver():
    def func(x, y):
        res = np.exp(x)
        return res

    solver = SolveODE(func, (0, 1), 1, maxstep=0.05, method='RK45')
    solver.solve()

    assert np.allclose(solver.sol, np.exp(solver.t))


def test_interplotion(armageddon):
    h = [2100, 30000, 41000, 71000, 81000]
    planet = armageddon.Planet(atmos_func='tabular')
    rho = planet.rhoa(h)
    res = np.array([9.96407070e-01, 1.80119000e-02, 3.32643241e-03,
                    6.41309249e-05, 1.33195940e-05])
    assert np.allclose(rho, res)


def test_calculate_energy(planet, result):

    energy = planet.calculate_energy(result=result)

    assert type(energy) is pd.DataFrame

    for key in ('velocity', 'mass', 'angle', 'altitude',
                'distance', 'radius', 'time', 'dedz'):
        assert key in energy.columns


def test_analyse_outcome(planet, outcome):

    assert type(outcome) is dict

    for key in ('outcome', 'burst_peak_dedz', 'burst_altitude',
                'burst_distance', 'burst_energy'):
        assert key in outcome.keys()


def test_airburst_case(planet):
    input = {'radius': 9.75,
             'velocity': 19e3,
             'density': 3300.,
             'strength': 2e6,
             'angle': 20.,
             'init_altitude': 100e3}

    result = planet.solve_atmospheric_entry(**input)
    result = planet.calculate_energy(result)
    outcome = planet.analyse_outcome(result)

    ans = [91.49271571624418, 31425.334946666408,
           194027.08322546008, 361.10530125690303]

    assert outcome["outcome"] == "Airburst"
    out = [outcome["burst_peak_dedz"], outcome["burst_altitude"],
           outcome["burst_distance"], outcome["burst_energy"]]

    assert np.allclose(out, ans)


def test_cratering_case(planet):
    input = {'radius': 50,
             'velocity': 20e3,
             'density': 7000.,
             'strength': 2e7,
             'angle': 60,
             'init_altitude': 100e3}

    result = planet.solve_atmospheric_entry(**input)
    result = planet.calculate_energy(result)
    outcome = planet.analyse_outcome(result)

    ans = [19663.838902840696, 0.0,
           57792.93316160227, 116758.92480298487]

    assert outcome["outcome"] == "Cratering"

    out = [outcome["burst_peak_dedz"], outcome["burst_altitude"],
           outcome["burst_distance"], outcome["burst_energy"]]

    assert np.allclose(out, ans)


def test_damage_zones(armageddon):

    outcome = {'burst_peak_dedz': 1000.,
               'burst_altitude': 9000.,
               'burst_distance': 90000.,
               'burst_energy': 6000.,
               'outcome': 'Airburst'}

    blat, blon, damrad = armageddon.damage_zones(outcome, 55.0, 0.,
                                                 135., [27e3, 43e3])

    assert type(blat) is float
    assert type(blon) is float
    assert type(damrad) is list
    assert len(damrad) == 2


def test_great_circle_distance(armageddon):

    pnts1 = np.array([[54.0, 0.0], [55.0, 1.0], [54.2, -3.0]])
    pnts2 = np.array([[55.0, 1.0], [56.0, -2.1], [54.001, -0.003]])
    pnts4 = np.array([[55.0, 1.0], [56.0, -2.1]])

    data = np.array([[1.28580537e+05, 2.59579735e+05, 2.25409117e+02],
                    [0.00000000e+00, 2.24656571e+05, 1.28581437e+05],
                    [2.72529953e+05, 2.08175028e+05, 1.96640630e+05]])

    data2 = [[0.]]
    data3 = np.array([1.286e+05, 6.378e+04])
    data5 = (pnts1.shape[0], pnts4.shape[0])

    dist = armageddon.great_circle_distance(pnts1, pnts2)
    dist2 = armageddon.great_circle_distance([55.0, 1.0], [55.0, 1.0])
    dist3 = armageddon.great_circle_distance(
        [[54.0, 0.0], [55, 0.0]], [55, 1.0])
    dist4 = type(armageddon.great_circle_distance(
        np.array([[54.0, 0.0], [55.0, 1.0], [54.2, -3.0]]),
        np.array([[55.0, 1.0], [56.0, -2.1], [54.001, -0.003]])))
    dist5 = armageddon.great_circle_distance(pnts1, pnts4).shape

    assert np.allclose(data, dist, rtol=1.0e-4)  # Simple test from docstring
    # check that distance is really zero
    assert np.isclose(data2, dist2, rtol=1.0e-4)
    # have to do this otherwise will be too precise
    assert np.allclose(data3, dist3, rtol=1.0e1)
    assert (dist4 == np.ndarray)  # Check that the ouput is a numpy array
    assert (dist5 == data5)  # Check that the ouput array is of size n x m


def test_get_postcodes_by_radius(loc):
    # tests on corner cases
    assert loc.get_postcodes_by_radius((0, 0), [0], True) == [[]]
    assert loc.get_postcodes_by_radius((0, 0), [0]) == [[]]

    postcodes = loc.get_postcodes_by_radius(
        (51.4981, -0.1773), [150, 300, 10000])
    assert len(postcodes[0]) == 9
    assert len(postcodes[1]) == 94
    assert len(postcodes[2]) == 76385


def test_get_population_by_postcodes(loc):

    data1 = [[19., 19., 19., 19.]]
    data2 = [[2283]]
    data3 = np.array([[27., 19.], [43., 38.]])
    data4 = list
    data5 = [[2283.0], [9856.0]]

    # Checking no space in the Input
    Answer1 = loc.get_population_of_postcode([['SW72AZ', 'SW72BT',
                                               'SW72BU', 'SW72DD']])
    # Checking Space in Input
    Answer4 = loc.get_population_of_postcode([['SW7 2AZ', 'SW7 2BT',
                                               'SW7 2BU', 'SW7 2DD']])
    # Checking PostCode Sector with Space
    Answer2 = loc.get_population_of_postcode([['SW7 2']], True)
    # Checking list of list Postcode Sector with Space
    Answer5 = loc.get_population_of_postcode([['SW7 2'], ['HD7 5']], True)
    # Checking list of list Postcode Sector with no-Space
    Answer8 = loc.get_population_of_postcode([['SW72'], ['HD75']], True)
    # Checking List of Lists with no-Space
    Answer3 = loc.get_population_of_postcode([['EX172BS', 'SW72BU'],
                                             ['RM154DU', 'HD75PG']])
    # Checking List of Lists with Space
    Answer6 = loc.get_population_of_postcode([['EX17 2BS', 'SW7 2BU'],
                                              ['RM15 4DU', 'HD7 5PG']])

    Answer7 = type(loc.get_population_of_postcode([['EX172BS'], ['HD75PG']]))

    # tests on corner cases
    assert loc.get_population_of_postcode([[]], True) == [0]
    assert loc.get_population_of_postcode([['SW75']], True) == [[5058.0]]
    assert loc.get_population_of_postcode([[]]) == [0]

    # Check implementation with no space in the postcode
    assert np.allclose(data1, Answer1)
    # Check that if the input includes a space, it still works
    assert np.allclose(data1, Answer4)

    # Check that the Boolean application works with PC Sectors
    assert np.allclose(data2, Answer2)
    # Check for that it also works with space inside the PC string
    assert np.allclose(data5, Answer5)
    assert np.allclose(data5, Answer8)

    # Check that implementation of List of list is correct
    assert np.allclose(data3, Answer3)
    assert np.allclose(data3, Answer6)  # Check it also works with spacing

    assert (data4 == Answer7)  # Assert that output is a list of lists


def test_locator_postcodes(loc):

    latlon = (52.2074, 0.1170)

    result = loc.get_postcodes_by_radius(latlon, [0.2e3, 0.1e3])
    result2 = loc.get_postcodes_by_radius(latlon, [9e1, 2e2])

    data = [['CB2 1TP', 'CB2 1TQ'],
            ['CB2 1SU', 'CB2 1TA',
             'CB2 1TB', 'CB2 1TJ',
             'CB2 1TP', 'CB2 1TQ',
             'CB2 1TR', 'CB2 1TW',
             'CB2 1TY', 'CB2 1UA',
             'CB2 1UF', 'CB2 1UR',
             'CB2 3HY', 'CB2 3JU',
             'CB2 3JX', 'CB2 3LL',
             'CB2 3LR', 'CB2 3LS',
             'CB5 8AD']]

    assert type(result) is list  # Checking if the output is a list of lists
    if len(result) > 0:
        for element in result:
            assert type(element) is list

    assert (result2 == data)
    # One more test to see if the size of the result is of 7:
    for element in result2:
        assert len(element[0]) == 7


def test_locator_sectors(loc):

    latlon = (52.2074, 0.1170)

    # Checking that output is Sectors
    result = loc.get_postcodes_by_radius(latlon, [3.0e3, 1.5e3], True)

    data_sectors = [['CB1 0', 'CB1 1', 'CB1 2', 'CB1 3',
                     'CB1 7', 'CB1 8', 'CB1 9', 'CB2 0',
                     'CB2 1', 'CB2 3', 'CB2 7', 'CB2 8',
                     'CB215', 'CB223', 'CB236', 'CB237',
                     'CB243', 'CB249', 'CB3 0', 'CB3 1',
                     'CB3 9', 'CB4 1', 'CB4 2', 'CB4 3',
                     'CB5 8'],
                    ['CB1 0', 'CB1 1', 'CB1 2', 'CB2 1',
                     'CB2 3', 'CB2 7', 'CB3 0', 'CB3 9',
                     'CB4 1', 'CB4 2', 'CB4 3', 'CB5 8']]

    assert type(result) is list  # Checking if the output is a list of lists
    if len(result) > 0:
        for element in result:
            assert type(element) is list

    # One more test to see if the size of the result is of 5:
    for element in result:
        assert len(element[0]) == 5

    assert (result == data_sectors)
