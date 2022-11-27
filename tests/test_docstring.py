import doctest
from armageddon import solver, damage, locator, mapping


def test_solver_docstrings():
    assert doctest.testmod(solver).failed == 0, \
        'Failed docstring tests in armageddon.solver :('


def test_damage_docstrings():
    assert doctest.testmod(damage).failed == 0, \
        'Failed docstring tests in armageddon.damage :('


def test_locator_docstrings():
    assert doctest.testmod(locator).failed == 0, \
        'Failed docstring tests in armageddon.locator :('


def test_mapping_docstrings():
    assert doctest.testmod(mapping).failed == 0, \
        'Failed docstring tests in armageddon.mapping :('
