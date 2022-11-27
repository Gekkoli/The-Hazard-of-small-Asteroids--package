import pandas as pd
import numpy as np
from scipy.optimize import root
from armageddon.locator import PostcodeLocator
from collections import Counter
from armageddon.mapping import plot_impact_risk_circles

__all__ = ['damage_zones', 'impact_risk']


def damage_zones(outcome, lat, lon, bearing, pressures):
    """
    Calculate the latitude and longitude of the surface zero location and the
    list of airblast damage radii (m) for a given impact scenario.

    Parameters
    ----------

    outcome: Dict
        the outcome dictionary from an impact scenario
    lat: float
        latitude of the meteoroid entry point (degrees)
    lon: float
        longitude of the meteoroid entry point (degrees)
    bearing: float
        Bearing (azimuth) relative to north of meteoroid trajectory (degrees)
    pressures: float, arraylike
        List of threshold pressures to define airblast damage levels

    Returns
    -------

    blat: float
        latitude of the surface zero point (degrees)
    blon: float
        longitude of the surface zero point (degrees)
    damrad: arraylike, float
        List of distances specifying the blast radii
        for the input damage levels

    Examples
    --------

    >>> import armageddon
    >>> import numpy
    >>> outcome = {'burst_altitude': 8e3, 'burst_energy': 7e3, \
                   'burst_distance': 90e3, 'burst_peak_dedz': 1e3, \
                   'outcome': 'Airburst'}
    >>> blat, blon, damrad = armageddon.damage_zones(outcome, 52.79, \
        -2.95, 135, pressures=[1e3, 3.5e3, 27e3, 43e3])
    >>> numpy.round(blat, 3), numpy.round(blon, 3), numpy.round(damrad, 3)
    (52.214, -2.016, array([115971.317,  42628.367,   9575.214,   5835.983]))
    """

    earth_radius = 6371000
    lat, bearing = np.deg2rad(lat), np.deg2rad(bearing)
    blat = np.arcsin(np.sin(lat) *
                     np.cos(outcome['burst_distance']/earth_radius) +
                     np.cos(lat) *
                     np.sin(outcome['burst_distance']/earth_radius) *
                     np.cos(bearing))

    blon = np.arctan((np.sin(bearing) *
                      np.sin(outcome['burst_distance']/earth_radius) *
                      np.cos(lat)) /
                     (np.cos(outcome['burst_distance']/earth_radius) -
                      np.sin(lat)*np.sin(blat)))
    blat = np.rad2deg(blat)
    blon = np.rad2deg(blon) + lon

    damrad = solve(outcome['burst_energy'],
                   outcome['burst_altitude'], pressures)

    return float(blat), float(blon), damrad


def solve(energy, altitude, pressures):
    radius = []
    for pressure in pressures:
        def f(r, E=energy, z=altitude, p=pressure):
            return 3.14*(10**11)*((r**2+z**2)/E**(2/3))**(-1.3) + \
                1.8*(10**7)*((r**2+z**2)/E**(2/3))**(-0.565) - p
        a = root(f, 120e3)
        radius.append(a.x[0])
    return radius


fiducial_means = {'radius': 35, 'angle': 45, 'strength': 1e7,
                  'density': 3000, 'velocity': 19e3,
                  'lat': 53.0, 'lon': -2.5, 'bearing': 115.}
fiducial_stdevs = {'radius': 1, 'angle': 1, 'strength': 5e6,
                   'density': 500, 'velocity': 1e3,
                   'lat': 0.025, 'lon': 0.025, 'bearing': 0.5}


def impact_risk(planet, means=fiducial_means, stdevs=fiducial_stdevs,
                pressure=27.e3, nsamples=100, sector=True):
    """
    Perform an uncertainty analysis to calculate the risk for each affected
    UK postcode or postcode sector

    Parameters
    ----------
    planet: armageddon.Planet instance
        The Planet instance from which to solve the atmospheric entry

    means: dict
        A dictionary of mean input values for the uncertainty analysis. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    stdevs: dict
        A dictionary of standard deviations for each input value. This
        should include values for ``radius``, ``angle``, ``strength``,
        ``density``, ``velocity``, ``lat``, ``lon`` and ``bearing``

    pressure: float
        A single pressure at which to calculate the damage zone for each impact

    nsamples: int
        The number of iterations to perform in the uncertainty analysis

    sector: logical, optional
        If True (default) calculate the risk for postcode sectors, otherwise
        calculate the risk for postcodes

    Returns
    -------
    risk: DataFrame
        A pandas DataFrame with columns for postcode (or postcode sector) and
        the associated risk. These should be called ``postcode`` or ``sector``,
        and ``risk``.
    """

    pl = PostcodeLocator()
    all_keys = list(means.keys())
    matrix = []
    for key in all_keys:
        arr = np.random.normal(means[key],
                               stdevs[key],
                               nsamples)
        np.random.shuffle(arr)
        matrix.append(arr.tolist())

    matrix = np.array(matrix).T
    nsamples_postcodes = []
    lat_lon_radii_list = []
    for r, a, s, d, v, lat, lon, bearing in matrix:
        result_df = planet.solve_atmospheric_entry(r, v, d, s, a)
        cal_enrgy_df = planet.calculate_energy(result_df)
        outcome = planet.analyse_outcome(cal_enrgy_df)

        blat, blon, diaran = damage_zones(outcome, lat,
                                          lon, bearing, [pressure])
        centre = (blat, blon)
        if (radius := diaran[0]) > 0:
            lat_lon_radii_list.append((blat, blon, radius))
        nsamples_postcodes.extend(
            pl.get_postcodes_by_radius(centre, diaran, sector)[0])

    postcode_count_dict = Counter(nsamples_postcodes)
    postcodes = list(postcode_count_dict.keys())
    counts = np.array(list(postcode_count_dict.values()))
    postcode_probability = counts / nsamples

    if not postcodes:
        risk = []
    else:
        postcode_population = np.array(
            [pl.get_population_of_postcode
             ([[x]], sector)[0][0] for x in postcodes])
        risk = [(prob*popu) for (prob, popu)
                in zip(postcode_probability, postcode_population)]

    risk_map = plot_impact_risk_circles(lat_lon_radii_list)
    risk_map.save('risk_map.html')

    if sector:
        return pd.DataFrame({'sector': postcodes, 'risk': risk}).\
            sort_values(by='risk', ascending=False)
    else:
        return pd.DataFrame({'postcode': postcodes, 'risk': risk}).\
            sort_values(by='risk', ascending=False)
