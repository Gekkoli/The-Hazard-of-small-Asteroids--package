"""Module dealing with postcode information."""
import os
import numpy as np
import pandas as pd

__all__ = ['PostcodeLocator', 'great_circle_distance']


def great_circle_distance(latlon1, latlon2):
    """
    Calculate the great circle distance (in metres) between pairs of
    points specified as latitude and longitude on a spherical Earth
    (with radius 6371 km).

    Parameters

    ----------

    latlon1: arraylike
        latitudes and longitudes of first point (as [n, 2] array for n points)
    latlon2: arraylike
        latitudes and longitudes of second point (as [m, 2] array for m points)

    Returns
    -------

    numpy.ndarray
        Distance in metres between each pair of points (as an n x m array)

    Examples
    --------

    >>> print(great_circle_distance([[54.0, 0.0], [55, 0.0]], [55, 1.0]))
    [[128580.53670808]
     [ 63778.24657475]]

    >>> print(great_circle_distance([[54.0, 0.], [55, 0.]], \
        [[55, 1.], [54., 0.]]))
    [[128580.53670808      0.        ]
     [ 63778.24657475 111194.92664456]]
    """
    latlon1 = np.array(latlon1).reshape((-1, 1, 2)) * np.pi/180
    latlon2 = np.array(latlon2).reshape((1, -1, 2)) * np.pi/180
    # Rp * distance (Haversine formula)
    return 6371_000 * 2 * np.arcsin(np.sqrt(np.sin(
        (np.abs(latlon1[:, :, 0] - latlon2[:, :, 0]))/2)**2
        + np.cos(latlon1[:, :, 0]) * np.cos(latlon2[:, :, 0])
        * np.sin((np.abs(latlon1[:, :, 1] - latlon2[:, :, 1]))/2)**2))


class PostcodeLocator(object):
    """Class to interact with a postcode database file."""

    def __init__(self,
                 postcode_file=os.sep.join((os.path.dirname(__file__), '..',
                                            'resources',
                                            'full_postcodes.csv')),
                 census_file=os.sep.join((os.path.dirname(__file__), '..',
                                          'resources',
                                          'population_by_postcode_sector.csv')
                                         ),
                 norm=great_circle_distance):
        """
        Parameters
        ----------

        postcode_file : str, optional
            Filename of a .csv file containing geographic
            location data for postcodes.

        census_file :  str, optional
            Filename of a .csv file containing census data by postcode sector.

        norm : function
            Python function defining the distance between points in
            latitude-longitude space.

        """
        self.postcode_df = pd.read_csv(postcode_file,
                                       usecols=['Postcode',
                                                'Latitude',
                                                'Longitude'])
        self.rm_space_func = lambda postcode: ''.join(postcode.split())
        self.postcode_df['Postcode'] = self.postcode_df['Postcode'].apply(
            self.rm_space_func)
        self.census_df = pd.read_csv(census_file, usecols=[1, 4])
        self.census_df['geography'] = self.census_df['geography'].apply(
            self.rm_space_func)
        self.norm = norm
        self.postcode_filler = {5: '  ', 6: ' ', 7: ''}
        self.__add_0populations_postcode()
        self.census_df.set_index('geography', inplace=True)

    def __add_0populations_postcode(self):
        all_sector_set = set(
            self.postcode_df['Postcode'].apply(lambda x: x[:-2]).unique())
        exist_sector_set = set(self.census_df['geography'])
        remains_sectors = all_sector_set - exist_sector_set
        remains_census_df = pd.DataFrame(
            (remains_sectors, np.zeros(len(remains_sectors))),
            index=self.census_df.columns).T
        self.census_df = pd.concat((self.census_df, remains_census_df), axis=0)

    def get_postcodes_by_radius(self, X, radii, sector=False):
        """
        Return (unit or sector) postcodes within specific distances of
        input location.

        Parameters
        ----------
        X : arraylike
            Latitude-longitude pair of centre location
        radii : arraylike
            array of radial distances from X
        sector : bool, optional
            if true return postcode sectors, otherwise postcode units

        Returns
        -------
        list of lists
            Contains the lists of postcodes closer than the elements
            of radii to the location X.


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), [0.13e3])
        [['SW7 2AZ', 'SW7 2BT', 'SW7 2BU', 'SW7 2DD', 'SW7 5HF', \
'SW7 5HG', 'SW7 5HQ']]

        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), \
            [0.4e3, 0.2e3], True)
        [['SW7 1', 'SW7 2', 'SW7 3', 'SW7 4', 'SW7 5', 'SW7 9'], \
['SW7 1', 'SW7 2', 'SW7 3', 'SW7 4', 'SW7 5', 'SW7 9']]
        """
        postcodes_list = []
        self.postcode_df['Distance'] = self.norm(
            self.postcode_df[['Latitude', 'Longitude']], X)
        for radius in radii:
            postcodes = self.postcode_df[
                    self.postcode_df['Distance'] <= radius]['Postcode']
            if sector:
                postcodes = postcodes.apply(
                    lambda x:
                    f'{x[:-3]}{self.postcode_filler[len(x[:-2])+2]}{x[-3]}'
                    ).unique().tolist()
            else:
                postcodes = postcodes.apply(
                    lambda x:
                    f'{x[:-3]}{self.postcode_filler[len(x)]}{x[-3:]}'
                    ).to_list()
            postcodes_list.append(postcodes)
        return postcodes_list

    def get_population_of_postcode(self, postcodes, sector=False):
        """
        Return populations of a list of postcode units or sectors.

        Parameters
        ----------
        postcodes : list of lists
            list of postcode units or postcode sectors
        sector : bool, optional
            if true return populations for postcode sectors,
            otherwise returns populations for postcode units

        Returns
        -------
        list of lists
            Contains the populations of input postcode units or sectors


        Examples
        --------

        >>> locator = PostcodeLocator()
        >>> locator.get_population_of_postcode([['SW7 2AZ', 'SW7 2BT', \
            'SW7 2BU', 'SW7 2DD']])
        [[19.0, 19.0, 19.0, 19.0]]

        >>> locator.get_population_of_postcode([['SW7  2']], True)
        [[2283.0]]
        """
        populations = []
        if sector:
            for codes in postcodes:
                if len(codes) == 0:
                    populations.append(0)
                    continue
                codes = pd.Series(codes).apply(self.rm_space_func)
                sector_populations = np.array(
                    self.census_df.loc[codes].values.ravel(), dtype=float)
                populations.append(np.round(sector_populations).tolist())
            return populations

        units_count = self.postcode_df['Postcode'].apply(
            lambda x: x[:-2]).value_counts()
        for codes in postcodes:
            if len(codes) == 0:
                populations.append(0)
                continue
            codes = pd.Series(codes).apply(lambda x: ''.join(x.split())[:-2])
            unit_populations = self.census_df.loc[codes].values.ravel()
            unit_populations = np.array(
                unit_populations, dtype=float) / units_count[codes]
            populations.append(np.round(unit_populations).tolist())
        return populations
