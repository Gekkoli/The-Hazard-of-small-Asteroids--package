import folium
from folium.plugins import HeatMap
from armageddon.locator import PostcodeLocator

__all__ = ['plot_circle', 'plot_impact_risk_circles']


def plot_circle(lat, lon, radii, map_obj=None, **kwargs):
    """
    Plot circles on a map (creating a new folium map instance if necessary).

    Parameters
    ----------

    lat: float
        latitude of circle to plot (degrees)
    lon: float
        longitude of circle to plot (degrees)
    radii: arraylike
        array of radial distances of circle to plot (m)
    map_obj: folium.Map
        existing map object

    Returns
    -------

    Folium map object

    Examples
    --------

    >>> type(plot_circle(52.79, -2.95, [1e3], map_obj=None))
    <class 'folium.folium.Map'>
    """

    if map_obj is None:
        map_obj = folium.Map(
            location=[lat, lon], control_scale=True, zoom_start=8)

    locator = PostcodeLocator()
    locator.postcode_df['Populations'] = locator.get_population_of_postcode(
        [locator.postcode_df['Postcode']])[0]
    lat_lon_populations_df = locator.postcode_df.drop(columns=['Postcode'])

    HeatMap(lat_lon_populations_df).add_to(map_obj)

    folium.Marker(
        location=[lat, lon],
        icon=None,
        tooltip='({:.3f}, {:.3f})'.format(lat, lon)
    ).add_to(map_obj)

    colors = ['purple', 'red', 'DarkOrange', 'yellow']
    for r, c in zip(radii, colors):
        folium.Circle([lat, lon], r, tooltip="radius={:.3f}".format(r),
                      color=c, fill=True, fill_color=c,
                      fillOpacity=0.7, **kwargs).add_to(map_obj)

    return map_obj


def plot_impact_risk_circles(lat_lon_radii_list, map_obj=None, **kwargs):
    """
    Plot a circle on a map (creating a new folium map instance if necessary).

    Parameters
    ----------

    lat_lon_radii_list: list of lists
        latitude-longitude-radius of circle to plot (degrees, degrees, m)
    map_obj: folium.Map
        existing map object

    Returns
    -------

    Folium map object

    Examples
    --------

    >>> type(plot_impact_risk_circles([(52.79, -2.95, 1e3)], map_obj=None))
    <class 'folium.folium.Map'>
    """

    if map_obj is None:
        lat, lon = lat_lon_radii_list[0][:-1]
        map_obj = folium.Map(
            location=[lat, lon], control_scale=True, zoom_start=8)

    for lat, lon, r in lat_lon_radii_list:
        folium.Circle([lat, lon], r, tooltip="radius={:.3f}".format(r),
                      color='blue', fill=True, fill_color='lightblue',
                      fillOpacity=0.99, **kwargs).add_to(map_obj)

    return map_obj
