import armageddon

#######################
#   Airburst Solver   #
#######################

# Initialise the Planet class
earth = armageddon.Planet()

# Solve the atmospheric entry problem for a given set of input parameters
result = earth.solve_atmospheric_entry(radius=35, angle=45,
                                       strength=1e7, density=3000,
                                       velocity=19e3)

# Calculate the kinetic energy lost per unit altitude and add it
# as a column to the result dataframe
result = earth.calculate_energy(result)

# Determine the outcomes of the impact event
outcome = earth.analyse_outcome(result)

#####################
#   Damage Mapper   #
#####################

# Calculate the blast location and damage radius for several pressure levels
pressures = [1e3, 3.5e3, 27e3, 43e3]
blast_lat, blast_lon, damage_rad = armageddon.damage_zones(outcome,
                                                           lat=53, lon=-2.5,
                                                           bearing=115.,
                                                           pressures=pressures)

# Plot a circle to show the limit of the lowest damage level
damage_map = armageddon.plot_circle(blast_lat, blast_lon, damage_rad[0])
damage_map.save("damage_map.html")

# The PostcodeLocator tool
locator = armageddon.PostcodeLocator()

# Find the postcodes in the damage radii
postcodes = locator.get_postcodes_by_radius((blast_lat, blast_lon),
                                            radii=damage_rad)

# Find the population in each postcode
population = locator.get_population_of_postcode(postcodes)

# Print the number of people affected in each damage zone
for pop, zone in zip(population, pressures):
    print("{:5.0f} Pa: {:8.0f}".format(zone, sum(pop)))

# Alternatively find the postcode sectors in the damage radii,
# and populations of the sectors
sectors = locator.get_postcodes_by_radius((blast_lat, blast_lon),
                                          radii=damage_rad, sector=True)
population_sector = locator.get_population_of_postcode(sectors, sector=True)
