import os

import numpy as np
import scipy.interpolate as si

import pandas as pd

__all__ = ['Planet']


class SolveODE:

    def __init__(self, func, tspan, y0, args=None,
                 method='RK45', tol=1e-6, stopping=None,
                 initstep=0.05, maxstep=0.1):
        """
        Constructor for the SolveODE class.

        Parameters
        ----------
        func : callable(t, y, ...)
            Return the RHS of the ODE system
        tspan : tuple (float, float)
            The interval of integration
        y0 : array
            Initial state vector
        args: tuple
            Arguments for the RHS function
        method: string
            To select the solver. Options 'Euler', 'RK4', 'RK45'
        tol : float
            Relative Error tolerance for adaptive stepping
        stopping: callback function list
            Functions which return true when stopping condition reached
        initstep: float
            step size to start integration. Remains same in Euler and RK4,
            Changes in RK45
        maxstep : float
            max allowed step size for adaptive stepping in RK45

        Attributes
        ----------
        h : float
            Current time step
        t : list
            list of solver times
        sol : list
            list of solution at corresponding times
        stopping : list
            list containing stopping conditions

        Methods
        -------
        solve() : Executes the solver as mentioned in method and stores output

        Examples
        --------
        >>> import numpy as np
        >>> from armageddon.solver import SolveODE
        >>> def func(t, y):
        ...     return np.cos(t)
        ...
        >>> solver = SolveODE(func, (0, np.pi/4), 0.0)
        >>> solver.solve()
        >>> np.round(solver.t, decimals=5)
        array([0.     , 0.02656, 0.05766, 0.15766, 0.25766, 0.35766, 0.45766,
               0.55766, 0.65766, 0.75766, 0.7854 ])
        >>> np.round(solver.sol, decimals=5)
        array([0.     , 0.02655, 0.05763, 0.15701, 0.25482, 0.35008, 0.44185,
               0.5292 , 0.61127, 0.68722, 0.70711])

        """
        self.rhs = func
        self.tspan = tspan
        self.y0 = y0
        if args is None:
            self.args = ()
        else:
            self.args = args
        self.method = method
        self.h = initstep
        self.hmax = maxstep
        self.tol = tol
        self.stop = stopping
        self.t = []
        self.sol = []
        # Only for analysis. Comment out after that
        # self.err = []

    # Euler Integration
    def __euler(self):
        t0 = self.tspan[0]
        tf = self.tspan[1]
        self.t.append(t0)
        self.sol.append(self.y0)
        if self.h >= (tf - t0):
            return
        tk = t0
        yk = self.y0

        while tk < tf:
            yk = yk + self.rhs(tk, yk, *self.args)*self.h
            tk = tk + self.h

            self.t.append(tk)
            self.sol.append(yk)

            # Checking for stopping conditions
            if self.stop is not None:
                stopCondition = np.array([stopcriterion(tk, yk)
                                          for stopcriterion in self.stop])
                if stopCondition.any():
                    return

    # RK4 implementation
    def __rk4(self):

        t0 = self.tspan[0]
        tf = self.tspan[1]
        self.tol = self.tol/(tf - t0)
        self.t.append(t0)
        self.sol.append(self.y0)
        if self.h >= (tf - t0):
            return
        tk = t0
        yk = self.y0

        while tk < tf:
            k1 = self.rhs(tk, yk, *self.args) * self.h
            k2 = self.rhs(tk + self.h/2.0,
                          yk + k1/2.0, *self.args) * self.h
            k3 = self.rhs(tk + self.h/2.0,
                          yk + k2/2.0, *self.args) * self.h
            k4 = self.rhs(tk + self.h,
                          yk + k3, *self.args) * self.h

            yk = yk + (k1 + 2*k2 + 2*k3 + k4)/6
            tk = tk + self.h

            self.t.append(tk)
            self.sol.append(yk)

            # Checking for stopping conditions
            if self.stop is not None:
                stopCondition = np.array([stopcriterion(tk, yk)
                                          for stopcriterion in self.stop])
                if stopCondition.any():
                    return

    # RK45 implementation
    def __rk45(self):

        """
        Solve the system of ODEs using RK45 Method

        Called internally by the solve method and not exposed to the user
        Populates time and y using adaptive step size

        Reference
        ---------
        Implementation as described in
        Meysam Mahooti (2022). Runge-Kutta-Fehlberg (RKF45)
        """

        # Constants
        err_pow = 0.25
        min_delta = 0.124
        max_delta = 4.0

        t0 = self.tspan[0]
        tf = self.tspan[1]
        self.t.append(t0)
        self.sol.append(self.y0)
        if self.h >= (tf - t0):
            return

        tk = t0
        yk = self.y0

        # To count unsuccesful step size updates
        tries = 0

        while tk < tf:

            k1 = self.rhs(tk, yk, *self.args) * self.h
            k2 = self.rhs(tk + self.h/4.0,
                          yk + k1/4.0, *self.args) * self.h
            k3 = self.rhs(tk + 3.0*self.h/8.0,
                          yk + (3.0*k1/32.0 +
                                9.0*k2/32.0), *self.args) * self.h
            k4 = self.rhs(tk + 12.0*self.h/13.0,
                          yk + (1932.0*k1/2197.0 - 7200.0*k2/2197.0 +
                                7296.0*k3/2197.0), *self.args) * self.h
            k5 = self.rhs(tk + self.h,
                          yk + (439.0*k1/216.0 - 8.0*k2 +
                                3680.0*k3/513.0 - 845.0*k4/4104.0),
                          *self.args) * self.h
            k6 = self.rhs(tk + self.h/2.0,
                          yk + (-8.0*k1/27.0 + 2.0*k2 -
                                3544.0*k3/2565.0 + 1859.0*k4/4104.0 -
                                11.0*k5/40.0), *self.args) * self.h

            # Error between 5th and 4th order
            r = (k1/360 - 128*k3/4275 - 2197*k4/75240 +
                 k5/50 + 2*k6/55)

            normy = np.linalg.norm(yk)
            if normy == 0:
                normy = self.tol

            err = np.linalg.norm(r)/normy

            if err == 0:
                delta = max_delta
            else:
                delta = 0.84*(self.tol/err)**err_pow

            if (err <= self.tol):
                tries = 0
                tk = tk + self.h
                yk = yk + (25.0*k1/216.0 + 1408.0*k3/2565.0 +
                           2197.0*k4/4104.0 - 1.0*k5/5.0)

                self.t.append(tk)
                self.sol.append(yk)
                # self.err.append(err)

                # Checking for stopping conditions
                if self.stop is not None:
                    stopCondition = np.array([stopcriterion(tk, yk)
                                              for stopcriterion in self.stop])
                    if stopCondition.any():
                        return

            else:
                tries += 1

            if delta <= min_delta:
                delta = min_delta
            elif delta >= max_delta:
                delta = max_delta

            self.h = delta * self.h

            if self.h > self.hmax:
                self.h = self.hmax

            if (tk != tf) and (tk+self.h) > tf:
                self.h = tf-tk

            if tries > 20:
                print("Relax Tolerances")
                return

    def solve(self):
        if self.method == 'Euler':
            self.__euler()
        if self.method == 'RK4':
            self.__rk4()
        if self.method == 'RK45':
            self.__rk45()
        self.t = np.array(self.t)
        self.sol = np.array(self.sol)


class Planet():
    """
    The class called Planet is initialised with constants appropriate
    for the given target planet, including the atmospheric density profile
    and other constants
    """

    def __init__(self, atmos_func='exponential',
                 atmos_filename=os.sep.join((os.path.dirname(__file__), '..',
                                             'resources',
                                             'AltitudeDensityTable.csv')),
                 Cd=1., Ch=0.1, Q=1e7, Cl=1e-3, alpha=0.3,
                 Rp=6371e3, g=9.81, H=8000., rho0=1.2):
        """
        Set up the initial parameters and constants for the target planet

        Parameters
        ----------
        atmos_func : string, optional
            Function which computes atmospheric density, rho, at altitude, z.
            Default is the exponential function rho = rho0 exp(-z/H).
            Options are 'exponential', 'tabular' and 'constant'

        atmos_filename : string, optional
            Name of the filename to use with the tabular atmos_func option

        Cd : float, optional
            The drag coefficient

        Ch : float, optional
            The heat transfer coefficient

        Q : float, optional
            The heat of ablation (J/kg)

        Cl : float, optional
            Lift coefficient

        alpha : float, optional
            Dispersion coefficient

        Rp : float, optional
            Planet radius (m)

        rho0 : float, optional
            Air density at zero altitude (kg/m^3)

        g : float, optional
            Surface gravity (m/s^2)

        H : float, optional
            Atmospheric scale height (m)

        """

        # Input constants
        self.Cd = Cd
        self.Ch = Ch
        self.Q = Q
        self.Cl = Cl
        self.alpha = alpha
        self.Rp = Rp
        self.g = g
        self.H = H
        self.rho0 = rho0
        self.atmos_filename = atmos_filename
        self.atmos_func = atmos_func

        try:
            # set function to define atmoshperic density
            if atmos_func == 'exponential':
                self.rhoa = lambda x: self.rho0*np.exp(-x/self.H)
            elif atmos_func == 'tabular':
                self.densityTable = np.genfromtxt(self.atmos_filename)
                Plinear1 = si.interp1d(self.densityTable[:, 0],
                                       self.densityTable[:, 1],
                                       'cubic',
                                       fill_value=(self.densityTable[0, 1],
                                                   self.densityTable[-1, 1]),
                                       bounds_error=False)
                self.rhoa = lambda x: Plinear1(x)
            elif atmos_func == 'constant':
                self.rhoa = lambda x: rho0
            else:
                raise NotImplementedError(
                    "atmos_func must be 'exponential', 'tabular' or 'constant'"
                    )
        except NotImplementedError:
            print("atmos_func {} not implemented yet.".format(atmos_func))
            print("Falling back to constant density atmosphere for now")
            self.rhoa = lambda x: rho0

    def solve_atmospheric_entry(
            self, radius, velocity, density, strength, angle,
            init_altitude=100e3, dt=0.05, radians=False, method='RK45'):
        """
        Solve the system of differential equations for a given impact scenario

        Parameters
        ----------
        radius : float
            The radius of the asteroid in meters

        velocity : float
            The entery speed of the asteroid in meters/second

        density : float
            The density of the asteroid in kg/m^3

        strength : float
            The strength of the asteroid (i.e. the maximum pressure it can
            take before fragmenting) in N/m^2

        angle : float
            The initial trajectory angle of the asteroid to the horizontal
            By default, input is in degrees. If 'radians' is set to True, the
            input should be in radians

        init_altitude : float, optional
            Initial altitude in m

        dt : float, optional
            The output timestep, in s

        radians : logical, optional
            Whether angles should be given in degrees or radians. Default=False
            Angles returned in the dataframe will have the same units as the
            input

        Returns
        -------
        Result : DataFrame
            A pandas dataframe containing the solution to the system.
            Includes the following columns:
            'velocity', 'mass', 'angle', 'altitude',
            'distance', 'radius', 'time'
        """
        def pallas_rhs(t, y):

            """ Fuction to return the RHS of the asteroid equations
            Paramters
            ---------
            t : scalar time
                Current time during integration

            y : np.array
                Statevector array

            Note
            ----
            y[0] = v        Velocity of the asteroid
            y[1] = m        Mass of the asteroid
            y[2] = theta    Angle in the entry plane
            y[3] = z        Altitude
            y[4] = x        Range
            y[5] = r        Radius of the asteroid
            """

            dy = np.zeros_like(y)

            rhoa = self.rhoa(y[3])
            A = np.pi*y[5]**2

            dy[0] = -0.5*self.Cd*rhoa*A*y[0]**2/y[1] + self.g*np.sin(y[2])
            dy[1] = -0.5*self.Ch*rhoa*A*y[0]**3/self.Q
            dy[2] = self.g*np.cos(y[2])/y[0] -\
                0.5*self.Cl*rhoa*A*y[0]/y[1] - y[0]*np.cos(y[2])/(self.Rp+y[3])
            dy[3] = -y[0]*np.sin(y[2])
            dy[4] = y[0]*np.cos(y[2])/(1 + y[3]/self.Rp)

            Yt = self.rhoa(y[3])*y[0]**2

            if Yt <= strength:
                dy[5] = 0.0
            else:
                dy[5] = np.sqrt(3.5*self.alpha*rhoa/density)*y[0]

            return dy

        def pallas_simple(t, y):

            """ Fuction to return the RHS of the asteroid equations
            Paramters
            ---------
            t : scalar time
                Current time during integration

            y : np.array
                Statevector array

            Note
            ----
            y[0] = v        Velocity of the asteroid
            y[1] = m        Mass of the asteroid
            y[2] = theta    Angle in the entry plane
            y[3] = z        Altitude
            y[4] = x        Range
            y[5] = r        Radius of the asteroid
            """

            dy = np.zeros_like(y)

            rhoa = self.rhoa(y[3])
            A = np.pi*y[5]**2

            dy[0] = -0.5*self.Cd*rhoa*A*y[0]**2/y[1]
            dy[1] = 0.0
            dy[2] = 0.0
            dy[3] = -y[0]*np.sin(y[2])
            dy[4] = y[0]*np.cos(y[2])

            return dy

        # Defining stopping conditions
        def stopVelocity(t, y):
            if y[0] < 0.:
                # print("Velocity is 0")
                return True
            else:
                return False

        def stopMass(t, y):
            if y[1] < 0.:
                # print("Mass is 0")
                return True
            else:
                return False

        def stopHeight(t, y):
            if y[3] < 0.:
                # print("Height is 0")
                return True
            else:
                return False

        def stopRange(t, y):
            if np.abs(y[4]) > 2*self.Rp:
                # print("Range is greater than earth diameter")
                return True
            else:
                return False

        def stopH(t, y):
            if y[3] > init_altitude:
                return True
            else:
                return False

        # Enter your code here to solve the differential equations
        t0 = 0.0
        tf = 1e10
        v0 = velocity
        m0 = (density*4*np.pi*radius**3)/3.0
        if radians:
            theta0 = angle
        else:
            theta0 = angle*np.pi/180.0
        z0 = init_altitude
        x0 = 0.0
        r0 = radius
        y0 = [v0, m0, theta0, z0, x0, r0]

        stoppingConditions = [stopVelocity, stopMass, stopHeight,
                              stopRange, stopH]

        pallas = SolveODE(pallas_rhs, (t0, tf), y0,
                          method=method,
                          stopping=stoppingConditions,
                          tol=1e-6, maxstep=dt)

        pallas.solve()

        time_out = np.arange(t0, pallas.t[-1]+dt, dt)
        velocity_out = np.interp(time_out, pallas.t, pallas.sol[:, 0])
        mass_out = np.interp(time_out, pallas.t, pallas.sol[:, 1])
        angle_out = np.interp(time_out, pallas.t, pallas.sol[:, 2])
        if not radians:
            angle_out = angle_out * 180.0/np.pi
        altitude_out = np.interp(time_out, pallas.t, pallas.sol[:, 3])
        distance_out = np.interp(time_out, pallas.t, pallas.sol[:, 4])
        radius_out = np.interp(time_out, pallas.t, pallas.sol[:, 5])

        return pd.DataFrame({'velocity': velocity_out,
                             'mass': mass_out,
                             'angle': angle_out,
                             'altitude': altitude_out,
                             'distance': distance_out,
                             'radius': radius_out,
                             'time': time_out}, index=range(len(time_out)))

    def calculate_energy(self, result):
        """
        Function to calculate the kinetic energy lost per unit altitude in
        kilotons TNT per km, for a given solution.

        Parameters
        ----------
        result : DataFrame
            A pandas dataframe with columns for the velocity, mass, angle,
            altitude, horizontal distance and radius as a function of time

        Returns : DataFrame
            Returns the dataframe with additional column ``dedz`` which is the
            kinetic energy lost per unit altitude

        """

        # Replace these lines with your code to add the dedz column to
        # the result DataFrame

        def central_diff(arr):
            dfor = np.zeros_like(arr)
            dback = np.zeros_like(arr)
            dfor[:-1] = arr[1:]
            dfor[-1] = arr[-1]
            dback[1:] = arr[:-1]
            dback[0] = arr[0]
            return dfor-dback

        kineticEnergy = (0.5*result["mass"]*result["velocity"]**2).values
        dz = central_diff(result["altitude"].values)*1e-3
        dE = central_diff(kineticEnergy)/(4.184*1e12)

        dEdz = dE/dz

        result = result.copy()
        result.insert(len(result.columns),
                      'dedz', dEdz)

        return result

    def analyse_outcome(self, result):
        """
        Inspect a pre-found solution to calculate the impact and airburst stats

        Parameters
        ----------
        result : DataFrame
            pandas dataframe with velocity, mass, angle, altitude, horizontal
            distance, radius and dedz as a function of time

        Returns
        -------
        outcome : Dict
            dictionary with details of the impact event, which should contain
            the key:
                ``outcome`` (which should contain one of the
                following strings: ``Airburst`` or ``Cratering``),
            as well as the following 4 keys:
                ``burst_peak_dedz``, ``burst_altitude``,
                ``burst_distance``, ``burst_energy``
        """
        def linearRoot(t1, y1, t2, y2):
            return (y2*t1 - y1*t2)/(y2 - y1)

        def linearPredict(t, t1, y1, t2, y2):
            return ((y2-y1)/(t2-t1))*t + (y1*t2-y2*t1)/(t2-t1)

        burst_peak_dedz = result.dedz.max()
        idx = result.dedz.idxmax()
        burst_altitude = result.loc[idx, "altitude"]

        if np.isclose(burst_altitude, 0.0):
            outcome = 'Cratering'

            burst_peak_dedz = result.loc[idx, "dedz"]
            burst_altitude = 0.0
            burst_distance = result.loc[idx, "distance"]
            Eh = 0.5*result.loc[idx, "mass"]*result.loc[idx, "velocity"]**2
            E0 = 0.5*result.loc[0, "mass"]*result.loc[0, "velocity"]**2
            burst_energy = max((E0-Eh), Eh)/(4.184*1e12)

        elif burst_altitude < 0.0:
            outcome = 'Cratering'

            # Use two datapoints for linear interpolation and figure out
            # the time where h = 0. Use this to calculate other required values

            t2 = result.loc[idx, "time"]
            t1 = result.loc[idx-1, "time"]

            h2 = result.loc[idx, "altitude"]
            h1 = result.loc[idx-1, "altitude"]

            x2 = result.loc[idx, "distance"]
            x1 = result.loc[idx-1, "distance"]

            e2 = result.loc[idx, "dedz"]
            e1 = result.loc[idx-1, "dedz"]

            v2 = result.loc[idx, "velocity"]
            v1 = result.loc[idx-1, "velocity"]

            m2 = result.loc[idx, "mass"]
            m1 = result.loc[idx-1, "mass"]

            th0 = linearRoot(t1, h1, t2, h2)
            xh0 = linearPredict(th0, t1, x1, t2, x2)
            eh0 = linearPredict(th0, t1, e1, t2, e2)
            vh0 = linearPredict(th0, t1, v1, t2, v2)
            mh0 = linearPredict(th0, t1, m1, t2, m2)

            burst_peak_dedz = eh0
            burst_altitude = 0.0
            burst_distance = xh0
            Eh = 0.5*mh0*vh0**2
            E0 = 0.5*result.loc[0, "mass"]*result.loc[0, "velocity"]**2
            burst_energy = max((E0-Eh), Eh)/(4.184*1e12)

        else:
            outcome = 'Airburst'

            burst_peak_dedz = result.loc[idx, "dedz"]
            burst_altitude = result.loc[idx, "altitude"]
            burst_distance = result.loc[idx, "distance"]
            Eh = 0.5*result.loc[idx, "mass"]*result.loc[idx, "velocity"]**2
            E0 = 0.5*result.loc[0, "mass"]*result.loc[0, "velocity"]**2
            burst_energy = max((E0-Eh), Eh)/(4.184*1e12)

        outcome = {'outcome': outcome,
                   'burst_peak_dedz': burst_peak_dedz,
                   'burst_altitude': burst_altitude,
                   'burst_distance': burst_distance,
                   'burst_energy': burst_energy}
        return outcome


if __name__ == "__main__":
    import doctest
    test = doctest.testmod()
    assert test[0] == 0
