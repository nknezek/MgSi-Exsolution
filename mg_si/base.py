import numpy as np

const_yr_to_sec = 3.154e7 # seconds

class Parameters(object):
    def __init__(self, source):
        self.source = source

class Layer(object):
    '''
    The layer base class defines the geometry of a spherical shell within
    a planet.
    '''

    def __init__(self, inner_radius, outer_radius, params={}):
        self.set_boundaries(inner_radius, outer_radius)
        self.params = params

    def set_boundaries(self, inner_radius, outer_radius):
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.thickness = outer_radius - inner_radius

        assert self.thickness >= 0.0, "Ri={0:.1f} km, Ro={1:.1f} km".format(self.inner_radius / 1e3,
                                                                            self.outer_radius / 1e3)

        self.inner_surface_area = 4.0 * np.pi * self.inner_radius ** 2.
        self.outer_surface_area = 4.0 * np.pi * self.outer_radius ** 2.

        self.volume = 4.0 / 3.0 * np.pi * (self.outer_radius ** 3. - self.inner_radius ** 3.)

