import meep as mp
import numpy as np

class MeepCartesianRevolution() :
    def __init__(self, zs = [-10,10], rs = [5,5], **kwargs):
        self.zs = zs
        self.rs = rs

    @property
    def xsize(self):
        return float(self.rs.max())

    @property
    def ysize(self):
        return float(self.rs.max())

    @property
    def zsize(self):
        return float(self.zs.max() - self.zs.min())

    @property
    def zs(self):
        return self._zs

    @zs.setter
    def zs(self, zs):

        # require at least two elements
        if len(zs) < 2:
            raise ValueError("zs must have at least 2 elements")

        # check dimension of array
        if type(zs) is np.array and len(zs.shape) != 1 :
            raise ValueError("zs must be 1 dimensional")

        # if a list upgrade to np.array
        if type(zs) is list :
            zs = np.array(zs)

        self._zs = zs

    @property
    def rs(self):
        return self._rs

    @rs.setter
    def rs(self, rs):

        # require at least two elements
        if len(rs) < 2 :
            raise ValueError("rs must have at least 2 elements")

        # check dimension of array
        if type(rs) is np.array and len(rs.shape) != 1 :
            raise ValueError("rs must be 1 dimensional")

        # if a list upgrade to np.array
        if type(rs) is list :
            rs = np.array(rs)

        self._rs = rs

class MeepCylindricalRevolution() :
    def __init__(self):
        pass
