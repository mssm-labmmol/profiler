#
# This file is part of the profilerTools suite (see
# https://github.com/mssm-labmmol/profiler).
#
# Copyright (c) 2020 mssm-labmmol
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
from .multiprofile import multiProfile
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_squared_error
from .readopts import MaskLooper
from .geometry_factor import calculate_geometry_factor


class GSLMultifitMixin:
    """This is a custom Mixin that may extend the functionality of
    LinearRegression with a method ``fit_gsl`` that relies on the GSL C API."""
    def fit_gsl(self, *args, **kwargs):
        import gslpyridge 
        self.coef_ = gslpyridge.fit(*args, **kwargs)


class ProfilerLinearRegression (LinearRegression, GSLMultifitMixin):
    def fit(self, reg=False, *args, **kwargs):
        """This wrapper calls the GSL implementation when ``reg`` argument is
        True."""
        if reg:
            return super().fit_gsl(*args, **kwargs)
        else:
            # In all other cases, rely on the usual fit.
            #return super().fit_gsl(lamb=0, *args, **kwargs)
            return super().fit(*args, **kwargs)


class LinearParameterOptimizer:
    """This class implements the LLS solution of the parameter optimization
    problem. It is actually a driver class to construct the LLS problem in
    matrix form and solve it.
    """
    def __init__(self, mp):
        """Creates a LinearParameterOptimizer instance.

        :param mp: A :class:`~profilerTools.multiprofile.multiProfile` object
                   that is responsible for a proper communication between this
                   Optimizer and the molecular-mechanics calculations.
        """
        self.mp = mp
        self.Sm = None # This is constructed below.
        self._build_matrix_S()
        self.model = None # This is initialized only when fitting.


    def _build_matrix_S(self):
        """Constructs the S matrix containing the isolated contribution of the
        energy-shift parameters. The number of rows of S is the total number of
        configurations, and the number of columns is the number of systems. An
        entry :math:`s_{ij}` is 1 if the i-th configuration belongs to the j-th
        system, and 0 otherwise.
        
        It sets `self.Sm`, a :class:`np.ndarray` with shape (M, N), where M is
        the total number of configurations and N is the number of systems.
        """
        confs = self.mp.calculateNumberOfConfigurations()
        # Start with the first line---this is '1' in the first column, followed
        # by as much '0's as necessary reach the number of systems.
        self.Sm = np.zeros((1, len(confs)))
        self.Sm[0, 0] = 1
        # Then stack lines as they are needed.
        for i, nconf in enumerate(confs):
            new_row = np.zeros((1, len(confs)))
            new_row[0, i] = 1
            if (i == 0):
                # For the first system, skip the first row since it is filled
                # during initialization.
                start = 1
            else:
                start = 0
            for conf in range(start, nconf):
                self.Sm = np.concatenate((self.Sm, new_row), axis=0)


    def _build_matrix_A(self):
        """Constructs the A matrix containing the geometric factors for the
        fitting. The number of rows of A is the number of configurations, and
        the number of columns is the number of optimized parameters (not
        considering the energy-shift parameters)

        :returns Am: A :class:`np.ndarray` with shape (M, N), where M is the
                     number of configurations and N is the number of optimized
                     parameters.

        """
        confs = self.mp.calculateNumberOfConfigurations()
        npars = self.mp.getNumberOfOptimizableParameters()
        k_phi_pairs = self.mp.getIndexConverter().get_k_phi_pairs()

        Am = np.zeros((np.sum(confs), npars))

        for k, (conf, profile) in enumerate(self.mp.getConfigurationsAndProfiles()):
            for p in range(npars):
                type_idx, selstring = self.mp.getIndexConverter().global_to_type(p)

                isTorsional = self.mp.getIndexConverter().global_is_torsional(p)

                dof_atom_idxs = profile.getAtomIdxsForOptParameter(
                    type_idx, isTorsional)

                Am[k,p] += np.sum(
                    calculate_geometry_factor(
                        conf, selstring, dof_atom_idxs, self.mp.dihType))

        # Fix geometry factors in case of optimization of force constant only.
        for k, (conf, profile) in enumerate(self.mp.getConfigurationsAndProfiles()):
            for kp, phi in k_phi_pairs:
                if phi is None:
                    # This means that the optimization is of force constant
                    # only, so the geometry factors must be corrected to
                    # account for the phase shifts initialized in the
                    # multiprofile.
                    ktp, name = self.mp.getIndexConverter().global_to_type(kp)
                    km = int(name.split('_')[1]) - 1
                    # It suffices to take one dihedral from the list of
                    # dihedrals associated with `kp', because the same
                    # correction factor multiplies each of them.
                    dih = profile.indexList.get(ktp)[0]

                    kidx = profile.mmCalc.idx_to_OptIdx(dih)
                    corr_factor = np.cos(np.radians(profile.mmCalc.optDihedralTerms.phi[km, kidx]))
                    Am[k, kp]  *= corr_factor

        return Am


    def _build_matrix_X(self):
        """Constructs the linear-system geometric matrix."""
        return np.concatenate((self._build_matrix_A(), self.Sm), axis=1)

    
    def _build_Y(self, target_data):
        """Constructs the Y vector by subtracting the MM energies from the
        target data.

        :param target_data: A :class:`np.ndarray` with shape (K, 1), where K is
                            the total number of configurations. It contains the
                            target data for fitting.
        """
        return (np.reshape(target_data, (-1,1)) - self.mp.getNonoptEnergy())

    def _calc_S_reg_center(self, Y):
        """Calculates regularization center for the energy-shift parameters
        based on the Y vector."""
        avgs = []
        confs = self.mp.calculateNumberOfConfigurations()
        start = 0
        stop  = 0
        for nconf in confs:
            stop += nconf
            avgs.append(-np.mean(Y[start:stop,0]))
            start += nconf
        avgs = np.array(avgs)
        return avgs

    
    def _fixup_reg_center(self, Y, reg_center):
        """Returns a padded version of `reg_center` as a column-vector."""
        if reg_center is None:
            return None
        S_reg_center = self._calc_S_reg_center(Y)
        ext_reg_center = np.concatenate((reg_center, S_reg_center), axis=0)
        return ext_reg_center.reshape(-1, 1)

    
    def _add_reg_to_Y(self, Y, X, reg_center):
        """Modify in-place the numpy array `Y` so that it includes the
        contribution of the regularization center."""
        if reg_center is None:
            return
        ext_reg_center = self._fixup_reg_center(Y, reg_center)
        Y -= np.matmul(X, ext_reg_center)

        
    def _prepare_Ldiag(self, Y, reg_center):
        return 1.0 / self._fixup_reg_center(Y, reg_center)


    def _add_reg_to_coefs(self, Y, reg_center):
        """Adds regularization center to current coefficients."""
        if reg_center is None:
            return
        ext_reg_center = self._fixup_reg_center(Y, reg_center)
        self.coef_ += ext_reg_center

        
    def get_number_of_parameters(self):
        """Returns the number of parameters in optimization."""
        return self.mp.getNumberOfOptimizableParameters()


    def get_parameters(self):
        """Returns the current parameters."""
        return self.coef_.flatten()[:-self.mp.getNumberOfSystems()]

    
    def fit(self, target_data, wei=None, reg_center=None):
        """Fits a linear (or Ridge) regression model to the target data,
        possibly using custom weights.

        :param target_data: A :class:`np.ndarray` with shape (K, 1), where K is
                            the total number of configurations. It contains the
                            target data for fitting.

        :param wei: A :class:`np.ndarray` with shape (K, 1), where K is the
                    total number of configurations. It contains the weights
                    attributed to each data point. If `None`, all weights are 1.

        :param reg_center: A :class:`np.ndarray` with shape (N,), where N is the
                           number of optimized parameters. It contains the point
                           in parameter space relative to which the
                           regularization distance is computed. If `None`, a
                           regular linear regression, without regularization, is
                           carried out.
        """
        if wei is None:
            wei = np.zeros(target_data.shape)
            wei.fill(1)

        self.model = ProfilerLinearRegression()

        X = self._build_matrix_X()
        Y = self._build_Y(target_data)

        self._add_reg_to_Y(Y, X, reg_center)

        if reg_center is None:
            self.model.fit(X=X, y=Y, sample_weight=wei)
        else:
            ldiag = self._prepare_Ldiag(Y, reg_center)
            self.model.fit(X=X, y=Y, sample_weight=wei, reg=True, ldiag=ldiag)

        self.coef_ = self.model.coef_

        self._add_reg_to_coefs(Y, reg_center)

        
class LLS_SC:
    """
    Class responsible for carrying out the self-consistent cycles of linear
    regression.
    """
    class LLS_SC_Logbook:
        """A nested class that serves to store the results of a LLS_SC run. It
        is essentially a dictionary with specific keys, and each element is a
        list of the quantity identified by the key along the self-consistent
        cycles."""
        def __init__(self):
            self.clear()

        def add_to_quantity(self, key, value):
            self.data[key] = value

        def clear(self):
            self.data = dict()
            self.data['wrmsd'] = []
            self.data['parameters'] = []


    def __init__(self, multiprofile):
        self.regressor = LinearParameterOptimizer(multiprofile)
        self.mp = multiprofile
        self.logbook = self.LLS_SC_Logbook()


    def run(self, max_cycles, max_dpar, target_data, wei=None, reg_center=None):
        """Runs the self-consistent scheme for at most `max_cycles` or until the
        maximum relative change in the values of the parameters is at most
        `max_dpar`.

        :param max_cycles: (int) Maximum number of cycles.

        :param max_dpar: (float) Tolerance for maximum relative change in the
                         values of the parameters.

        :param target_data: A :class:`np.ndarray` with shape (K, 1), where K is
                            the total number of configurations. It contains the
                            target data for fitting.

        :param wei: A :class:`np.ndarray` with shape (K, 1), where K is the
                    total number of configurations. It contains the weights
                    attributed to each data point. If `None`, all weights are 1.

        :param reg_center: A :class:`np.ndarray` with shape (N,), where N is the
                           number of optimized parameters. It contains the point
                           in parameter space relative to which the
                           regularization distance is computed. If `None`, a
                           regular linear regression, without regularization, is
                           carried out.

        :returns: `True` if the parameters converged before the maximum number
                  of cycles was reached; `False` otherwise.

        """
        self.logbook.clear()
        npars  = self.regressor.get_number_of_parameters()
        prev_pars = np.zeros((npars,))
        curr_pars = np.zeros((npars,))
        prev_pars.fill(1.0)
        curr_pars.fill(10.0)
        cycles = 0
        while ((cycles < max_cycles) and
               (np.abs(np.max((curr_pars-prev_pars)/prev_pars)) >= max_dpar)):
            print(f"LLS-SC: Running cycle {cycles+1}... ", end="")

            self.regressor.fit(target_data, wei, reg_center)
            
            prev_pars = np.copy(curr_pars)
            curr_pars = self.get_parameters()

            self.mp.setOptimizableParameters(slice(npars), curr_pars)

            if self.mp.areThereUnphysicalParameters():
                raise ValueError("Couldn't obtain physical parameters.")

            self.mp.minimizeProfiles()
            
            # Update Logbook
            wrmsd = self.mp.rmsdToData()
            self.logbook.add_to_quantity('parameters', curr_pars)
            self.logbook.add_to_quantity('wrmsd', wrmsd)

            cycles += 1

            print("Done.")
            #print(f"Parameters:")
            #for p in curr_pars:
            #    print("%14.6e" % p)
            print(f"WRMSD:")
            print("%14.6e" % wrmsd)

        return not (cycles == max_cycles)

    @staticmethod
    def hr2standard(A_m, B_m):
        if A_m is None:
            return None, B_m
        if B_m is None:
            return A_m, None
        phi = np.arctan2(B_m, A_m)
        k = A_m/np.cos(phi) # this form is sign-aware
        return k, np.degrees(phi)

    def get_parameters(self):
        """Returns the optimized parameters, converting from Hopkins-Roitberg to
        standard form if necessary."""
        pars = self.regressor.get_parameters()
        k_phi_pairs = self.mp.getIndexConverter().get_k_phi_pairs()
        for k, phi in k_phi_pairs:
            if k is None:
                A_m = None
            else:
                A_m = pars[k]
            if phi is None:
                B_m = None
            else:
                B_m = pars[phi]
            k_new, phi_new = self.hr2standard(A_m, B_m)
            if k is not None:
                pars[k] = k_new
            if phi is not None:
                pars[phi] = phi_new
        return pars
