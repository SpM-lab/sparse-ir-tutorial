from copy import deepcopy
import numpy as np


class Solver:
    """ FLEX solver
    An object of this class should be immutable.
    In other words, calling any member function must not modify the object.
    """
    def __init__(self, basis_set, nk1, nk2):
        self.basis_set = basis_set
        self.nk1, self.nk2, self.nk = nk1, nk2, nk1 * nk2
        self.iw0_f = **
        self.iw0_b = **
        self.k1, self.k2 = **
        self.iwn_f = ***
        self.ek = * 
        self.ek_ = *
        self.T = **
    
    def smpl(self, statistics):
        """ Return sampling objects for given statistics """
        smpl_tau = {"F": self.basis_set.smpl_tau_f, "B": self.basis_set.smpl_tau_b}[statistics]
        smpl_wn = {"F": self.basis_set.smpl_wn_f, "B": self.basis_set.smpl_wn_b}[statistics]
        return smpl_tau, smpl_wn
    
    def tau_to_wn(self, statistics, obj_tau):
        """ From tau to wn """
        smpl_tau, smpl_wn = self.smpl(statistics)
        obj_tau = obj_tau.reshape((smpl_tau.tau.size, self.nk1, self.nk2))
        obj_l = smpl_tau.fit(obj_tau, axis=0)
        return smpl_wn.evaluate(obj_l, axis=0).reshape((smpl_wn.wn.size, self.nk))
    
    def k_to_r(self, obj_k):
        """ Fourier transfom from k space to real space """
        pass

    def r_to_k(self, obj_r):
        pass
        
    
    def gkio_calc(self, mu, sigma):
        return 1/(self.iwn_f  - (self.ek - mu) - sigma)

    def grit_calc(self, gkio):
        # do something
    
    def ckio_calc(self, sigma):
        # do something
        return ckio

    def solve(self, U, simga_init, maxiter, sfc_tol):
        """
        Solve FLEX for a given U.
        The target accuracy can be controlled through the maxiter and scf_tol options.

        sigma_init: Initial guess for the self-energy
        maxiter: Max number of iterations
        """
        # First compute ckio and call renormalization_loop if needed
        ckio = ***
        if np.abs(ckio).max() >= 1:
            sigma = self.renormalization_loop(...)
        
        # Do actual SCF loops

        return sigma # Returning only simga may be enough because other quanties can be reconstructed from it.
    
    def renormalization_loop(self, n, U_final, sigma_init, U_it_max=100):
        """
        Renormalizatoin loop
        """
        U = U_final
        sigma = sigma_init.copy()
        mu = self.mu_calc(n, sigma)
        gkio = self.gkio_calc(mu, sigma)
        for U_it in range(U_it_max):
            # Renormalize U and compute ckio
            if U_final * np.amax(np.abs(ckio)) < 1:
                break

        return sigma
