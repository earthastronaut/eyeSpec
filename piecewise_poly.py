"""A class for representing and fitting piecewise polynomial functions with and without
regularity constraints (such as requiring continuity between adjacent polynomial knots.
"""

import numpy as np
lna = np.linalg
poly1d = np.poly1d
import matplotlib.pyplot as plt
#Legendre = np.polynomial.legendre.Legendre

class Centered_Poly:
    
    def __init__(self, coefficients, center):
        self.poly = poly1d(coefficients)
        self.center = center
    
    def __call__(self, xdat):
        return self.poly(xdat-self.center)
    
    def deriv(self):
        return Centered_Poly(self.poly.deriv().c, self.center)

class Binning:

    def __init__(self, bins):
        self.bins = bins
        self.lb = bins[0]
        self.ub = bins[-1]
        self.last_bin = bins[0], bins[1]
        self.last_bin_idx = 0
        #TODO build an efficient binary search tree structure for the bins

    def get_bin_index(self, xvec):
        xv = np.array(xvec)
        out_idxs = np.zeros(len(xv.flat))
        for x_idx in xrange(len(xv.flat)):
            llb, lub = self.last_bin
            if llb <= xv[x_idx] <= lub:
                out_idxs[x_idx] = self.last_bin_idx
            elif self.lb < xv[x_idx] <= self.ub:
                for test_idx in xrange(len(self.bins)-1):
                    clb, cub = self.bins[test_idx], self.bins[test_idx+1]
                    if clb < xv[x_idx] <= cub:
                        out_idxs[x_idx] = test_idx
                        self.last_bin_idx = test_idx
                        self.last_bin = clb, cub
                        break
            elif xv[x_idx] == self.lb:
                out_idxs[x_idx] = 0 #make the global lower bound inclusive
            else:
                out_idxs[x_idx] = np.nan
        return out_idxs

class Piecewise_Polynomial:
    
    def __init__(self, coefficients, control_points, centers = None, bounds = (float("-inf"), float("inf")), fill_value = np.nan):
        """represents a piecewise polynomial function which transitions from one polynomial
        to the next at the control points.
        coefficients should be an (m, n) array
        m is the number of polynomial pieces == len(control_points) + 1
        n is the order of the polynomial pieces
        
        The function takes on values which are determined by the polynomial coefficients
        with the highest order terms coming first and each polynomail being centered around
        either the corresponding value in the centers array if it is passed as an argument
        By default the center is chosen as the midpoint of its two bounding points. 
        If one of the current bounding points is + or -infinity the other bounding point 
        is taken as the "center" of that polynomial bin
        
        Example:
        coefficients = np.array([[3, 2], [1, 0], [-1, -1]]) control_points = [5, 6]
        and bounds = (-float('inf'), 8)
        
        because the centers are 
        
        would be evaluated at a point x < 5 as 
        y = 3*(x-5) + 2
        
        and at a point 5 < x < 6 
        
        y = 1*(x-4.5) + 0
        
        and at a point 6 < x < 8
        
        y = -1*(x-7) + -1
        
        TODO: make this example an interactive session
        """
        self.coefficients = coefficients
        self.bounds = bounds
        self.control_points = control_points
        n_polys, poly_order = coefficients.shape
        self.poly_order = poly_order
        self.ncp = len(control_points)  
        boundary_points = np.zeros(self.ncp+2)
        boundary_points[0] = bounds[0]
        boundary_points[-1] = bounds[1]
        boundary_points[1:-1] = control_points
        self.binning = Binning(boundary_points)
        self.n_polys = n_polys
        self.fill_value = fill_value
        if centers == None:
            self.centers = np.zeros(n_polys)
            #set the centers in such a way to allow for infinite bounds
            for center_idx in range(n_polys):
                lb = boundary_points[center_idx]
                ub = boundary_points[center_idx+1]
                if lb == float("-inf"):
                    lb = boundary_points[center_idx+1]
                if ub == float("inf"):
                    ub = boundary_points[center_idx]
                self.centers[center_idx] = 0.5*(lb+ub)
        else:
            self.centers = centers
        self.poly_list = []
        for poly_idx in range(n_polys):
            self.poly_list.append(Centered_Poly(coefficients[poly_idx], self.centers[poly_idx]))

    def __call__(self, x_in):
        output = np.zeros(x_in.shape)
        poly_idxs = self.binning.get_bin_index(x_in)
        output[np.isnan(poly_idxs)] = self.fill_value
        for p_idx in xrange(self.n_polys):
            pmask = poly_idxs == p_idx
            output[pmask] = self.poly_list[p_idx](x_in[pmask])
        return output
    
def fit_piecewise_polynomial(x, y, order, control_points, bounds = (float("-inf"), float("inf")), regularity_constraints = None, centers = None):
    pp_gen = Regularity_Constrained_Piecewise_Polynomial_Basis(order, control_points=control_points, bounds = bounds, regularity_constraints = regularity_constraints, centers = centers)
    gbasis = pp_gen.get_basis(x)
    n_polys = len(control_points) + 1
    n_coeffs = order+1
    out_coeffs = np.zeros((n_polys, n_coeffs))
    fit_coeffs = np.linalg.lstsq(gbasis.transpose(), y)[0]
    for basis_idx in xrange(pp_gen.n_basis):
        c_coeffs = pp_gen.basis_coefficients[basis_idx].reshape((n_polys, n_coeffs))
        out_coeffs += c_coeffs*fit_coeffs[basis_idx]
    return Piecewise_Polynomial(out_coeffs, control_points, centers, bounds)
    

class Regularity_Constrained_Piecewise_Polynomial_Basis:

    def __init__(self, poly_order, control_points, centers = None, regularity_constraints = None, bounds = (float("-inf"), float("inf"))):
        self.bounds = bounds
        self.control_points = control_points
        self.poly_order = poly_order
        self.ncp = len(control_points)
        if regularity_constraints == None:
            self.regularity_constraints = np.ones((poly_order, self.ncp), dtype = bool)
        else:
            self.regularity_constraints = regularity_constraints
        boundary_points = np.zeros(self.ncp+2)
        boundary_points[0] = bounds[0]
        boundary_points[-1] = bounds[1]
        boundary_points[1:-1] = control_points
        self.binning = Binning(boundary_points)
        n_polys = self.ncp+1
        self.n_polys = n_polys
        if centers == None:
            self.centers = np.zeros(n_polys)
            #set the centers in such a way to allow for infinite bounds
            for center_idx in range(n_polys):
                lb = boundary_points[center_idx]
                ub = boundary_points[center_idx+1]
                if lb == float("-inf"):
                    lb = boundary_points[center_idx+1]
                if ub == float("inf"):
                    ub = boundary_points[center_idx]
                self.centers[center_idx] = 0.5*(lb+ub)
        else:
            self.centers = centers
        poly_basis_list = [[] for i in range(n_polys)]
        for poly_i in range(n_polys):
            #cdomain = (self.boundary_points[poly_i], self.boundary_points[poly_i+1])
            for comp_i in range(poly_order+1):
                comp_vec = np.zeros((poly_order+1))
                comp_vec[comp_i] = 1.0
                #poly_basis_list[poly_i].append(Legendre(comp_vec, domain = cdomain)) 
                poly_basis_list[poly_i].append(Centered_Poly(comp_vec, self.centers[poly_i]))

        #generate the constraint matrix
        #nrows = self.poly_order*self.ncp
        nrows = np.sum(self.regularity_constraints)
        constraint_matrix = np.zeros((nrows, (self.poly_order+1)*self.n_polys))
        constraint_number = 0
        nco, ncp = self.regularity_constraints.shape
        for control_i in range(ncp):
            c_control_point = self.control_points[control_i]
            l_basis = poly_basis_list[control_i] #left basis functions
            r_basis = poly_basis_list[control_i+1] #right basis functions
            for constraint_order in range(nco):
                if not self.regularity_constraints[constraint_order, control_i]:
                    continue
                fp_coeff_idx = control_i*(self.poly_order+1)
                sp_coeff_idx = (control_i+1)*(self.poly_order+1)
                #print "cp", control_i, "sp i", sp_coeff_idx
                for coefficient_i in range(self.poly_order+1):
                    lreg_coeff = l_basis[coefficient_i](c_control_point)
                    rreg_coeff = r_basis[coefficient_i](c_control_point)
                    constraint_matrix[constraint_number, fp_coeff_idx+coefficient_i] = lreg_coeff
                    constraint_matrix[constraint_number, sp_coeff_idx+coefficient_i] = -rreg_coeff
                #go up to the next order constraint by taking the derivative of our basis functions
                constraint_number += 1
                l_basis = [cpoly.deriv() for cpoly in l_basis]
                r_basis = [cpoly.deriv() for cpoly in r_basis]
        self.constraint_matrix = constraint_matrix
        u, s, v = lna.svd(self.constraint_matrix, full_matrices=True)
        self.n_basis = (self.poly_order+1)*self.n_polys-nrows
        self.basis_coefficients = np.zeros((self.n_basis, self.n_polys, self.poly_order+1))
        self.basis_polys = [[] for bi in range(self.n_basis)]
        for basis_i in range(self.n_basis):
            for poly_i in range(self.n_polys):
                coeff_lb = (self.poly_order+1)*poly_i
                coeff_ub = coeff_lb + self.poly_order+1
                ccoeffs = v[-(basis_i+1)][coeff_lb:coeff_ub]
                self.basis_coefficients[basis_i, poly_i] = ccoeffs
                self.basis_polys[basis_i].append(Centered_Poly(ccoeffs, self.centers[poly_i]))
            
    def get_basis(self, xvec):
        poly_idxs = self.binning.get_bin_index(xvec)
        out_basis = np.zeros((self.n_basis, len(xvec)))
        for basis_idx in xrange(self.n_basis):
            for poly_idx in xrange(self.n_polys):
                xmask = poly_idxs == poly_idx 
                cx = xvec[xmask]
                out_basis[basis_idx][xmask] = self.basis_polys[basis_idx][poly_idx](cx)
        return out_basis

# if __name__ == "__main__":
#     test_x = np.linspace(-1, 1, 4000)
#     test_y = test_x * 2 - test_x**2 + 3.14
    
#     ppol = fit_piecewise_polynomial(test_x, test_y, 3, np.array([0.5, 0.5]))
#     fit_y = ppol(test_x)
#     if np.sum(np.abs(test_x-test_y)) <= 1e-10:
#         print "PASSED exact fit test"
#     else:
#         print "FAILED exact fit test"
    
    
