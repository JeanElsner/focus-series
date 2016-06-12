#!/usr/bin/env python

__project__  = 'focus-series'
__author__   = 'Jean Elsner'
__version__  = '1.0.0'
__email__    = 'jean.elsner@googlemail.com'
__homepage__ = 'https://github.com/JeanElsner/focus-series'

knots = ([-8.19125032,  -8.19125032,  -8.19125032,  -8.19125032, 
          16.28125048,  16.28125048,  16.28125048,  16.28125048], 
         [14.16073385,  14.16073385,  14.16073385,  14.16073385, 
          66.94437713,  66.94437713,  66.94437713,  66.94437713])

coeffs = [2.53674674,  2.59941779,  2.44225603,  2.7434,
          2.22122555,  2.14744304,  2.36902851,  2.17264242, 
          2.086511,    2.11510915,  2.15300977,  2.19045497,
          1.67541365,  1.69895363,  1.63709489,  1.73183233]

def _knot_vector_slopes(u_knots, v_knots, coeffs, k, l):
    u_slope = eval_surf_spline(u_knots[-1], v_knots[0], u_knots, v_knots, coeffs, k, l)
    u_slope -= eval_surf_spline(u_knots[0], v_knots[0], u_knots, v_knots, coeffs, k, l)
    u_slope /= (u_knots[-1] - u_knots[0])

    v_slope = eval_surf_spline(u_knots[0], v_knots[-1], u_knots, v_knots, coeffs, k, l)
    v_slope -= eval_surf_spline(u_knots[0], v_knots[0], u_knots, v_knots, coeffs, k, l)
    v_slope /= (v_knots[-1] - v_knots[0])
    return u_slope, v_slope

def eval_surf_spline (u, v, u_knots, v_knots, coeffs, k = 4, l = 4):
    r"""Evaluates the B-Spline surface at (u,v).
    
    Evaluates the surface as defined by the knot vectors
    and the coefficients at the point (u,v).
    
    Parameters
    ----------
    u, v : float
        Bivariate coordinates to evaluate.
    u_knots, v_knots : array_like
        Knot vectors as returned by the interpolation routine.
    coeffs : array_like
        Coefficients corresponding to the knot vectors.
    k, l : int
        Order of the spline.
    
    Returns
    -------
    val : float
        Value of the surface at (u,v).
    """
    if not min(u_knots) <= u <= max(u_knots) or not min(v_knots) <= v <= max(v_knots):
        u_nn = max(min(u_knots), min(max(u_knots), u))
        v_nn = max(min(v_knots), min(max(v_knots), v))
        p = eval_surf_spline(u_nn, v_nn, u_knots, v_knots, coeffs, k, l)
        u_slope, v_slope = _knot_vector_slopes(u_knots, v_knots, coeffs, k, l)
        
        return p + u_slope*(u-u_nn) + v_slope*(v-v_nn)
    
    m = len(u_knots) - k
    n = len(v_knots) - l
    val = 0
    for i in range(m):
        for j in range(n):
            u_basis = get_bspline_basis(u, i, k, u_knots)
            v_basis = get_bspline_basis(v, j, l, v_knots)
            val += u_basis * v_basis * coeffs[j+i*n]
    return val

def get_bspline_basis(t, i, k, knots):
    r"""Calculates the basis function for a given knot vector.
    
    Given a knot vector, the basis function as needed for a
    B-Spline is evaluated at position t.
    
    Parameters
    ----------
    t : float
        The point the basis function will be evaluated at.
    i, k : int
        The indices of the basis function.
    knots : array_like
        An array_like object used as the knot vector.
    
    Returns
    -------
    basis : float
        The basis function evaluated at t.
    """
    if k == 1:
        if knots[i] <= t <= knots[i+1]:
            return 1
        else:
            return 0
    t_i     = knots[i]
    t_i_1   = knots[i+1]
    t_i_k   = knots[i+k]
    t_i_k_1 = knots[i+k-1]
    basis = 0
    if (t_i_k_1 - t_i) != 0:
        basis += (t - t_i)/(t_i_k_1 - t_i)*get_bspline_basis(t, i, k-1, knots)
    if (t_i_k - t_i_1) != 0:
        basis += (t_i_k - t)/(t_i_k - t_i_1)*get_bspline_basis(t, i+1, k-1, knots)
    return basis

import argparse
 
parser = argparse.ArgumentParser(description='Focus interpolation based on best focus data.')
parser.add_argument('-t', '--temperature', help='Temperature to evaluate', required=True)
parser.add_argument('-zd' ,'--zenith-distance', help='Zenith-distance to evaluate', required=True)
args = parser.parse_args()

print(eval_surf_spline(float(args.temperature), float(args.zenith_distance), knots[0], knots[1], coeffs))