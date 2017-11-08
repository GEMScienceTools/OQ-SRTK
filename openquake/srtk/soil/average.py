# =============================================================================
#
# Copyright (C) 2010-2017 GEM Foundation
#
# This file is part of the OpenQuake's Site Response Toolkit (OQ-SRTK)
#
# OQ-SRTK is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# OQ-SRTK is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# Author: Valerio Poggi
#
# =============================================================================
"""
The module contains several functions to derive and manipulate average
soil parameters, such as travel-time average velocity and site kappa.
"""

import numpy as _np
import scipy.optimize as _spo


# =============================================================================

def depth_weighted_average(thickness, soil_param, depth):
    """
    Compute the weighted average of a soil property at
    arbitrary depth.

    :param numpy.array tickness:
        array of layer's thicknesses in meters (half-space is 0.)

    :param numpy.array soil_param:
        array of soil properties (e.g. slowness, density)

    :param float depth:
        averaging depth in meters

    :return float mean_param:
        the weighted mean of the given soil property
    """

    mean_param = 0
    total_depth = 0

    for tk, sp in zip(thickness[:-1], soil_param[:-1]):
        if (tk + total_depth) < depth:
            mean_param += tk*sp/depth
        else:
            mean_param += (depth - total_depth)*sp/depth
            break
        total_depth += tk

    # Check for the half-space
    if total_depth == _np.sum(thickness[:-1]):
        mean_param += (depth - total_depth)*soil_param[-1]/depth

    return mean_param


# =============================================================================

def traveltime_average_velocity(thickness, s_velocity, depth=30):
    """
    The function calculates the travel-time average (harmonic mean)
    velocity at arbitrary depth (e.g. the widespread Vs30).

    :param numpy.array tickness:
        array of layer's thicknesses in meters (half-space is 0.)

    :param numpy.array s_velocity:
        array of layer's shear-wave velocities in m/s

    :param float depth:
        averaging depth in meters; if depth is not specified,
        depth is fixed to 30m

    :return float mean_velocity:
        the average velocity in m/s
    """

    # Converting velocity to slowness
    slowness = 1./s_velocity

    # Harmonic averaging is done in slowness
    mean_slowness = depth_weighted_average(thickness, slowness, depth)

    # Back to velocity
    mean_velocity = 1./mean_slowness

    return mean_velocity


# =============================================================================

def compute_site_kappa(thickness, s_velocity, s_quality, depth=[]):
    """
    This function calucalte the site attenuation parameter Kappa(0)
    for a given soil profile at arbitrary depth.

    :param numpy.array tickness:
        array of layer's thicknesses in meters (half-space is 0.)

    :param numpy.array s_velocity:
        array of layer's shear-wave velocities in m/s

    :param numpy.array s_quality:
        array of layer's shear-wave quality factors (adimensional)

    :param float depth:
        averaging depth in meters; if depth is not specified,
        the last layer interface is used instead

    :return float kappa0:
        the site attenuation parameter kappa(0) in seconds
    """

    # If depth not given, using the whole profile
    if not depth:
        depth = _np.sum(thickness)

    # Kappa vector
    layer_kappa = depth/(s_velocity*s_quality)

    # Computing the average kappa of the profile
    kappa0 = depth_weighted_average(thickness, layer_kappa, depth)

    return kappa0


# =============================================================================

def quarter_wavelength_velocity(thickness, s_velocity, frequency):
    """
    This function solves the quarter-wavelength problem (Boore 2003)
    and return the frequency-dependent average velocity.

    :param numpy.array tickness:
        array of layer's thicknesses in meters (half-space is 0.)

    :param numpy.array s_velocity:
        array of layer's shear-wave velocities in m/s

    :param numpy.array frequency:
        array of frequencies in Hz for the calculation

    :return numpy.array qwl_depth:
        array of averaging depths

    :return numpy.array qwl_velocity:
        array of quarter-wavelength average velocities
    """

    # Initialisation
    freq_num = len(frequency)
    slowness = 1./s_velocity
    qwl_depth = _np.zeros(freq_num)
    qwl_velocity = _np.zeros(freq_num)

    for nf in range(freq_num):

        # Upper depth bound for the search
        ubnd = _np.max(1./(4.*frequency[nf]*slowness))

        # Input arguments for the search function
        args = (thickness, slowness, frequency[nf])

        # Compute the quarter-wavelength depth
        qwl_depth[nf] = _spo.fminbound(_qwl_fit_func, 0., ubnd, args)

        # Computing average soil property at the qwl-depth
        qwl_velocity[nf] = 1./depth_weighted_average(thickness,
                                                     slowness,
                                                     qwl_depth[nf])

    return qwl_depth, qwl_velocity


# =============================================================================

def _qwl_fit_func(search_depth, thickness, slowness, frequency):
    """
    Internal function to recursively search for the qwl depth.
    """

    qwl_slowness = depth_weighted_average(thickness,
                                          slowness,
                                          search_depth)

    # Misfit is computed as a simple L1 norm
    misfit = _np.abs(search_depth - (1./(4.*frequency*qwl_slowness)))

    return misfit
