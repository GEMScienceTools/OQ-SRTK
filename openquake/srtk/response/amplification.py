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
Collection of methods to compute seismic amplification from a
one-dimensional soil column.
"""

import numpy as _np


# =============================================================================

def impedance_amplification(top_vs, top_dn, ref_vs=[], ref_dn=[], inc_ang=0.):
    """
    This function calculates the amplification due to a single seismic
    impedance contrast as formalized in Joyner et al. (1981) and
    Boore (2013).
    Site (top) parameters can be either scalars or arrays, as for
    the case of quarter-wavelenght approximation. If no reference
    is provided, the last value of the site array is used.
    Default angle of incidence is vertical.

    :param float or numpy.array top_vs:
        topmost shear-wave velocity in m/s

    :param float or numpy.array top_dn:
        topmost density in kg/m3

    :param float or numpy.array ref_vs:
        lowermost (reference) shear-wave velocity in m/s

    :param float or numpy.array ref_dn:
        lowermost (reference) density in kg/m3

    :param float inc_ang:
        angle of incidence relative to the vertical direction in degree

    :return float or array imp_amp:
        amplification due to the seismic impedance contrast
    """

    # Setting the reference, if not provided
    if not _np.sum(ref_vs):
        ref_vs = top_vs[-1] if _np.size(top_vs) > 1. else top_vs

    if not _np.sum(ref_dn):
        ref_dn = top_dn[-1] if _np.size(top_dn) > 1. else top_dn

    # Computing square-root impedance amplification
    imp_amp = _np.sqrt((ref_dn*ref_vs)/(top_dn*top_vs))

    if inc_ang > 0.:
        # Converting incident angle from degrees to radiants
        inc_ang *= _np.pi/180.

        # Effective angle of incidence computed using Snell's law
        eff_ang = _np.arcsin((top_vs/ref_vs)*_np.sin(inc_ang))

        # Correcting for non-vertical incidence
        imp_amp *= _np.sqrt(_np.cos(inc_ang)/_np.cos(eff_ang))

    return imp_amp
