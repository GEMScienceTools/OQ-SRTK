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
Plotting utilities for data visualisation
"""

import numpy as _np
import matplotlib.pyplot as _plt

# =============================================================================
# Figure settings

FIG_SIZE_MODEL = (4, 6)
FIG_SIZE_RESPONSE = (8, 4)

MODEL_LABELS = {'hl': 'Thickness (m)',
                'vp': 'P-wave velocity (m/s)',
                'vs': 'S-wave velocity (m/s)',
                'dn': 'Density (kg/m3)',
                'qp': 'P-wave quality factor',
                'qs': 'S-wave quality factor'}

RESPONSE_LABELS = {'shtf': 'SH-wave amplification',
                   'kappa': 'Attenuation',
                   'qwl': 'Qwl amplification'}


# =============================================================================

def plot_models(site1d, key='vs', hold=False):
    """
    """

    if not hold:
        _plt.figure(figsize=FIG_SIZE_MODEL)

    for mod in site1d.model:
        plot_profile(mod, key, hold=True, show=False)

    _profile_decoration(key)


# =============================================================================

def plot_profile(model, key='vs', color='r', hold=False, show=True):
    """
    """

    K = model.geo[key]
    H = _np.cumsum(model.geo['hl'][:-1])

    layer_num = len(K)

    X = []
    Y = []
    for nl in range(layer_num):
        if nl == 0:
            X += [K[0], K[0]]
            Y += [0, H[0]]
            X += [K[0], K[1]]
            Y += [H[0], H[0]]
        else:
            if nl != (layer_num - 1):
                X += [K[nl], K[nl]]
                Y += [H[nl-1], H[nl]]
                X += [K[nl], K[nl+1]]
                Y += [H[nl], H[nl]]
            else:
                X += [K[nl], K[nl]]
                Y += [H[nl-1], H[nl-1]*1.4]

    if not hold:
        _plt.figure()

    _plt.plot(X, Y, color=color, linewidth=3)

    if show:
        _profile_decoration(key)


# =============================================================================

def _profile_decoration(key):
    """
    """

    _plt.grid('on')
    _plt.gca().invert_yaxis()
    _plt.xlabel(MODEL_LABELS[key])
    _plt.ylabel(MODEL_LABELS['hl'])
    _plt.draw_all()
    _plt.show(block=False)
