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

# =============================================================================

def sh_transfer_function(freq, hl, vs, dn, qs=[], inc_ang=0.):
    """
    Compute the SH-wave transfer function using Knopoff formalism
    (implicit layer matrix scheme). Calculation can be done for an
    arbitrary angle of incidence (0-90), with or without anelastic
    attenuation.
    NOTE: the implicit scheme is simple to understand and to implement,
    but is also computationally intensive and memory consuming.
    For the future, an explicit (recursive) scheme must be implemented.

    :param numpy.array hl:

    :param numpy.array vs:

    :param numpy.array dn:

    :param numpy.array qs:

    :param float or numpy.array freq:

    :param numpy.array freq:

    """

    # Precision of the complex type
    CTP = 'complex128'

    # Check for single frequency value
    if not isinstance(freq, list):
        freq = _np.array([freq])

    # Model size
    lnum = len(hl)
    fnum = len(freq)

    # Variable recasting to numpy complex
    hl = _np.array(hl, dtype=CTP)
    vs = _np.array(vs, dtype=CTP)
    dn = _np.array(dn, dtype=CTP)

    # Attenuation using complex velocities
    if qs:
        qs = _np.array(qs, dtype=CTP)
        vs *= ((2.*qs*1j)/(2.*qs*1j-1.))

    # Angular frequency conversion
    angf = 2.*_np.pi*freq

    # -------------------------------------------------------------------------
    # Computing angle of propagation within layers

    iD = _np.zeros((lnum, 1), dtype=CTP)
    iCORE = _np.zeros((lnum, lnum), dtype=CTP)

    iD[0] = _np.sin(inc_ang)
    iCORE[0,-1] = 1.

    for nl in range(lnum-1):
      iCORE[nl+1,nl] = 1./vs[nl]
      iCORE[nl+1,nl+1] = -1./vs[nl+1]

    iA = _np.linalg.solve(iCORE, iD)
    iS = _np.arcsin(iA)

    # -------------------------------------------------------------------------
    # Elastic parameters

    # Lame Parameters : shear modulus
    mu = dn*(vs**2.)

    # Horizontal slowness
    ns = _np.zeros((lnum, 1), dtype=CTP)
    for nl in range(lnum):
        ns[nl] = _np.cos(iS[nl])/vs[nl]

    # -------------------------------------------------------------------------
    # Data vector initialisation

    # Layer's amplitude vector
    A = _np.zeros((lnum*2, 1), dtype=CTP)

    # Input motion vector
    D = _np.zeros((lnum*2, 1), dtype=CTP)
    D[-1] = 1.

    # Layer matrix
    CORE = _np.zeros((lnum*2, lnum*2), dtype=CTP)

    # Ouput layer's amplitude matrix
    AMAT = _np.zeros(fnum, dtype=CTP)

    # -------------------------------------------------------------------------
    # Loop over frequencies

    for nf in range(fnum):

        # Reinitialise the layer matrix
        CORE *= 0.

        # Free surface constraints
        CORE[0, 0] = 1.
        CORE[0, 1] = -1.

        # Interface constraints
        for nl in range(lnum-1):

            row = (nl*2)+1
            col = nl*2

            expDSA = _np.exp(1j*angf[nf]*ns[nl]*hl[nl])
            expUSA = _np.exp(-1j*angf[nf]*ns[nl]*hl[nl])

            CORE[row, col+0] = expDSA[0]
            CORE[row, col+1] = expUSA[0]
            CORE[row, col+2] = -1.
            CORE[row, col+3] = -1.

            CORE[row+1, col+0] =  mu[nl][0]*ns[nl][0]*expDSA[0]
            CORE[row+1, col+1] = -mu[nl][0]*ns[nl][0]*expUSA[0]
            CORE[row+1, col+2] = -mu[nl+1][0]*ns[nl+1][0]
            CORE[row+1, col+3] =  mu[nl+1][0]*ns[nl+1][0]

        # Input motion constraints
        CORE[-1, -1] = 1.

        # Solving linear system of layer's displacement amplitudes
        try:
            A = _np.linalg.solve(CORE, D)
        except:
            A[:] = _np.nan

        # Adding displacement amplitude to the output matrix
        AMAT[nf] = A

    return AMAT

