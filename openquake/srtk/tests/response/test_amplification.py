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

import unittest
import numpy as np

from openquake.srtk.response.amplification import impedance_amplification


# =============================================================================

class ImpedanceAmplificationTestCase(unittest.TestCase):

    def check_amplification(self,
                            top_vs,
                            top_dn,
                            ref_vs,
                            ref_dn,
                            inc_angle,
                            expected_result,
                            tolerance=0.):

        computed_result = impedance_amplification(top_vs,
                                                  top_dn,
                                                  ref_vs,
                                                  ref_dn,
                                                  inc_angle)
        print computed_result

        if not isinstance(computed_result, np.ndarray):
            self.assertAlmostEquals(computed_result,
                                    expected_result,
                                    delta=tolerance)
        else:
            for cr, er, in zip(computed_result, expected_result):
                self.assertAlmostEquals(cr,
                                        er,
                                        delta=tolerance)

    def test_no_interface(self):
        """
        testing a model with no seismic contrast
        """

        self.check_amplification(200., 1900., 200., 1900., 10., 1.)

    def test_one_interface_vertical(self):
        """
        testing a model with a single seismic impedance contrast
        """

        self.check_amplification(200., 1900., 1500., 2500., 0., 3.14, 0.01)

    def test_one_interface_angle(self):
        """
        testing a model with a single seismic impedance contrast
        and non-vertical incidence angle
        """

        self.check_amplification(200., 1900., 1500., 2500., 45., 2.65, 0.01)

    def test_multiple_interfaces_vertical(self):
        """
        testing a model with an array of seismic impedance contrasts
        """

        self.check_amplification(np.array([200., 800., 2000.]),
                                 np.array([1900., 2100., 2500.]),
                                 2000.,
                                 2500.,
                                 0.,
                                 np.array([3.63, 1.72, 1.]),
                                 0.01)
