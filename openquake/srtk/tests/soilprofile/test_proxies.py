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

from openquake.srtk.soilprofile import proxies


# =============================================================================

class DepthAverageTestCase(unittest.TestCase):

    def check_average(self,
                      thickness,
                      soil_param,
                      depth,
                      expected_result,
                      tolerance=0.):

        average = proxies.depth_weighted_average(thickness,
                                                 soil_param,
                                                 depth)
        self.assertAlmostEqual(average,
                               expected_result,
                               delta=tolerance)

    def test_one_layer(self):
        """
        Case with only one layer (homogenous half space)
        """

        self.check_average(np.array([0]),
                           np.array([100.]),
                           50.,
                           100.)

    def test_three_layers_case_1(self):
        """
        Case where the averaging depth is within the first layer
        """

        self.check_average(np.array([50., 10., 0.]),
                           np.array([5., 10., 50.]),
                           25.,
                           5.)

    def test_three_layers_case_2(self):
        """
        Case where the averaging depth is within an arbitrary layer
        """

        self.check_average(np.array([10., 20., 0.]),
                           np.array([5., 10., 50.]),
                           20.,
                           7.5)

    def test_three_layers_case_3(self):
        """
        Case where the averaging depth is below the last layer
        """

        self.check_average(np.array([10., 20., 0.]),
                           np.array([5., 10., 50.]),
                           100.,
                           37.5)

    def test_three_layers_case_4(self):
        """
        Case where the averaging depth is at an interface
        """

        self.check_average(np.array([10., 20., 0.]),
                           np.array([5., 10., 50.]),
                           30.,
                           8.3333,
                           tolerance=0.001)


# =============================================================================

class TravelTimeAverageTestCase(unittest.TestCase):

    def check_average(self,
                      thickness,
                      s_velocity,
                      depth,
                      expected_result,
                      tolerance=0.):

        average = proxies.traveltime_average_velocity(thickness,
                                                      s_velocity,
                                                      depth)
        self.assertAlmostEqual(average,
                               expected_result,
                               delta=tolerance)

    def test_one_layer(self):
        """
        Case with only one layer (homogenous half space)
        """

        self.check_average(np.array([0]),
                           np.array([200.]),
                           50.,
                           200.)

    def test_three_layer_vs30(self):
        """
        Case with only one layer (homogenous half space)
        """

        self.check_average(np.array([5., 10., 10., 20., 0.]),
                           np.array([200., 300., 350., 650., 1200.]),
                           30.,
                           317.13,
                           tolerance=0.01)


# =============================================================================

class SiteKappaTestCase(unittest.TestCase):

    def check_kappa0(self,
                     thickness,
                     s_velocity,
                     s_quality,
                     depth,
                     expected_result,
                     tolerance=0.):

        kappa0 = proxies.compute_site_kappa(thickness,
                                            s_velocity,
                                            s_quality,
                                            depth)
        self.assertAlmostEqual(kappa0,
                               expected_result,
                               delta=tolerance)

    def test_one_layer(self):
        """
        Case with only one layer (homogenous half space)
        """

        self.check_kappa0(np.array([0]),
                          np.array([100.]),
                          np.array([10.]),
                          50.,
                          0.05)

    def test_three_layers(self):
        """
        Case where the averaging depth is below the last layer
        """

        self.check_kappa0(np.array([10., 50., 0.]),
                          np.array([100., 500., 1000.]),
                          np.array([10., 20., 100.]),
                          100.,
                          0.01539,
                          tolerance=0.00001)
