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
Module containing the database classes to handle site information.
"""

import numpy as _np


# =============================================================================

class Model(object):
    """
    Base class to store a single site model, including the vertical
    soil profile, derived engineering parameters and amplification.
    """

    geo_keys = ['hl','vp','vs','dn','qp','qs']
    eng_keys = ['vz','qwl','k0','gc']
    amp_keys = ['shtf','imp','attf']

    #--------------------------------------------------------------------------

    def __init__(self):
        self._par_init()
        self._eng_init()
        self._amp_init()

    def _par_init(self):
        self.geo = {}
        for K in self.geo_keys: self.geo[K] = []

    def _eng_init(self):
        self.eng = {}
        for K in self.eng_keys: self.eng[K] = []

    def _amp_init(self):
        self.amp = {}
        for K in self.amp_keys: self.amp[K] = []

    #--------------------------------------------------------------------------

    def add_layer(self, data, index=-1):
        """
        Method to add a single layer to the 1d soil profile
        at arbitrary location.

        :param listo or dictionary data:
            Data can be a list of values sorted according to key list,
            or a dictionary with the corresponding keys
        :param int index:
            Index of the position along the profile where the layer
            should be added. Use -1 for the last layer (default).
        """

        index = int(index)

        # Case: List
        if isinstance(data, list):
            for I, K in enumerate(self.geo_keys):
                if index < 0: index = len(self.geo[K])
                if I < len(data):
                    self.geo[K].insert(index, data[I])
                else:
                    self.geo[K].insert(index, None)

        # Case: Dictionary
        if isinstance(data, dict):
            for K in self.geo_keys:
                if index < 0: index = len(self.geo[K])
                if K in data.keys():
                    self.geo[K].insert(index, data[K])
                else:
                    self.geo[K].insert(index, None)

    #--------------------------------------------------------------------------

    def del_layer(self, index=-1):
        """
        Method to remove a single layer from te 1d soil profile
        ar arbitrary location

        :param int index:
            Index of the position along the profile where the layer
            should be removed. Use -1 for the last layer (default).
        """

        del self.geo[int(index)]


# =============================================================================

class Site1D(object):
    """
    Base class for a single one-dimensional site.
    It contains the collection of soil models and the
    corresponding derived parameters.
    """

    def __init__(self, id=None, x=None, y=None, z=None):

      self.head = {}
      self.head['id'] = Id
      self.head['x'] = X
      self.head['y'] = Y
      self.head['z'] = Z

      self.model = []

    #--------------------------------------------------------------------------

    def add_model(self, model=[], index=-1):
      """
      Add a soil model to the site database
      """

      index = int(index)
      if index < 0: index = len(self.model)

      if mod:
        self.model.insert(index, model)
      else:
        self.model.insert(index, Model())

    #--------------------------------------------------------------------------

    def del_model(self, index=-1):
      """
      Remove a model from the site database.
      """

      del self.model[int(index)]

