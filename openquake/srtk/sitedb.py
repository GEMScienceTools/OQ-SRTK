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
import openquake.srtk.soil.average as _avg
import openquake.srtk.utils as _ut

# =============================================================================
# Constants

GEO_KEYS = ['hl','vp','vs','dn','qp','qs']
ENG_KEYS = ['vsz','qwl','kappa','class']
AMP_KEYS = ['shtf','aimp','attf']

# Precision for decimal rounding
DECIMALS = 4

# =============================================================================

class Model(object):
    """
    Base class to store a single site model, including the vertical
    soil profile, derived engineering parameters and amplification.
    """

    #--------------------------------------------------------------------------

    def __init__(self):

        self._geo_init()
        self._eng_init()
        self._amp_init()

    def _geo_init(self):
        self.geo = {}
        for K in GEO_KEYS: self.geo[K] = _np.array([])

    def _eng_init(self):
        self.eng = {}
        for K in ENG_KEYS: self.eng[K] = _np.array([])

    def _amp_init(self):
        self.amp = {}
        for K in AMP_KEYS: self.amp[K] = _np.array([])

    #--------------------------------------------------------------------------

    def add_layer(self, data, index=-1):
        """
        Method to add a single layer to the 1d soil profile
        at arbitrary location.

        :param list or dictionary data:
            Data can be a list of values sorted according to key list,
            or a dictionary with the corresponding keys

        :param int index:
            Index of the position along the profile where the layer
            should be added. Use -1 for the last layer (default).
        """

        index = int(index)

        # Case: List
        if isinstance(data, list):
            for I, K in enumerate(GEO_KEYS):
                if index < 0: index = len(self.geo[K])
                if I < len(data):
                    self.geo[K] = _np.insert(self.geo[K], index, data[I])
                else:
                    self.geo[K] = _np.insert(self.geo[K], index, _np.nan)

        # Case: Dictionary
        if isinstance(data, dict):
            for K in GEO_KEYS:
                if index < 0: index = len(self.geo[K])
                if K in data.keys():
                    self.geo[K] = _np.insert(self.geo[K], index, data[K])
                else:
                    self.geo[K] = _np.insert(self.geo[K], index, _np.nan)

    #--------------------------------------------------------------------------

    def del_layer(self, index=-1):
        """
        Method to remove a single layer from te 1d soil profile
        ar arbitrary location

        :param int index:
            Index of the position along the profile where the layer
            should be removed. Use -1 for the last layer (default).
        """

        for K in GEO_KEYS:
            self.geo[K] = _np.delete(self.geo[K], int(index))

    #--------------------------------------------------------------------------
    def from_file(self, ascii_file, header=[],
                                    skip=0,
                                    comment='#',
                                    delimiter=','):
        """
        Method to parse soil properties from tabular ascii file;
        arbitrary formatting is allowed

        :param string ascii_file:
            Input model file. Default format is:
                hl,vp,vs,dn
                10,300,200,1900
                10,500,300,1900

        :param list header:
            List of header keys, to be used when not
            available within the input file

        :param int skip:
            Number of intitial lines to skip;
            default value is 0

        :param char or string comment:
            String to mark comments (which are not parsed);
            default value is the hash character

        :param char delimiter:
            Character separator between data fields;
            default value is comma
        """

        # Delete any previous model
        self._geo_init()

        # Open input ascii file
        try:
            f = open(ascii_file, 'r')

        except:
            print('Error: Wrong file or file path')
            return

        else:
            # Ignore initial line(s) if necessary
            for _ in range(0, skip):
                next(f)

            for line in f:
                line = line.strip()

                # Skip comments
                if line[0] != comment:
                    data = line.split(delimiter)

                    # Import header and data
                    if not header:
                        header = data
                    else:
                        layer = {k:float(d) for k,d in zip(header,data)}
                        self.add_layer(layer)

            f.close()


# =============================================================================

class Site1D(object):
    """
    Base class for a single one-dimensional site.
    It contains the collection of soil models and the
    corresponding derived parameters.
    """

    def __init__(self, id=None, x=None, y=None, z=None):

        self.head = {}
        self.head['id'] = id
        self.head['x'] = x
        self.head['y'] = y
        self.head['z'] = x

        self.model = []

        self._eng_init()

    def _eng_init(self):
        self.eng = {}
        for K in ENG_KEYS: self.eng[K] = _np.array([])

    #--------------------------------------------------------------------------

    def add_model(self, model=[], index=-1):
        """
        Add a single soil model to the site database

        :param Model model:
            The model the be added; if not specificed,
            an empty model is added

        :param int index:
            Index of where to include the model in the database;
            use -1 to append (default)
        """

        index = int(index)
        if index < 0: index = len(self.model)

        if model:
            self.model.insert(index, model)
        else:
            self.model.insert(index, Model())

    #--------------------------------------------------------------------------

    def del_model(self, index=-1):
        """
        Remove a model from the site database

        :param int index:
            Index of the model to be removed from the database;
            use -1 for the last position (default)
        """

        del self.model[int(index)]

    #--------------------------------------------------------------------------

    def read_model(self, ascii_file, header=[],
                                     skip=0,
                                     comment='#',
                                     delimiter=',',
                                     index=-1,
                                     owrite=False):
        """
        Method to parse soil properties from a single tabular
        ascii file or a list of files; arbitrary formatting is allowed
        
        The method is essentially a wrapper of the read_model
        method of the Model() class, from whom it inherits the
        input paramters (header, skip, ...)

        :param string or list ascii_file:
            Single input model file or list of files

        :param list header:
            List of header keys, to be used when not
            available within the input file

        :param int skip:
            Number of intitial lines to skip;
            default value is 0

        :param char or string comment:
            String to mark comments (which are not parsed);
            default value is the hash character

        :param char delimiter:
            Character separator between data fields;
            default value is comma

        :param int index:
            Index of where to include the model in the database;
            use -1 to append (default)

        :param boolean owrite:
            Flag to enable model overwriting; in this case,
            the index is that of the model to be overwritten
        """

        if not isinstance(ascii_file, list):
            ascii_file = [ascii_file]

        for af in ascii_file:
            model = Model()
            model.from_file(af, header, skip, comment, delimiter)

            if owrite:
                self.model[index] = model
            else:
                self.add_model(model, index)

    #--------------------------------------------------------------------------

    def traveltime_average_velocity(self, depth=30.):
        """
        Compute and store travel-time average velocity at a given depth.
        Multiple depths are also allowed (as list). Default is Vs30.

        :param float or list depth:
            Calculation depth
        """

        if not isinstance(depth, list):
            depth = [depth]

        for mod in self.model:
            mod.eng['vsz'] = {}
            for z in depth:
                vz = _avg.traveltime_average_velocity(mod.geo['hl'],
                                                      mod.geo['vs'],
                                                      depth=z)
                mod.eng['vsz'][z] = _ut.a_round(vz, DECIMALS)

        # Perform statistics
        self.eng['vsz'] = {}
        for z in depth:
            data = [mod.eng['vsz'][z] for mod in self.model]
            mn, sd = _ut.log_stat(data)
            self.eng['vsz'][z] = (_ut.a_round(mn, DECIMALS),
                                  _ut.a_round(sd, DECIMALS))


