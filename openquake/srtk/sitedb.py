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
import openquake.srtk.response.amplification as _amp
import openquake.srtk.utils as _ut

# =============================================================================
# Constants & initialisation variables

# Precision for decimal rounding
DECIMALS = 6

GEO_KEYS = ['hl','vp','vs','dn','qp','qs']
ENG_KEYS = ['vsz','qwl','kappa','class']
AMP_KEYS = ['shtf','qwl','kappa']

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
            data can be a list of values sorted according to key list,
            or a dictionary with the corresponding keys

        :param int index:
            index of the position along the profile where the layer
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
            index of the position along the profile where the layer
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
            input model file. Default format is:
                hl,vp,vs,dn
                10,300,200,1900
                10,500,300,1900

        :param list header:
            list of header keys, to be used when not
            available within the input file

        :param int skip:
            number of intitial lines to skip;
            default value is 0

        :param char or string comment:
            string to mark comments (which are not parsed);
            default value is the hash character

        :param char delimiter:
            character separator between data fields;
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

        self.freq = []
        self.model = []
        self.mean = Model()

    #--------------------------------------------------------------------------

    def add_model(self, model=[], index=-1):
        """
        Add a single soil model to the site database

        :param Model model:
            the model the be added; if not specificed,
            an empty model is added

        :param int index:
            index of where to include the model in the database;
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
            single input model file or list of files

        :param list header:
            list of header keys, to be used when not
            available within the input file

        :param int skip:
            number of intitial lines to skip;
            default value is 0

        :param char or string comment:
            string to mark comments (which are not parsed);
            default value is the hash character

        :param char delimiter:
            character separator between data fields;
            default value is comma

        :param int index:
            index of where to include the model in the database;
            use -1 to append (default)

        :param boolean owrite:
            flag to enable model overwriting; in this case,
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

    def traveltime_velocity(self, depth=30.):
        """
        Compute and store travel-time average velocity at a given depth.
        Multiple depths are also allowed (as list). Default is Vs30.

        :param float or list depth:
            calculation depth
        """

        if not isinstance(depth, list):
            depth = [depth]

        for mod in self.model:
            mod.eng['vsz'] = {}
            for z in depth:
                vz = _avg.traveltime_velocity(mod.geo['hl'],
                                              mod.geo['vs'],
                                              depth=z)
                mod.eng['vsz'][z] = _ut.a_round(vz, DECIMALS)

        # Perform statistics (log-normal)
        self.mean.eng['vsz'] = {}
        for z in depth:
            data = [mod.eng['vsz'][z] for mod in self.model]
            mn, sd = _ut.log_stat(data)
            self.mean.eng['vsz'][z] = (_ut.a_round(mn, DECIMALS),
                                       _ut.a_round(sd, DECIMALS))

    #--------------------------------------------------------------------------

    def compute_soil_class(self, code='EC8'):
        """
        Compute geotechnical classification according to specified
        building code. Default is EC8 (missing special classes).

        :param string code:
            the reference building code for the classification
            (default EC8)
        """

        try:
            # Class from the mean Vs30 model
            vs30 = self.mean.eng['vsz'][30.][0]
            self.mean.eng['class'] = _avg.gt_soil_class(vs30, code)

            # Class of each profile
            for mod in self.model:
                vs30 = mod.eng['vsz'][30.]
                mod.eng['class'] = _avg.gt_soil_class(vs30, code)

        except:
            print 'Error: Vs30 must be calculated first'
            return

    #--------------------------------------------------------------------------

    def frequency_axis(self, fmin, fmax, fnum, log=True):
        """
        Compute a linear or logarithmic frequency axis.

        :param float fmin:
            minimum frequency

        :param float fmax:
            maximum frequency

        :param int fnum:
            number of frequencies

        :param boolean log:
            switch between linear or logarithmic spacing
            (default is logarithmic)
        """

        self.freq = _amp.frequency_axis(fmin, fmax, fnum, log)


    def _check_frequency(self):
        """
        Internal: check if the frequency axis has been instantiated
        """

        if not _np.sum(self.freq):
            raise ValueError('Frequency axis must be first instantiated')

    #--------------------------------------------------------------------------

    def quarter_wavelength_average(self):
        """
        Compute quarter-wavelength parameters (velocity and density)
        and store them into the site database.
        """

        self._check_frequency()

        for mod in self.model:

            # Compute average velocity
            qwl_par = _avg.quarter_wavelength_average(mod.geo['hl'],
                                                      mod.geo['vs'],
                                                      mod.geo['dn'],
                                                      self.freq)

            mod.eng['qwl'] = {}
            mod.eng['qwl']['z'] = _ut.a_round(qwl_par[0], DECIMALS)
            mod.eng['qwl']['vs'] = _ut.a_round(qwl_par[1], DECIMALS)
            mod.eng['qwl']['dn'] = _ut.a_round(qwl_par[2], DECIMALS)

        # Perform statistics (log-normal)
        self.mean.eng['qwl'] = {}
        for key in ['z','vs','dn']:
            data = [mod.eng['qwl'][key] for mod in self.model]
            mn, sd = _ut.log_stat(data)
            self.mean.eng['qwl'][key] = (_ut.a_round(mn, DECIMALS),
                                         _ut.a_round(sd, DECIMALS))


    #--------------------------------------------------------------------------

    def quarter_wavelength_amplification (self, vs_ref=[],
                                                dn_ref=[],
                                                inc_ang=0.):
        """
        Compute the amplification as impedance contrast of
        quarter-wavelength parameters. Aribitrary reference
        can be provided, otherwise the last layer of the model
        is used. Angle of incidence is optional.

        :param float or numpy.array ref_vs:
            lowermost (reference) shear-wave velocity in m/s

        :param float or numpy.array ref_dn:
            lowermost (reference) density in kg/m3

        :param float inc_ang:
            angle of incidence in degrees, relative to the vertical
            (default is vertical incidence)
        """

        for mod in self.model:

            if _ut.is_empty(vs_ref):
                vs_ref = mod.geo['vs'][-1]

            if _ut.is_empty(dn_ref):
                dn_ref = mod.geo['dn'][-1]

            qwl_amp = _amp.impedance_amplification(mod.eng['qwl']['vs'],
                                                   mod.eng['qwl']['dn'],
                                                   vs_ref,
                                                   dn_ref,
                                                   inc_ang)

            mod.amp['qwl'] = _ut.a_round(qwl_amp, DECIMALS)

        # Perform statistics (log-normal)
        data = [mod.amp['qwl'] for mod in self.model]
        mn, sd = _ut.log_stat(data)
        self.mean.amp['qwl'] = (_ut.a_round(mn, DECIMALS),
                                _ut.a_round(sd, DECIMALS))


    #--------------------------------------------------------------------------

    def compute_site_kappa(self, depth=[]):
        """
        Compute the Kappa parameter directly from the site model
        using harmonic averaging. Caculation depth can be speficied,
        otherwise the depth of the last layer interface is used.

        :param float depth:
            averaging depth in meters (optional)
        """

        for mod in self.model:

            # Compute kappa attenuation
            kappa = _avg.compute_site_kappa(mod.geo['hl'],
                                            mod.geo['vs'],
                                            mod.geo['qs'],
                                           depth)

            mod.eng['kappa'] = _ut.a_round(kappa, DECIMALS)

        # Perform statistics (normal)
        data = [mod.eng['kappa'] for mod in self.model]
        mn, sd = _ut.lin_stat(data)
        self.mean.eng['kappa'] = (_ut.a_round(mn, DECIMALS),
                                  _ut.a_round(sd, DECIMALS))


    #--------------------------------------------------------------------------

    def attenuation_decay(self):
        """
        Compute the frequency-dependent attenuation function
        for a given site Kappa (0).
        """

        self._check_frequency()

        for mod in self.model:

            # Compute attenuation decay
            att_fun = _amp.attenuation_decay(self.freq, mod.eng['kappa'])

            mod.amp['kappa'] = _ut.a_round(att_fun, DECIMALS)

        # Perform statistics (log-normal)
        data = [mod.amp['kappa'] for mod in self.model]
        mn, sd = _ut.log_stat(data)
        self.mean.amp['kappa'] = (_ut.a_round(mn, DECIMALS),
                                  _ut.a_round(sd, DECIMALS))


    #--------------------------------------------------------------------------

    def sh_transfer_function(self, inc_ang=0., elastic=False, complex=False):
        """
        Compute the complex SH-wave transfer function at the
        surface for outcropping rock reference conditions.
        Calculation can be done for an arbitrary angle of
        incidence (0-90), elastic or anelastic.

        For more options (e.g. calculation at arbitrary depth)
        see the generic sh_transfer_function method in the
        amplification module

        Statistic is performed linearly on complex spectra
        (to check!)

        :param float inc_ang:
            angle of incidence in degrees, relative to the
            vertical (default is vertical incidence)

        :param boolean elastic:
            switch between elastic and anelastic calculation
            (default is anelastic)

        :param boolean complex:
            switch to output real (abs) or complex spectra
            (note that type of statistic is affected)
        """

        self._check_frequency()

        for mod in self.model:

            qs = mod.geo['qs'] if not elastic else None

            # Compute transfer function
            dis_mat = _amp.sh_transfer_function(self.freq,
                                                mod.geo['hl'],
                                                mod.geo['vs'],
                                                mod.geo['dn'],
                                                qs,
                                                inc_ang,
                                                0)

            if complex:
                mod.amp['shtf'] = dis_mat[0]/2
            else:
                mod.amp['shtf'] = _np.abs(dis_mat[0])/2

        # Perform statistics (normal on complex)
        data = [mod.amp['shtf'] for mod in self.model]
        if complex:
            self.mean.amp['shtf'] = _ut.lin_stat(data)
        else:
            self.mean.amp['shtf'] = _ut.log_stat(data)


    #--------------------------------------------------------------------------

    def resonance_frequency(self):
        """
        Identify resonance frequencies on an amplification spectrum.
        Note that frequency on the average spectrum are identified
        directly without performing statistic on the single models.
        """

        for mod in self.model:
            mod.amp['fn'] = _amp.resonance_frequency(self.freq,
                                                     mod.amp['shtf'])

        fn  = _amp.resonance_frequency(self.freq,
                                       self.mean.amp['shtf'][0])
        self.mean.amp['fn'] = fn

