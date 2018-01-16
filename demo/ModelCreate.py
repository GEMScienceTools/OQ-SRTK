import numpy as np
import matplotlib.pyplot as plt

from openquake.srtk import sitedb

#--------------------------------------------------------------
# Manually creating a site model using dictionaries


# Case 1
mod = sitedb.Model()

mod.add_layer({'hl': 10.,
               'vp': 300.,
               'vs': 200.,
               'dn': 1900.})

mod.add_layer({'hl': 0.,
               'vp': 1000.,
               'vs': 800.,
               'dn': 2100.})

# Case 2
mod = sitedb.Model()

mod.add_layer({'hl': 10.,
               'vp': 300.,
               'vs': 200.,
               'dn': 1900.,
               'qp': 50.,
               'qs': 20.})

mod.add_layer({'hl': 10.,
               'vp': 500.,
               'vs': 300.,
               'dn': 1900.,
               'qp': 50.,
               'qs': 20.})

mod.add_layer({'hl': 0.,
               'vp': 1000.,
               'vs': 800.,
               'dn': 2100.,
               'qp': 100.,
               'qs': 50.})

#--------------------------------------------------------------
# Manually creating a site model using lists

# Case 1
mod = sitedb.Model()
mod.add_layer([10., 300., 200., 1900.])
mod.add_layer([0., 1000., 800., 2100.])

# Case 2

mod = sitedb.Model()
mod.add_layer([10., 300., 200., 1900., 50., 20.])
mod.add_layer([10., 500., 300., 1900., 50., 20.])
mod.add_layer([0., 1000., 800., 2100., 100., 50.])

#--------------------------------------------------------------
# Import site models

site = sitedb.Site1D(0,0,0)
site.read_model(['data/site01.csv', 'data/site02.csv'])
