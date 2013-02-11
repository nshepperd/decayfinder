#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys

if len(sys.argv) != 3:
    # the program name and the two arguments
    # stop the program and print an error message
    sys.exit('usage: {} <output root filename> <number of events>'.format(sys.argv[0]))

import os, sys
from basf2 import *
from collections import namedtuple

rootFileName = sys.argv[1]
nOfEvents = int(sys.argv[2])
logFileName = rootFileName + '.log'

sys.stdout = open(logFileName, 'w')

# ---------------------------------------------------------------
# EvtGen
evtgeninput = register_module('EvtGenInput')
evtgeninput.param('boost2LAB', True)
# specify number of events to be generated in job
evtmetagen = register_module('EvtMetaGen')
evtmetagen.param('EvtNumList', [nOfEvents])  # process nOfEvents events
evtmetagen.param('RunList', [1])  # from run number 1
evtmetagen.param('ExpList', [1])  # and experiment number 1

progress = register_module('Progress')
gearbox = register_module('Gearbox')


# ---------------------------------------------------------------
# Add all modules to the main path
main = create_path()
main.add_module(evtmetagen)
main.add_module(evtgeninput)
main.add_module(progress)
main.add_module(gearbox)

# # c++ module attempt
# from descriptor import decode, to_json
# analysis = register_module('DecayFinder')
# analysis.param('pattern', to_json(decode('B0 => X- (e+ || mu+) (nu_e || nu_mu)')))
# main.add_module(analysis)

from decay_filter import DecayFilterModule
filter = DecayFilterModule()
# filter.set_decay('B0 -> (D*- -> D- pi0) e+ nu_e')
filter.set_decay('B0 --> X- [e+,mu+,tau+] [nu_e, nu_mu, nu_tau] ...')
main.add_module(filter)

process(main)
