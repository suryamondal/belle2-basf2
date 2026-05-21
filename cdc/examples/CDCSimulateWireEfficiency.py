#!/usr/bin/env python3

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

"""
Full-chain CDC wire-efficiency validation: simulate, reconstruct, and plot.

Runs e+e- -> B0 B0bar events through simulation and the default
reconstruction chain.  CDCDQMSim fills the wire-by-wire tracking efficiency
histogram (CDC/hTrackingWireEff) without requiring trigger or HLT simulation.
After processing, CDCWireEfficiencyCanvas.py reads the DQM file and writes
the 4-panel wire efficiency canvas.

Settings: exp=0, run=0, default global tag, default reconstruction chain.

Output:
    dqm_cdc.root             -- DQM histograms (CDCDQMSim output)
    wire_eff.pdf / .root     -- 4-panel wire efficiency canvas

Usage:
    basf2 CDCSimulateWireEfficiency.py [-- --events 200 --output wire_eff]
"""

import argparse
import os
import subprocess
import sys
import basf2 as b2
from simulation import add_simulation
from reconstruction import add_reconstruction
from generators import add_evtgen_generator

parser = argparse.ArgumentParser()
parser.add_argument('--events', type=int, default=200,
                    help='number of events to generate (default: 200)')
parser.add_argument('--dqm-output', default='dqm_cdc.root',
                    help='DQM ROOT file written by CDCDQMSim (default: dqm_cdc.root)')
parser.add_argument('--output', default='wire_eff',
                    help='canvas output base name without extension (default: wire_eff)')
parser.add_argument('--label', default='Belle II Simulation',
                    help='primary label drawn on the canvas')
parser.add_argument('--sublabel', default="B^{0}#bar{B^{0}} (no Bkg)",
                    help='secondary label drawn on the canvas')
args, _ = parser.parse_known_args()

# ── Phase 1: simulate, reconstruct, fill DQM histograms ──────────────────────

b2.set_log_level(b2.LogLevel.WARNING)
b2.set_random_seed(0)

main = b2.create_path()

main.add_module('EventInfoSetter', expList=[0], runList=[0], evtNumList=[args.events])
main.add_module('Progress')
main.add_module('HistoManager', histoFileName=args.dqm_output)
main.add_module('Gearbox')
main.add_module('Geometry')

add_evtgen_generator(path=main, finalstate='mixed')
add_simulation(main)
main.add_module('SimulateEventLevelTriggerTimeInfo')

add_reconstruction(main)

main.add_module('CDCDQMSim')

b2.process(main)
print(b2.statistics)

# ── Phase 2: produce wire efficiency canvas from DQM file ─────────────────────

script = os.path.join(os.path.dirname(__file__), 'CDCWireEfficiencyCanvas.py')
cmd = [
    'basf2', script, '--',
    '--input',    args.dqm_output,
    '--output',   args.output,
    '--label',    args.label,
    '--sublabel', args.sublabel,
]
result = subprocess.run(cmd, check=False)
sys.exit(result.returncode)
