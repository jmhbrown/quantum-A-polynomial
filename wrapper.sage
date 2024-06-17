#!/usr/bin/env python
# coding: utf-8

import argparse
import logging
import json

parser = argparse.ArgumentParser(description='Compute (quantum!) A-polynomials')
parser.add_argument('-f', '--file', type=argparse.FileType('r'), help='JSON file with run options.')
parser.add_argument('-o', '--output', help='File where output is written.')
parser.add_argument('-l', '--log_level', default='info', choices=['debug','info','warning','error','critical'],help='Set the logging level.')
parser.add_argument('-t', '--timeout', default='30', type=int, help='Set the timeout duration in minutes.')
parser.add_argument('-r', '--randomize', default=0, type=int, help='Number of times we randomize the triangulation then re-compute.')

args = parser.parse_args()

options = json.load(args.file)

# Logging

if args.output is None: # log to the console.
    logging.basicConfig(level=args.log_level.upper(),format='%(asctime)s %(levelname)s - %(message)s')
else:
    logging.basicConfig(level=args.log_level.upper(),format='%(asctime)s %(levelname)s - %(message)s',filename=args.output)

logger = logging.getLogger(__name__)

sage.repl.load.load('main.sage',globals())

for run in options["runs"]:
    logger.info("\n\n")
    logger.info("---------------------------------\n Running: {0}".format(run))

    alarm(60*(args.timeout)) # set timeout length.


    t = QuantumAPolynomial(run['knot'],run['pachner_moves'])
    r = args.randomize
    while r > -1:
        try:
            t.compute_skein_module()
            if r > 0: # set up for the next run
                t.knot_comp.randomize()
                logger.info("---------------------------------\n Randomizing then re-running {}".format(run))
                logger.debug("Randomized triangulation: {}".format(QuantumAPolynomial.get_gluing_dict(t.knot_comp)))
            r += -1
        except AlarmInterrupt:
            r += -1
            continue
        else:
            # reset the alarm
            cancel_alarm()

