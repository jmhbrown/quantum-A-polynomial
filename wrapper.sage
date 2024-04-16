#!/usr/bin/env python
# coding: utf-8

import argparse
import logging
import json

parser = argparse.ArgumentParser(description='Compute (quantum!) A-polynomials')
parser.add_argument('-f', '--file', type=argparse.FileType('r'), help='JSON file with run options.')
parser.add_argument('-o', '--output', help='File where output is written.')
parser.add_argument('-l', '--log_level', default='info', choices=['debug','info','warning','error','critical'],help='Set the logging level.')

args = parser.parse_args()

options = json.load(args.file)


# Logging

if args.output is None: # log to the console.
    logging.basicConfig(level=args.log_level.upper(),format='%(levelname)s - %(message)s')
else:
    logging.basicConfig(level=args.log_level.upper(),format='%(levelname)s - %(message)s',filename=args.output)

logger = logging.getLogger(__name__)

sage.repl.load.load('main.sage',globals())

for run in options["runs"]:
    logger.info("Running: {}".format(run))

    t = QuantumAPolynomial(run['knot'],run['pachner_moves'])

    t.compute_skein_module()

