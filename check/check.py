#!/usr/bin/env python2

import sys

# Check for python version - requires 2.7 (and not 3)
if sys.version_info[:2] != (2,7):
    print "Requires Python 2.7"
    print "On BCP3 run:"
    print "module load languages/python-2.7.6"
    exit(1)

import argparse
import numpy as np


# Intermediate class to parse arguments
class InputParser(argparse.ArgumentParser):
    def __init__(self):
        super(InputParser, self).__init__(
            description="Testing script for HPC LBM coursework",
            fromfile_prefix_chars='@',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            )

        # % tolerance
        self.add_argument("--tolerance",
            nargs=1,
            default=[1],
            type=float,
            help="""Percentage tolerance to match against reference results""",
            action='store')

        # Reference results
        self.add_argument("--ref-av-vels-file",
            nargs=1,
            required=True,
            help="""reference av_vels results file""",
            action='store')

        self.add_argument("--ref-final-state-file",
            nargs=1,
            required=True,
            help="""reference final_state results file""",
            action='store')

        # Calculated results
        self.add_argument("--av-vels-file",
            nargs=1,
            required=True,
            help="""calculated av_vels results file""",
            action='store')

        self.add_argument("--final-state-file",
            nargs=1,
            required=True,
            help="""calculated final_state results file""",
            action='store')

parser = InputParser()
parsed_args = parser.parse_args()

def load_dat_files(av_vels_filename, final_state_filename):
    with open(av_vels_filename, "r") as av_vels_ref_file:
        with open(final_state_filename, "r") as final_state_ref_file:
            av_vels = np.loadtxt(av_vels_ref_file, usecols=[1])
            final_state = np.loadtxt(final_state_ref_file, usecols=[0, 1, 5])

            return av_vels, final_state

# Open reference and input files
av_vels_ref, final_state_ref = load_dat_files(parsed_args.ref_av_vels_file[0], parsed_args.ref_final_state_file[0])
av_vels_sim, final_state_sim = load_dat_files(parsed_args.av_vels_file[0], parsed_args.final_state_file[0])

# Make sure the coordinates are in the right order
if np.any(final_state_ref[:,0:2] != final_state_sim[:,0:2]):
    print "Final state files coordinates were not the same"
    exit(1)

# Make sure the av_vels have the same number of steps
if av_vels_ref.size != av_vels_sim.size:
    print "Different number of steps in av_vels files"
    exit(1)

def get_diff_values(ref_vals, sim_vals):
    # Get the differences between the original and reference results
    diff = ref_vals - sim_vals
    diff_pcnt = 100.0*(diff/(ref_vals - diff))

    max_diff_step = np.argmax(np.abs(diff_pcnt))

    diffs = {
        "max_diff_step" : max_diff_step,
        "max_diff" : diff[max_diff_step],
        "max_diff_pcnt" : diff_pcnt[max_diff_step],
        "sim_val" : sim_vals[max_diff_step],
        "ref_val" : ref_vals[max_diff_step],
        "total" : np.sum(np.abs(diff)),
    }

    return diffs

def print_diffs(format_strings, format_dict):
    for s in format_strings:
        print s.format(**format_dict)

av_vels_diffs = get_diff_values(av_vels_ref, av_vels_sim)
av_vels_strings = [
    "Total difference in av_vels : {total:.12E}",
    "Biggest difference (at step {max_diff_step:d}) : {max_diff:.12E}",
    "  {sim_val:.12E} vs. {ref_val:.12E} = {max_diff_pcnt:.2g}%",
]

print_diffs(av_vels_strings, av_vels_diffs)

print

final_state_diffs = get_diff_values(final_state_ref[:,2], final_state_sim[:,2])
final_state_strings = [
    "Total difference in final_state : {total:.12E}",
    "Biggest difference (at coord ({jj:d},{ii:d})) : {max_diff:.12E}",
    av_vels_strings[2],
]

# We want the location of the biggest difference
max_diff_loc = int(final_state_diffs["max_diff_step"])
final_state_diffs["jj"] = int(final_state_sim[max_diff_loc,0])
final_state_diffs["ii"] = int(final_state_sim[max_diff_loc,1])

print_diffs(final_state_strings, final_state_diffs)

print

# Find out if either of them failed
final_state_failed = (not np.isfinite(final_state_diffs["max_diff_pcnt"])) or (np.abs(final_state_diffs["max_diff_pcnt"]) > parsed_args.tolerance[0])
av_vels_failed = (not np.isfinite(av_vels_diffs["max_diff_pcnt"])) or (np.abs(av_vels_diffs["max_diff_pcnt"]) > parsed_args.tolerance[0])

if final_state_failed:
    print "final state failed check"
if av_vels_failed:
    print "av_vels failed check"

# Return 1 on failure
if final_state_failed or av_vels_failed:
    exit(1)
else:
    print "Both tests passed!"
    exit(0)

