#!/usr/bin/env python

import argparse

from rna_gpf.base_pairing import BasePairing
from rna_gpf.mapping_functions import gp_mapper, nussinov


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--subopt_map", help="Suboptimal sets")
    parser.add_argument("-p", "--prop", help="mfe propensity scores")

    

    args = parser.parse_args()
    
