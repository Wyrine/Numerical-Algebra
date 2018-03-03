#!/usr/local/bin/python3

import sys
import numpy as np


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: ./conditionNumber.py hilbertDimensions")
        sys.exit(1)
    n = int(sys.argv[1])
