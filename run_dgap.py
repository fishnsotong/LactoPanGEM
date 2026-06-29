#!/usr/bin/env python3
"""Compatibility wrapper for the maintained dgap.py completion CLI.

The original run_dgap.py was a one-off driver hard-coded to /home/omidard
paths and one template model. Keep this entry point for older notes, but route
all real work through dgap.py.
"""

from dgap import main


if __name__ == "__main__":
    main()
