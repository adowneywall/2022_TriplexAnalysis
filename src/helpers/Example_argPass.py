#!/usr/bin/env python3

import sys
import os

print('cmd entry:', sys.argv)

print(sys.argv[1])

print(os.listdir(sys.argv[1]))