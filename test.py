#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 23:51:17 2019

@author: david
"""

import regelungstechnik as rt
import numpy as np

T = 1

I = rt.I(T=T)
PT2 = rt.PT2(omega=2 * np.pi, D=0)

step = rt.StepResponse([I, PT2], ["$I$", "$PT_2$"], duration=4, N=1024)
step.save(filename="test.png")