#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 23:51:17 2019

@author: david
"""

import regelungstechnik as rt

I = rt.I(T=0.5)

step = rt.StepResponse([I], ["I"], duration=3)
step.save(filename="test.png")