"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.core import stats as pcs

import numpy as np

import unittest

sqx = [float(n)/20.0 for n in range(21)]
sqy = [np.sqrt(float(n)/20.0) for n in range(21)]
areasq = 3.0/8.0

#TO DO
class TestPISACovStats(unittest.TestCase):
    def test_bezier_fit_1(self):
        
        sqdata = [sqx, sqy]
        
        t, B, derB, der2B, LR, K, AUC = pcs.bezier_parametrization(sqdata)
        
        
        expected_subint = [[5, 5]]

        calc_subint = cei._intervalise(_INTEGER_1).subint

        self.assertListEqual(calc_subint, expected_subint)