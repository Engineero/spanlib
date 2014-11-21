import unittest
import numpy as npy
import os, sys
sys.path.insert(0, '../lib')
from spanlib.analyzer import Analyzer
from spanlib.filler import Filler
from spanlib_extra import setup_data1, setup_data2, gensin2d


class TSF(unittest.TestCase):

    def test_pca_rec(self):
        """Test a PCA analysis and reconstruction with a gap inside"""
#        data = setup_data1(nt=50, nx=120)
        data = gensin2d(xnper=5, ynper=5)
        data[20:70, 50:120] = npy.ma.masked
        A = Analyzer(data)
        rec = A.pca_rec()
#        import pylab as P
#        P.subplot(211)
#        P.pcolor(data, vmin=data.min(), vmax=data.max())
#        P.subplot(212)
#        P.pcolor(rec, vmin=data.min(), vmax=data.max())
#        P.show()
        self.assertAlmostEqual(((data-rec)**2).sum(), 48.2021433535)


    def test_fill_simple(self):
        """Test hole filling and forecast estimation with a single variable"""
        # Init
        ref = setup_data1(nt=50, nx=120, xyfact=0)
        withholes = ref.copy()
#        bad = npy.random.random_integers(0, ref.size-1, ref.size/20)
#        withholes.ravel()[bad] = npy.ma.masked
        withholes[25:35, 50:60] = npy.ma.masked
#        withholes[200:215] = npy.ma.masked
#        withholes[200:415] = npy.ma.masked
#        withholes[480:] = npy.ma.masked
#        withholes[380:] = npy.ma.masked

        # Fill
#        F = Filler(withholes, loglevel='debug')
#        filtered = F.filtered
#        import pylab as P
#        P.subplot(311)
#        P.pcolor(ref, vmin=ref.min(), vmax=ref.max())
#        P.subplot(312)
#        P.pcolor(withholes, vmin=ref.min(), vmax=ref.max())
#        P.subplot(313)
#        P.pcolor(filtered, vmin=ref.min(), vmax=ref.max())
#        P.show()
        self.assertTrue(npy.allclose(filled.filtered.filled()[200:202,0],
            npy.array([10.61948956,  10.74003157])))

    def test_fill_double(self):
        """Test hole filling and forecast estimation with a pair of variables"""
        # Init
        nt = 500
        tmax = 80.
        ref = npy.ma.sin(npy.linspace(0., tmax, nt))+10. # period = 2pi
        ref = npy.ma.resize(ref, (3, nt)).T
        ref[:, 1] = npy.ma.masked
        withholes = ref.copy()
        withholes[200:215] = npy.ma.masked
        withholes[480:] = npy.ma.masked
        # Fill
        filled = Filler([withholes, withholes*100])
#        import pylab as P
#        P.plot(filled.filtered[0][:, 0], 'r')
#        P.plot(withholes[:, 0], 'b')
#        P.show()
        self.assertTrue(npy.allclose(filled.filtered[0].filled()[200:202,0],
            npy.array([10.61834159,  10.73914392])))

if __name__ == '__main__':
    unittest.main()
