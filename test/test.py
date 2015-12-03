#!/usr/bin/env python

import sys
sys.path.append("..")

import unittest
import math
import calculate_BDT

"""
              W    
              .   R
              .  .
              . .
        1ang  ..
    K . . . . D . . . . V
             ..
            . .
           .  .
          T   .
              C


    Octahedron like this.

"""

pdb_file      = "./data/pdb1a8e_onlyCA.ent"
bindres_files = {
    'D'   : "./data/pdb1a8e_onlyCA_D.bindres",
    'K'   : "./data/pdb1a8e_onlyCA_K.bindres",
    'DK'  : "./data/pdb1a8e_onlyCA_DK.bindres",
    'DW'  : "./data/pdb1a8e_onlyCA_DW.bindres",
    'all' : "./data/pdb1a8e_onlyCA_all.bindres"
}
pred_files    = {
    'D'   : "./data/pdb1a8e_onlyCA_D.pred",
    'K'   : "./data/pdb1a8e_onlyCA_K.pred",
    'T'   : "./data/pdb1a8e_onlyCA_T.pred",
    'V'   : "./data/pdb1a8e_onlyCA_V.pred",
    'R'   : "./data/pdb1a8e_onlyCA_R.pred",
    'W'   : "./data/pdb1a8e_onlyCA_W.pred",
    'C'   : "./data/pdb1a8e_onlyCA_C.pred",
    'DK'  : "./data/pdb1a8e_onlyCA_KD.pred",
    'WC'  : "./data/pdb1a8e_onlyCA_WC.pred",
    'all' : "./data/pdb1a8e_onlyCA_all.pred"
}


def score(d, d_zero):
    return 1.0 / (1 + (d / d_zero)**2)

def my_round(f):
    n = 6
    return int(f * (10**n))


d_zeros = [1.0, 2.0, 3.0, 5.0]
aminos = ['D', 'K', 'T', 'V', 'R', 'W', 'C']
class TestCalculate_BDT(unittest.TestCase):

    def test_calculate_BDT_allright(self):
        for i in xrange(len(d_zeros)):
            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['D'], pred_files['D'], d_zeros[i])
            self.assertEqual(BDT, 1.0)
            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['K'], pred_files['K'], d_zeros[i])
            self.assertEqual(BDT, 1.0)
            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['DK'], pred_files['DK'], d_zeros[i])
            self.assertEqual(BDT, 1.0)
            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['all'], pred_files['all'], d_zeros[i])
            self.assertEqual(BDT, 1.0)

    def test_calculate_BDT_singlepred(self):
        for i in xrange(len(d_zeros)):
            for amino in aminos:
                if amino == 'D':
                    continue
                BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['D'], pred_files[amino], d_zeros[i])
                correct_BDT = score(1, d_zeros[i])
                self.assertEqual(BDT, correct_BDT)

    def test_calculate_BDT_doublepred(self):
        for i in xrange(len(d_zeros)):
            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['D'], pred_files['DK'], d_zeros[i])
            N = 2.0
            correct_BDT = (1 + score(1, d_zeros[i])) / N
            self.assertEqual(BDT, correct_BDT)

            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['D'], pred_files['WC'], d_zeros[i])
            N = 2.0
            correct_BDT = 2 * score(1, d_zeros[i]) / N
            self.assertEqual(BDT, correct_BDT)

    def test_calculate_BDT_doublecorrect(self):
        for i in xrange(len(d_zeros)):
            for amino in aminos:
                BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['DK'], pred_files[amino], d_zeros[i])
                N = 2.0
                if amino in {'D', 'K'}:
                    correct_BDT = 1 / N
                else:
                    correct_BDT = score(1, d_zeros[i]) / N
                self.assertEqual(BDT, correct_BDT)

                BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['DW'], pred_files[amino], d_zeros[i])
                if amino in {'D', 'W'}:
                    correct_BDT = 1 / N
                else:
                    correct_BDT = score(1, d_zeros[i]) / N
                self.assertEqual(BDT, correct_BDT)

    def test_calculate_BDT_allpred(self):
        for i in xrange(len(d_zeros)):
            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['K'], pred_files['all'], d_zeros[i])
            N = 7.0
            correct_BDT = (1 + score(1, d_zeros[i]) + 4 * score(math.sqrt(2), d_zeros[i]) + score(2, d_zeros[i])) / N
            self.assertEqual(my_round(BDT), my_round(correct_BDT))

            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['D'], pred_files['all'], d_zeros[i])
            correct_BDT = (1 + 6 * score(1, d_zeros[i])) / N
            self.assertEqual(my_round(BDT), my_round(correct_BDT))

            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['DK'], pred_files['all'], d_zeros[i])
            correct_BDT = (2 * 1 + 5 * score(1, d_zeros[i])) / N
            self.assertEqual(my_round(BDT), my_round(correct_BDT))

            BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['DW'], pred_files['all'], d_zeros[i])
            correct_BDT = (2 * 1 + 5 * score(1, d_zeros[i])) / N
            self.assertEqual(my_round(BDT), my_round(correct_BDT))

    def test_calculate_BDT_allcorrect(self):
        for i in xrange(len(d_zeros)):
            for pred_key in pred_files.keys():
                BDT = calculate_BDT.calculate_BDT(pdb_file, bindres_files['all'], pred_files[pred_key], d_zeros[i])
                N = 7.0
                if len(pred_key) == 1:
                    correct_BDT = 1 / N
                elif len(pred_key) == 2:
                    correct_BDT = 2 / N
                elif pred_key == 'all':
                    continue
                self.assertEqual(BDT, correct_BDT)
        
if __name__ == "__main__":
    unittest.main()
