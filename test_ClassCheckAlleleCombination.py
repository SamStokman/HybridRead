"""
30-07-'19

This script contains unittests for the class CheckAlleleCombination from SelectHybridReads.py.


"""

import unittest
import SelectHybridReads


class TestCheckAlleleCombination(unittest.TestCase):
    """
    This class contains unittests

    """

    def setUp(self):
        """

        """
        self.read_consensus = '---CC**CCC---'

        self.indicator_string = ''

        self.allele_data = [['allele_A1','---CC--CCC---'],
                             ['allele_A2','---CC--TCC---'],
                             ['allele_B1','---CC--CCC---'],
                             ['allele_B2','---CC--CTT---'],
                             ['allele_C1','---CC--CGC---'],
                             ['allele_C2','---CC--CCC---']]

    def test_create_indicator_string(self):
        """
        Get the aligened sequences of the given alleles.

        """

        # Test case 1, second allele has a mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        self.assertEqual(Allele_test.create_indicator_string(), ['---CC--CCC---', '---CC--TCC---'])
        self.assertEqual(Allele_test.indicator_string, '-------Y-----')

        # Test case 2, first allele has a mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A2', 'allele_B1'], self.allele_data)
        self.assertEqual(Allele_test.create_indicator_string(), ['---CC--TCC---', '---CC--CCC---'])
        self.assertEqual(Allele_test.indicator_string, '-------X-----')

        # Test case 2, both alleles have a mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A2', 'allele_B2'], self.allele_data)
        self.assertEqual(Allele_test.create_indicator_string(), ['---CC--TCC---', '---CC--CTT---'])
        self.assertEqual(Allele_test.indicator_string, '-------XYY---')

        # Test case 3, alleles have mutual mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_B2', 'allele_C1'], self.allele_data)
        self.assertEqual(Allele_test.create_indicator_string(), ['---CC--CTT---', '---CC--CGC---', ])
        self.assertEqual(Allele_test.indicator_string, '--------MX---')

        # Test case 4, aligned sequences are identical (which should not be possible)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_B1'], self.allele_data)
        with self.assertRaises(ValueError):
            Allele_test.create_indicator_string()

    def test_check_indicative_SNPs(self):
        """

        """

        # Test case 1, indicator string has enough indicative mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XYYYY---'
        self.assertEqual(Allele_test.check_indicative_SNPs(), True)

        # Test case 2, indicator string has not enough indicative mismatches (Y's)
        Allele_test.indicator_string = '------XXX-----XY----'
        self.assertEqual(Allele_test.check_indicative_SNPs(), False)

        # Test case 3, indicator string has not enough indicative mismatches (X's)
        Allele_test.indicator_string = '------YYY-----XY----'
        self.assertEqual(Allele_test.check_indicative_SNPs(), False)

        # Test case 4, indicator string has not enough indicative mismatches (X's and Y's)
        Allele_test.indicator_string = '--XY---------'
        self.assertEqual(Allele_test.check_indicative_SNPs(), False)



    def test_check_mutual_SNPs(self):
        """

        """
        # Test case 1, no mutual mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XYYYY---'



        self.assertEqual(Allele_test.check_mutual_SNPs(), True)
        self.assertEqual(Allele_test.number_of_artefacts, 0)

        # Test case 2, one mutual mismatch (is accepted)
        Allele_test.indicator_string = '------XXM-----XYYYY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), True)
        self.assertEqual(Allele_test.number_of_artefacts, 1)

        # Test case 3, two mutual mismatches (is accepted)
        Allele_test.indicator_string = '------XXM-----XYMY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), True)
        self.assertEqual(Allele_test.number_of_artefacts, 3)

        # Test case 4, three mutual mismatches (is not accepted)
        Allele_test.indicator_string = '------XXM--M--XYMY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), False)
        self.assertEqual(Allele_test.number_of_artefacts, 6)


    def test_check_alternately_SNPs(self):
         
        # Test case 1, no alternately mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XYYYY---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,4,1,2,3,4])
        self.assertEqual(number_of_artefacts, 0)

        # Test case 2, 1 alternately mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XYX-----XX---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,1,1,2,3])
        self.assertEqual(number_of_artefacts, 2)

       # Test case 2, 1 alternately mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XX---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,4,5])
        self.assertEqual(number_of_artefacts, 0)


if __name__ == '__main__':
    unittest.main()



