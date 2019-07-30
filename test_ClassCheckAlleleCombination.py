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
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXM-----XYYYY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), True)
        self.assertEqual(Allele_test.number_of_artefacts, 1)

        # Test case 3, two mutual mismatches (is accepted)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXM-----XYMY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), True)
        self.assertEqual(Allele_test.number_of_artefacts, 2)

        # Test case 4, three mutual mismatches (is not accepted)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXM--M--XYMY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), False)
        self.assertEqual(Allele_test.number_of_artefacts, 3)


    def test_check_alternately_SNPs(self):
         
        # Test case 1, no alternately mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XYYYY---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,4,1,2,3,4])
        self.assertEqual(number_of_artefacts, 0)

        # Test case 2, 1 alternately mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '-----X-XYX-----XX---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,1,1,2,3])
        self.assertEqual(number_of_artefacts, 1)

        # Test case 3, 0 alternately mismatches, only 'X's'
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XX---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,4,5])
        self.assertEqual(number_of_artefacts, 0)

        # Test case 4, too many alternately mismatches, 3
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '-----X-XYX--XX---YXXYYYYXYYYY---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, None)
        self.assertEqual(number_of_artefacts, 3)

        # Test case 5, two 'Y's'instead of 1, does not count as alternatively mismatch (correct = a more switches hybrid)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '-----X-XYYX-----XX---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,1,2,1,2,3])
        self.assertEqual(number_of_artefacts, 0)

        # Test case 6, starts with a single 'X' three 1's at the start
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XYX-----XXYYYY-'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,1,1,2,3,1,2,3,4])
        self.assertEqual(number_of_artefacts, 1)

         # Test case 7, two alternatively mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYXXXX-YYYXYYY-'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,1,1,2,3,4,1,2,3,1,1,2,3])
        self.assertEqual(number_of_artefacts, 2)

         # Test case 8, 3 alternatively mismatches next to each other
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '----XXXYX-Y-XX-'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, None)
        self.assertEqual(number_of_artefacts, 3)

        # Test case 9, 2 alternatively mismatches, three 1's in the middle of the count indicator list
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '----XXXYX-YY-XX-'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,1,1,1,2,1,2])
        self.assertEqual(number_of_artefacts, 2)

    def test_update_indicator_string(self):
        """

        """
        # Test case 1, no changes needed (no mismatches)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XYYYY---'
        count_indicator_list = [1,2,3,4,1,2,3,4]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '------XXX-----XYYYY---')

        # Test case 2, 1 alternatively mismatches 
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '-----X-XYX-----XX---'
        count_indicator_list = [1,2,1,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '-----X-X-X-----XX---')

        # Test case 3, 1 mutual mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-M---XYYYY---'
        count_indicator_list = [1,2,3,4,1,2,3,4]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '------XXX-----XYYYY---')

        # Test case 4, 2 alternatively mismatches 
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYXXXX-YYYXYYY-'
        count_indicator_list = [1,2,3,1,1,2,3,4,1,2,3,1,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '--XXX-XXXX-YYY-YYY-')

        # Test case 5, 1 alternatively mismatch and 1 mutual mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYMXXX-YYYYYY-'
        count_indicator_list = [1,2,3,1,1,2,3,4,1,2,3,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '--XXX--XXX-YYYYYY-')

        # Test case 6, 2 alternatively mismatches and 1 mutual mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYMXXX-YYYXYYY-'
        count_indicator_list = [1,2,3,1,1,2,3,1,2,3,1,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '--XXX--XXX-YYY-YYY-')

        # Test case 7, 2 alternatively mismatches and 2 mutual mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYMXXX-YYYXYMYY-'
        count_indicator_list = [1,2,3,1,1,2,3,1,2,3,1,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '--XXX--XXX-YYY-Y-YY-')

    def test_get_switches(self):
        """
        """

        # Test case 1, 1 switch, turn over region of length 1
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--XXX-YYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)

        self.assertEqual(nr_of_switches, 1)
        self.assertEqual(start_turn_pos, 10)
        self.assertEqual(end_turn_pos, 10)

        # Test case 2, 1 switch, turn over region of length 8
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--XXX------YYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)
        self.assertEqual(nr_of_switches, 1)
        self.assertEqual(start_turn_pos, 10)
        self.assertEqual(end_turn_pos, 15)
        
        # Test case 3, 1 switch, turn over region of length 0
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--XXXYYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)
        self.assertEqual(nr_of_switches, 1)
        self.assertEqual(start_turn_pos, 10)
        self.assertEqual(end_turn_pos, 9)

        # Test case 4, 3 switches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--YYXX------YYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)
        self.assertEqual(nr_of_switches, 3)
        self.assertEqual(start_turn_pos, None)
        self.assertEqual(end_turn_pos, None)

        # Test case 5, 5 switches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--YYXX---YX---YYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)
        self.assertEqual(nr_of_switches, 5)
        self.assertEqual(start_turn_pos, None)
        self.assertEqual(end_turn_pos, None)

if __name__ == '__main__':
    unittest.main()



