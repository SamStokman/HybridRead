"""
30-07-'19

This script contains 6 unittests for the class CheckAlleleCombination from SelectHybridReads.py.
The test can be ran with the bash command line: python3 test_ClassCheckAlleleCombination.py
"""

import unittest
import SelectHybridReads


class TestCheckAlleleCombination(unittest.TestCase):
    """
    This class contains unittests for the methods create_indicator_string(), check_indicative_SNPs(), check_mutual_SNPs(),
    check_alternately_SNPs(), update_indicator_string() and get_switches(). 
    """

    def setUp(self):
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
        The function return the aligened sequences of the given alleles without adjustments and the indicator string is created.
        If the first allele has a mismatch, an 'X' is added.
        If the second allele has a mismatch, an 'Y' is added.
        If both alleles have a mismatch, an 'M' is added (mutual mismatch).
        If the allele has a '-' or a match than '-' is added.
        """

        #Test case 1: second allele has a mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        self.assertEqual(Allele_test.create_indicator_string(), ['---CC--CCC---', '---CC--TCC---'])
        self.assertEqual(Allele_test.indicator_string, '-------Y-----')

        #Test case 2: first allele has a mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A2', 'allele_B1'], self.allele_data)
        self.assertEqual(Allele_test.create_indicator_string(), ['---CC--TCC---', '---CC--CCC---'])
        self.assertEqual(Allele_test.indicator_string, '-------X-----')

        #Test case 2: both alleles have a mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A2', 'allele_B2'], self.allele_data)
        self.assertEqual(Allele_test.create_indicator_string(), ['---CC--TCC---', '---CC--CTT---'])
        self.assertEqual(Allele_test.indicator_string, '-------XYY---')

        #Test case 3: alleles have mutual mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_B2', 'allele_C1'], self.allele_data)
        self.assertEqual(Allele_test.create_indicator_string(), ['---CC--CTT---', '---CC--CGC---', ])
        self.assertEqual(Allele_test.indicator_string, '--------MX---')

        #Test case 4: aligned sequences are identical (which should not be possible)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_B1'], self.allele_data)
        with self.assertRaises(ValueError):
            Allele_test.create_indicator_string()

    def test_check_indicative_SNPs(self):
        """
        A simple method to check whether the combined alleles have enough indicative mismatches or not.
        If the indicator string has at least two X's and two Y's, boolean is True.
        If the indicator string has less than two X's or less than 2 Y's, boolean is False.
        """

        # Test case 1: indicator string has enough indicative mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XYYYY---'
        self.assertEqual(Allele_test.check_indicative_SNPs(), True)

        # Test case 2: indicator string has not enough indicative mismatches (Y's)
        Allele_test.indicator_string = '------XXX-----XY----'
        self.assertEqual(Allele_test.check_indicative_SNPs(), False)

        # Test case 3: indicator string has not enough indicative mismatches (X's)
        Allele_test.indicator_string = '------YYY-----XY----'
        self.assertEqual(Allele_test.check_indicative_SNPs(), False)

        # Test case 4; indicator string has not enough indicative mismatches (X's and Y's)
        Allele_test.indicator_string = '--XY---------'
        self.assertEqual(Allele_test.check_indicative_SNPs(), False)

    def test_check_mutual_SNPs(self):
        """
        A simple method to check whether the combined alleles have not too many mutual mismatches.
        If the indicator string has two M's or less, boolean is True.
        If the indiactor string has more than two M's, booelean is False. 
        """

        #Test case 1: no mutual mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XYYYY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), True)
        self.assertEqual(Allele_test.number_of_artefacts, 0)

        #Test case 2: one mutual mismatch (is accepted)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXM-----XYYYY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), True)
        self.assertEqual(Allele_test.number_of_artefacts, 1)

        #Test case 3: two mutual mismatches (is accepted)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXM-----XYMY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), True)
        self.assertEqual(Allele_test.number_of_artefacts, 2)

        #Test case 4: three mutual mismatches (is not accepted)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXM--M--XYMY---'
        self.assertEqual(Allele_test.check_mutual_SNPs(), False)
        self.assertEqual(Allele_test.number_of_artefacts, 3)


    def test_check_alternately_SNPs(self):
        """
        Checks for alternately mismatches (artefact type 2).
        If the indicator string contains two or less alternately mismatches, then count indicator is created.
        If the indicator string contains more than two alternately mismatches, then count indicator is None.
        The count indicator list contains ascending values (int), starting from 1, each new start indicates a 
        switch for 'X' to 'Y' or vice versa. Two 1's next to each other indicate an alternately mismatch.  
        The number of artefacts is the number of alternately mismatches.
        """
         
        #Test case 1: no alternately mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XYYYY---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,4,1,2,3,4])
        self.assertEqual(number_of_artefacts, 0)

        #Test case 2: 1 alternately mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '-----X-XYX-----XX---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,1,1,2,3])
        self.assertEqual(number_of_artefacts, 1)

        #Test case 3: 0 alternately mismatches, only 'X's'
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XX---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,4,5])
        self.assertEqual(number_of_artefacts, 0)

        #Test case 4: too many alternately mismatches, 3
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '-----X-XYX--XX---YXXYYYYXYYYY---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, None)
        self.assertEqual(number_of_artefacts, 3)

        #Test case 5: two 'Y's'instead of 1, does not count as alternatively mismatch (correct = a more switches hybrid)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '-----X-XYYX-----XX---'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,1,2,1,2,3])
        self.assertEqual(number_of_artefacts, 0)

        #Test case 6: starts with a single 'X' three 1's at the start
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XYX-----XXYYYY-'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,1,1,2,3,1,2,3,4])
        self.assertEqual(number_of_artefacts, 1)

        #Test case 7: two alternatively mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYXXXX-YYYXYYY-'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,1,1,2,3,4,1,2,3,1,1,2,3])
        self.assertEqual(number_of_artefacts, 2)

        #Test case 8: 3 alternatively mismatches next to each other
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '----XXXYX-Y-XX-'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, None)
        self.assertEqual(number_of_artefacts, 3)

        #Test case 9: 2 alternatively mismatches, three 1's in the middle of the count indicator list
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '----XXXYX-YY-XX-'
        count_indicator_list, number_of_artefacts = Allele_test.check_alternately_SNPs()
        self.assertEqual(count_indicator_list, [1,2,3,1,1,1,2,1,2])
        self.assertEqual(number_of_artefacts, 2)

    def test_update_indicator_string(self):
        """
        Here, the mutual mismatches and alternatively mismatches are dismissed if there are maximal two of each (thus four in total).
        If the indicator string contains an 'M' then it should be replaced a '-'.
        If the indicator string contains an alternatively mismatch (XYX or YXY) the single nucleotide which causes the mismatch is replaced 
        by a '-'
        """

        #Test case 1: no changes needed (no mismatches)
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-----XYYYY---'
        count_indicator_list = [1,2,3,4,1,2,3,4]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '------XXX-----XYYYY---')

        #Test case 2: 1 alternatively mismatches 
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '-----X-XYX-----XX---'
        count_indicator_list = [1,2,1,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '-----X-X-X-----XX---')

        #Test case 3: 1 mutual mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '------XXX-M---XYYYY---'
        count_indicator_list = [1,2,3,4,1,2,3,4]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '------XXX-----XYYYY---')

        #Test case 4: 2 alternatively mismatches 
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYXXXX-YYYXYYY-'
        count_indicator_list = [1,2,3,1,1,2,3,4,1,2,3,1,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '--XXX-XXXX-YYY-YYY-')

        #Test case 5: 1 alternatively mismatch and 1 mutual mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYMXXX-YYYYYY-'
        count_indicator_list = [1,2,3,1,1,2,3,4,1,2,3,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '--XXX--XXX-YYYYYY-')

        #Test case 6: 2 alternatively mismatches and 1 mutual mismatch
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYMXXX-YYYXYYY-'
        count_indicator_list = [1,2,3,1,1,2,3,1,2,3,1,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '--XXX--XXX-YYY-YYY-')

        #Test case 7: 2 alternatively mismatches and 2 mutual mismatches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        Allele_test.indicator_string = '--XXXYMXXX-YYYXYMYY-'
        count_indicator_list = [1,2,3,1,1,2,3,1,2,3,1,1,2,3]
        self.assertEqual(Allele_test.update_indicator_string(count_indicator_list), '--XXX--XXX-YYY-Y-YY-')

    def test_get_switches(self):
        """
        Switch definition: if the indicator string goes from X to Y or vice versa. With or without '-' in between the indicator characters (X and Y).
        If the indicator string contains more than 1 switch, the start_turn_pos and end_turn_pos are None (they do not exist). 
        If the indicator string contains 1 switch, and the number of positions ('-') between the X and Y (= turnover region length) is larger
        than 1. The start_turn_pos is the position of the first '-' after the first indicator character and end_turn_pos is the position of 
        the first '-' in front of the second indicator character.
        If the indicator string contains 1 switch, and turnover region length is 1, then the start_turn_pos is the position of the first '-',
        and end_turn_pos is the position of the first '-' in front of the second indicator character. Which is the same position now.
        If the indicator string contains 1 switch, and turnover region length is 0, the indicator characters are directly next to eachother. Then
        the start_turn_pos is the position of second indicator character and the end_turn_pos the position of first indicator character.
        """

        #Test case 1: 1 switch, turn over region of length 1
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--XXX-YYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)
        self.assertEqual(nr_of_switches, 1)
        self.assertEqual(start_turn_pos, 10)
        self.assertEqual(end_turn_pos, 10)

        #Test case 2: 1 switch, turn over region of length 8
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--XXX------YYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)
        self.assertEqual(nr_of_switches, 1)
        self.assertEqual(start_turn_pos, 10)
        self.assertEqual(end_turn_pos, 15)
        
        #Test case 3: 1 switch, turn over region of length 0
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--XXXYYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)
        self.assertEqual(nr_of_switches, 1)
        self.assertEqual(start_turn_pos, 10)
        self.assertEqual(end_turn_pos, 9)

        #Test case 4: 3 switches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--YYXX------YYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)
        self.assertEqual(nr_of_switches, 3)
        self.assertEqual(start_turn_pos, None)
        self.assertEqual(end_turn_pos, None)

        #Test case 5: 5 switches
        Allele_test = SelectHybridReads.CheckAlleleCombination(self.read_consensus, ['allele_A1', 'allele_A2'], self.allele_data)
        final_indicator_string = '--XXX--YYXX---YX---YYY-Y-YY-'
        nr_of_switches, start_turn_pos, end_turn_pos = Allele_test.get_switches(final_indicator_string)
        self.assertEqual(nr_of_switches, 5)
        self.assertEqual(start_turn_pos, None)
        self.assertEqual(end_turn_pos, None)

if __name__ == '__main__':
    unittest.main()



