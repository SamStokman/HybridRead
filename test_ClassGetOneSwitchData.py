"""
30-07-'19

This script contains 3 unittests for the class GetOneSwitchData from SelectHybridReads.py.
The test can be ran with the bash command line: python3 test_ClassGetOneSwitchData.py
"""

import unittest
import SelectHybridReads

class TestGetOneSwitchData(unittest.TestCase):
    """
    This class contains unittests for the methods get_read_position(), prep_for_turnover_position() and get_TO_position(). 
    """

    def setUp(self):
        self.allele1 = 'allele_A1'
        self.allele2 = 'allele_A2'
        self.allele_data = [['allele_A1','---CCCCCCCC--'],
                            ['allele_A2','---AAAAAAAA--'],
                            ['allele_B1','---CC--CCC---'],
                            ['allele_B2','---CC--CCC---'],
                            ['allele_C1','---CC--CCC---'],
                            ['allele_C2','---CC--CCC---']]

    def test_get_read_position(self):
        """
        This method only has a parsing purpose.
        It gets the first and last postion from the dict values for the given alleles and puts them in a string.
        If the there is just 1 position (rare event), then this position is used twice in the string.
        """
        
        #Test case 1: standard situation
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        read1_pos_dict = {'allele_A1': [2,3,4,5],
                          'allele_A2': [2,3,4,5],
                          'allele_B1': [2,3,4,5],
                          'allele_B2': [2,3,4,5],
                          'allele_C1': [2,3,4,5],
                          'allele_C2': [2,3,4,5]}
        read2_pos_dict = {'allele_A1': [2,3,4,5],
                          'allele_A2': [2,3,4,5],
                          'allele_B1': [2,3,4,5],
                          'allele_B2': [2,3,4,5],
                          'allele_C1': [2,3,4,5],
                          'allele_C2': [2,3,4,5]}
        position_read1_allele1, position_read2_allele1, position_read1_allele2, position_read2_allele2 = Switch_test.get_read_position(read1_pos_dict, read2_pos_dict)
        self.assertEqual(position_read1_allele1, '2-5')
        self.assertEqual(position_read2_allele1, '2-5')
        self.assertEqual(position_read1_allele2, '2-5')
        self.assertEqual(position_read2_allele2, '2-5')
                
        #Test case 2: 1 relative read position
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        read1_pos_dict = {'allele_A1': [2],
                          'allele_A2': [2]}
        read2_pos_dict = {'allele_A1': [8],
                          'allele_A2': [5]}
        position_read1_allele1, position_read2_allele1, position_read1_allele2, position_read2_allele2 = Switch_test.get_read_position(read1_pos_dict, read2_pos_dict)
        self.assertEqual(position_read1_allele1, '2-2')
        self.assertEqual(position_read2_allele1, '8-8')
        self.assertEqual(position_read1_allele2, '2-2')
        self.assertEqual(position_read2_allele2, '5-5')

    def test_prep_for_turnover_position(self):
        """
        Here, the tunover region in alignment for both alleles is created.
        The allele sequence is copied (nucleotides and '-'s) for the turnover region positions, the '-'s from the alignment are also included.
        If the turnover region length is 1 and the sequence for the allele is '-', then a 'Z' is used.
        If the turnover region length is 1, it does not matter if the allele has a nucleotide or '-', then a 'K' is used at the start_turn_pos 
        (the position of the second indicator character).
        """

        #Test case 1: standard situation 
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 4
        end_pos = 7
        allele_seq_list = ['---CCCCCCCC--', '---AAAAAAAA--']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)
        self.assertEqual(turn_over_region1_for_pos, '----CCCC-----')
        self.assertEqual(seq_list_allele1, [('allele_A1', '---CCCCCCCC--')])
        self.assertEqual(turn_over_region2_for_pos, '----AAAA-----')
        self.assertEqual(seq_list_allele2, [('allele_A2', '---AAAAAAAA--')])

        #The allele data was not adjusted for test cases 2 - 5,  and therefore that test is the same each case.
        #Test case 2: turnover region with large gap
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 4
        end_pos = 20
        allele_seq_list = ['ACCGAC----------AGGACCGCG', 'ACCGAC----------AGGACCGCG']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)
        self.assertEqual(turn_over_region1_for_pos, '----AC----------AGGAC----')
        self.assertEqual(turn_over_region2_for_pos, '----AC----------AGGAC----')

        #Test case 3: turnover region length 1
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 6
        end_pos = 6
        allele_seq_list = ['GGGGGGGGGG', 'TTTTTTTTTT']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)
        self.assertEqual(turn_over_region1_for_pos, '------G---')
        self.assertEqual(turn_over_region2_for_pos, '------T---')

        #Test case 4: turnover region length 1, no nucleotide at the turn over region position for first allele
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 6
        end_pos = 6
        allele_seq_list = ['GGGGGG-GGG', 'TTTTTTTTTT']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)
        self.assertEqual(turn_over_region1_for_pos, '------Z---')
        self.assertEqual(turn_over_region2_for_pos, '------T---')
    
        #Test case 5: turnover region length 0
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 6
        end_pos = 5
        allele_seq_list = ['GGGGGGGGGG', 'TTTTTT-TTT']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)
        self.assertEqual(turn_over_region1_for_pos, '------K---')
        self.assertEqual(turn_over_region2_for_pos, '------K---')

    def test_get_TO_position(self):
        """
        This method only has parsing purposes.
        It gets the first and last postion from the dict values for the given alleles and puts them in a string.
        The '-'s in the aligned turnover region are deleted, which result in the turnover region sequence.
        If the aligned turnover region contains a 'K', the sequence must be empty (the turnover region length is 0) and the 
        position_to_region is the first value in de TO_allele_dict.
        """
        
        #Test case 1: standard situation
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        TO_allele1_dict = {'allele_A1': [2,3,4,5]}
        TO_allele2_dict = {'allele_A2': [2,3,4,5]}
        turn_over_region1_for_pos = '--CCCC-----'
        turn_over_region2_for_pos = '--TTTT-----'
        position_to_region1, position_to_region2, turn_over_region1, turn_over_region2 = Switch_test.get_TO_position(TO_allele1_dict, TO_allele2_dict, turn_over_region1_for_pos, turn_over_region2_for_pos)
        self.assertEqual(position_to_region1, '2-5')
        self.assertEqual(position_to_region2, '2-5')
        self.assertEqual(turn_over_region1, 'CCCC')
        self.assertEqual(turn_over_region2, 'TTTT')

        #Test case 2: turnover region with length 1
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        TO_allele1_dict = {'allele_A1': [2,2]}
        TO_allele2_dict = {'allele_A2': [2,2]}
        turn_over_region1_for_pos = '--C------'
        turn_over_region2_for_pos = '--T------'
        position_to_region1, position_to_region2, turn_over_region1, turn_over_region2 = Switch_test.get_TO_position(TO_allele1_dict, TO_allele2_dict, turn_over_region1_for_pos, turn_over_region2_for_pos)
        self.assertEqual(position_to_region1, '2')
        self.assertEqual(position_to_region2, '2')
        self.assertEqual(turn_over_region1, 'C')
        self.assertEqual(turn_over_region2, 'T')

        #Test case 3: turnover region with length 0
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        TO_allele1_dict = {'allele_A1': [2,1]}
        TO_allele2_dict = {'allele_A2': [2,1]}
        turn_over_region1_for_pos = '--K------'
        turn_over_region2_for_pos = '--K------'
        position_to_region1, position_to_region2, turn_over_region1, turn_over_region2 = Switch_test.get_TO_position(TO_allele1_dict, TO_allele2_dict, turn_over_region1_for_pos, turn_over_region2_for_pos)
        self.assertEqual(position_to_region1, '2')
        self.assertEqual(position_to_region2, '2')
        self.assertEqual(turn_over_region1, '')
        self.assertEqual(turn_over_region2, '')

if __name__ == '__main__':
    unittest.main()


