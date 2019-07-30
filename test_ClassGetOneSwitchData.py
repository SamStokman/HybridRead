"""
30-07-'19

This script contains unittests for the class GetOneSwitchData from SelectHybridReads.py.


"""

import unittest
import SelectHybridReads

class TestGetOneSwitchData(unittest.TestCase):
    """
    This class contains unittests

    """

    def setUp(self):
        """

        """
        self.allele1 = 'allele_A1'
        self.allele2 = 'allele_A2'
        self.allele_data = [['allele_A1','---CCCCCCCC--'],
                            ['allele_A2','---AAAAAAAA--'],
                            ['allele_B1','---CC--CCC---'],
                            ['allele_B2','---CC--CCC---'],
                            ['allele_C1','---CC--CCC---'],
                            ['allele_C2','---CC--CCC---']]

    def test_get_read_position(self):
        
        # test case 1, standard values
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

    def test_prep_for_turnover_position(self):
        """

        """

        # Test case 1 
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 4
        end_pos = 7
        allele_seq_list = ['---CCCCCCCC--', '---AAAAAAAA--']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)

        self.assertEqual(turn_over_region1_for_pos, '----CCCC-----')
        self.assertEqual(seq_list_allele1, [('allele_A1', '---CCCCCCCC--')])
        self.assertEqual(turn_over_region2_for_pos, '----AAAA-----')
        self.assertEqual(seq_list_allele2, [('allele_A2', '---AAAAAAAA--')])

        # The allele data was not adjusted for test cases 2 - 5,  and therefore that test is the same each case.

        # Test case 2 
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 4
        end_pos = 20
        allele_seq_list = ['ACCGAC----------AGGACCGCG', 'ACCGAC----------AGGACCGCG']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)

        self.assertEqual(turn_over_region1_for_pos, '----AC----------AGGAC----')
        self.assertEqual(turn_over_region2_for_pos, '----AC----------AGGAC----')

        # Test case 3, turn over region length 1
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 6
        end_pos = 6
        allele_seq_list = ['GGGGGGGGGG', 'TTTTTTTTTT']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)

        self.assertEqual(turn_over_region1_for_pos, '------G---')
        self.assertEqual(turn_over_region2_for_pos, '------T---')

        # Test case 4, turn over region length 1, no nucleotide at the turn over region position for first allele
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 6
        end_pos = 6
        allele_seq_list = ['GGGGGG-GGG', 'TTTTTTTTTT']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)

        self.assertEqual(turn_over_region1_for_pos, '------Z---')
        self.assertEqual(turn_over_region2_for_pos, '------T---')
    
        # Test case 5, turn over region length 0
        Switch_test = SelectHybridReads.GetOneSwitchData(self.allele1, self.allele2)
        start_pos = 6
        end_pos = 5
        allele_seq_list = ['GGGGGGGGGG', 'TTTTTT-TTT']
        turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = Switch_test.prep_for_turnover_position(start_pos, end_pos, allele_seq_list, self.allele_data)

        self.assertEqual(turn_over_region1_for_pos, '------K---')
        self.assertEqual(turn_over_region2_for_pos, '------K---')

    def test_get_TO_position(self):
        """

        """
        
        # Test case 1
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

if __name__ == '__main__':
    unittest.main()


