"""
22-07-'19

A few unit test for SelectHybridReads.py


"""

import unittest
import SelectHybridReads

class TestRead(unittest.TestCase):

    def setUp(self):
        self.read_seq = 'CCCCC'
        self.read_aligned_seq = '---CC--CCC---'
        self.allele_data1 = [['allele_A1','---CC--CCC---'],
                             ['allele_A2','---CC--CCC---'],
                             ['allele_B1','---CC--CCC---'],
                             ['allele_B2','---CC--CCC---'],
                             ['allele_C1','---CC--CCC---'],
                             ['allele_C2','---CC--CCC---']]
               
        self.allele_data2 =[['allele_A1','---CT--CCC---'],
                            ['allele_A2','---CT--CCC---'],
                            ['allele_B1','---CT--CCC---'],
                            ['allele_B2','---CT--CCC---'],
                            ['allele_C1','---CT--CCC---'],
                            ['allele_C2','---CT--CCC---']]

    def test_check_alignment(self):

        allele_data = {}
        # Test case 1
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
       # result = Read_test.check_alignment()
        self.assertEqual(Read_test.check_alignment(), True)

        # Test case 2
        read_seq = 'CCACC'
        Read_test = SelectHybridReads.Read(read_seq, self.read_aligned_seq, allele_data)
        self.assertEqual(Read_test.check_alignment(), False)

        #Test case 3
        read_seq = ''
        Read_test = SelectHybridReads.Read(read_seq, self.read_aligned_seq, allele_data)
        with self.assertRaises(ValueError):
            Read_test.check_alignment()

    
    def test_apply_qv(self):

        allele_data = {}

        #Test case 1, low quality values
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
        read1_qv = '!!!!!'
        self.assertEqual(Read_test.apply_qv(read1_qv), '---NN--NNN---')
        
        #Test case 1, high quality values
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
        read1_qv = 'IIIII'
        self.assertEqual(Read_test.apply_qv(read1_qv), '---CC--CCC---')

        #Test case 1, high and low quality values
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
        read1_qv = 'I!!I!'
        self.assertEqual(Read_test.apply_qv(read1_qv), '---CN--NCN---')

        
    def test_check_read_artefacts(self):
        
        #Test case 1, no mismatches


        Read_alignment_after_first_check = '---CC--CCN---'
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, self.allele_data1)
        self.assertEqual(Read_test.check_read_artefacts(Read_alignment_after_first_check), '---CC--CCN---')

        #Test case 2, all alleles have 1 mismatches

        Read_alignment_after_first_check = '---CC--CCN---'

        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, self.allele_data2)
        self.assertEqual(Read_test.check_read_artefacts(Read_alignment_after_first_check), '---CN--CCN---')
      
        #Test case 3, not enough alleles
        allele_data =[['allele_A1','---CT--CCC---'],
                      ['allele_A2','---CT--CCC---'],
                      ['allele_B1','---CT--CCC---'],
                      ['allele_C2','---CT--CCC---']]

        Read_alignment_after_first_check = '---CC--CCN---'
        read_seq = 'CCCCC'
        read_aligned_seq = '---CC--CCC---'
        Read_test = SelectHybridReads.Read(read_seq, read_aligned_seq, allele_data)
        with self.assertRaises(ValueError):
            Read_test.check_read_artefacts(Read_alignment_after_first_check)
    
            
    def test_get_mismatches(self):

        # Test case 1 no mismatches

        Read_alignment_after_second_check = '---CC--CCN---'
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, self.allele_data1)
        R_mismatch_dict, R_mismatch_dict_ex = Read_test.get_mismatches(Read_alignment_after_second_check)
        self.assertDictEqual(R_mismatch_dict, {'allele_A1': [0],
                                               'allele_A2': [0],
                                               'allele_B1': [0],
                                               'allele_B2': [0],
                                               'allele_C1': [0],
                                               'allele_C2': [0]})
        # values: substitutions, insertions, deletions and total
        self.assertDictEqual(R_mismatch_dict_ex, {'allele_A1': [0,0,0,0],
                                                  'allele_A2': [0,0,0,0],
                                                  'allele_B1': [0,0,0,0],
                                                  'allele_B2': [0,0,0,0],
                                                  'allele_C1': [0,0,0,0],
                                                  'allele_C2': [0,0,0,0]})

        
        #Test case 2, all alleles have 1 mismatches
        Read_alignment_after_second_check = '---CC--CCN---'
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, self.allele_data2)
        R_mismatch_dict, R_mismatch_dict_ex = Read_test.get_mismatches(Read_alignment_after_second_check)
        self.assertDictEqual(R_mismatch_dict, {'allele_A1': [1],
                                               'allele_A2': [1],
                                               'allele_B1': [1],
                                               'allele_B2': [1],
                                               'allele_C1': [1],
                                               'allele_C2': [1]})

        self.assertDictEqual(R_mismatch_dict_ex, {'allele_A1': [1,0,0,1],
                                                  'allele_A2': [1,0,0,1],
                                                  'allele_B1': [1,0,0,1],
                                                  'allele_B2': [1,0,0,1],
                                                  'allele_C1': [1,0,0,1],
                                                  'allele_C2': [1,0,0,1]})
       # Test case random mismatches
        allele_data  = [['allele_A1','-CC----------'],
                        ['allele_A2','---TT--TTT---'],
                        ['allele_B1','CCCCCCCCCCCCC'],
                        ['allele_B2','GGGGGGGGGGGGG'],
                        ['allele_C1','-------CCG---'],  # check for __check_read_start
                        ['allele_C2','CCC--CC---CCC']]
       
        Read_alignment_after_second_check = '---CC--CCN---'
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
        R_mismatch_dict, R_mismatch_dict_ex = Read_test.get_mismatches(Read_alignment_after_second_check)
        self.assertDictEqual(R_mismatch_dict, {'allele_A1': [5],
                                               'allele_A2': [4],
                                               'allele_B1': [2],
                                               'allele_B2': [6],
                                               'allele_C1': [0],
                                               'allele_C2': [7]})

        self.assertDictEqual(R_mismatch_dict_ex, {'allele_A1': [0,0,5,5],
                                                  'allele_A2': [4,0,0,4],
                                                  'allele_B1': [0,2,0,2],
                                                  'allele_B2': [4,2,0,6],
                                                  'allele_C1': [0,0,0,0],
                                                  'allele_C2': [0,2,5,7]})

    def test_classmethod_for_non_read(cls):
        single_sequence_in_alignment = '---TTTTT---'
        single_allele_in_alignment = '-----CCC-----'
       
        test_classmethod_result = SelectHybridReads.Read.classmethod_for_non_read(single_sequence_in_alignment, single_allele_in_alignment)
        cls.assertEqual(test_classmethod_result.read_seq, 'TTTTT')
        cls.assertEqual(test_classmethod_result.read_aligned_seq, '---TTTTT---')
        cls.assertEqual(test_classmethod_result.allele_data, '-----CCC-----')


if __name__ == '__main__':
    unittest.main()
