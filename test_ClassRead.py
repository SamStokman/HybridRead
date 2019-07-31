"""
22-07-'19

This script contains 6 unittests for the class Read from SelectHybridReads.py.
The test can be ran with the bash command line: python3 test_ClassRead.py

"""

import unittest
import SelectHybridReads

class TestRead(unittest.TestCase):
    """
    This class contains unittests for the methods check_alignment(), apply_qv(), check_read_artefacts(),
    get_mismatches(), classmethod_for_non_read() and get_relative_position().
    """

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
        """
        True if the alignment of the read sequence is identical to the read sequence in alignment without the hyphen characters.
        False if the alignment of the read sequence is not identical to the read sequence in alignment without the hyphen characters. 
        ValueError raised if the read sequence is an empty string.
        """

        allele_data = {}
        # Test case 1: correct alignemt
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
        self.assertEqual(Read_test.check_alignment(), True)

        # Test case 2: incorrect alignment
        read_seq = 'CCACC'
        Read_test = SelectHybridReads.Read(read_seq, self.read_aligned_seq, allele_data)
        self.assertEqual(Read_test.check_alignment(), False)

        #Test case 3: empty alignment
        read_seq = ''
        Read_test = SelectHybridReads.Read(read_seq, self.read_aligned_seq, allele_data)
        with self.assertRaises(ValueError):
            Read_test.check_alignment()

    
    def test_apply_qv(self):
        """
        The quality value (qv) can be set by the user, therefore only the smallest and largest qv's are taken
        into account. Character '!' indicates a qv of 0 (lowest) and character 'I' indicates a qv if 40 (highest).
        If the qv is low, then the nucleotide must be replaced by an 'N'. If the qv is high then the nucleotide
        should not change.
        """

        allele_data = {}

        #Test case 1: low quality values
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
        read1_qv = '!!!!!'
        self.assertEqual(Read_test.apply_qv(read1_qv), '---NN--NNN---')
        
        #Test case 2: high quality values
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
        read1_qv = 'IIIII'
        self.assertEqual(Read_test.apply_qv(read1_qv), '---CC--CCC---')

        #Test case 3: high and low quality values
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
        read1_qv = 'I!!I!'
        self.assertEqual(Read_test.apply_qv(read1_qv), '---CN--NCN---')
        
    def test_check_read_artefacts(self):
        """
        If the aligned read has no mismatches, then the read alignment after the first check is identical to the input alignment.
        If all alleles have a mismatch, then the nucleotide at that postion is replaced by an 'N'
        If some alleles have a mismatch, then the read alignment after the first check is identical to the input alignment.
        If the data does not contain enough alleles, then an error is raised.
        """

        #Test case 1: no mismatches
        Read_alignment_after_first_check = '---CC--CCN---'
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, self.allele_data1)
        self.assertEqual(Read_test.check_read_artefacts(Read_alignment_after_first_check), '---CC--CCN---')

        #Test case 2: all alleles have 1 mismatches
        Read_alignment_after_first_check = '---CC--CCN---'
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, self.allele_data2)
        self.assertEqual(Read_test.check_read_artefacts(Read_alignment_after_first_check), '---CN--CCN---')
        
        #Test case 3: some alleles have 1 mismatches, but not all
        allele_data =[['allele_A1','---CT--CCC---'],
                      ['allele_A2','---CT--CCC---'],
                      ['allele_B1','---CC--CCC---'],
                      ['allele_B2','---CC--CCC---'],
                      ['allele_C1','---CT--CCC---'],
                      ['allele_C2','---CT--CCC---']]
        Read_alignment_after_first_check = '---CC--CCN---'
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, allele_data)
        self.assertEqual(Read_test.check_read_artefacts(Read_alignment_after_first_check), '---CC--CCN---')
    
        #Test case 4: not enough alleles
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
        """
        Substitutions, insertions and deletions are mismatches. 
        R_mismatch_dict should contain the total number of mismatches per allele.
        R_mismatch_dict_ex should contain substitutions, insertions, deletions and the total number of mismatches per alleles.
        """

        #Test case 1: no mismatches
        Read_alignment_after_second_check = '---CC--CCN---'
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, self.allele_data1)
        R_mismatch_dict, R_mismatch_dict_ex = Read_test.get_mismatches(Read_alignment_after_second_check)
        self.assertDictEqual(R_mismatch_dict, {'allele_A1': [0],
                                               'allele_A2': [0],
                                               'allele_B1': [0],
                                               'allele_B2': [0],
                                               'allele_C1': [0],
                                               'allele_C2': [0]})
        # dict values: substitutions, insertions, deletions and total
        self.assertDictEqual(R_mismatch_dict_ex, {'allele_A1': [0,0,0,0],
                                                  'allele_A2': [0,0,0,0],
                                                  'allele_B1': [0,0,0,0],
                                                  'allele_B2': [0,0,0,0],
                                                  'allele_C1': [0,0,0,0],
                                                  'allele_C2': [0,0,0,0]})

        #Test case 2: all alleles have 1 mismatches, only substitutions
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

        #Test case 3: more complex situations including substitutions, insertions and deletions
        allele_data  = [['allele_A1','-CC----------'],
                        ['allele_A2','---TT--TTT---'],
                        ['allele_B1','CCCCCCCCCCCCC'],
                        ['allele_B2','GGGGGGGGGGGGG'],
                        ['allele_C1','-------CCG---'],  # check for private method __check_read_start()
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
        """
        The input for the read consensus and turnover region are in alignment. The class' constructor requires the sequence
        without the alignemnt.

        The allele data should not change.
        The allele in alignment should not change.
        The read sequence is the allele in alignment minus the hyphen characters.
        """

        #Test case 1: normal situation
        single_sequence_in_alignment = '---TTTTT---'
        single_allele_in_alignment = '-----CCC-----'
       
        test_classmethod_result = SelectHybridReads.Read.classmethod_for_non_read(single_sequence_in_alignment, single_allele_in_alignment)
        cls.assertEqual(test_classmethod_result.read_seq, 'TTTTT')
        cls.assertEqual(test_classmethod_result.read_aligned_seq, '---TTTTT---')
        cls.assertEqual(test_classmethod_result.allele_data, '-----CCC-----')

    def test_get_relative_position(self):
        """
        The read position relative to the allele must be collected (the positions of the allele nucleotides
        that the read covers). 
        If the read has a nucleotide and the allele a '-' then the position should be excluded.
        """
        
        #Test case 1: normal situation
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, self.allele_data1)
        self.assertDictEqual(Read_test.get_relative_position(), {'allele_A1': [0,1,2,3,4],
                                                                 'allele_A2': [0,1,2,3,4],
                                                                 'allele_B1': [0,1,2,3,4],
                                                                 'allele_B2': [0,1,2,3,4],
                                                                 'allele_C1': [0,1,2,3,4],
                                                                 'allele_C2': [0,1,2,3,4]})

        #Test case 2: complex situations, including allele deletions 
        self.read_seq = 'CCCCC'
        self.read_aligned_seq = '---CC--CCC---'
        self.allele_data1 = [['allele_A1','-------CCC---'],
                             ['allele_A2','-CCCC--CCC---'],
                             ['allele_B1','-----CCCCCCCC'],
                             ['allele_B2','---TT--TTT---'],
                             ['allele_C1','---CCCCCCC---'],
                             ['allele_C2','C-C-C-C-C-C-C']]
        Read_test = SelectHybridReads.Read(self.read_seq, self.read_aligned_seq, self.allele_data1)
        self.assertDictEqual(Read_test.get_relative_position(), {'allele_A1': [0,1,2],
                                                                 'allele_A2': [2,3,4,5,6],
                                                                 'allele_B1': [2,3,4],
                                                                 'allele_B2': [0,1,2,3,4],
                                                                 'allele_C1': [0,1,2,3,4,5,6],
                                                                 'allele_C2': [2,3,4]})

        #Test case 3: read starts in front of allele
        self.read_seq = 'CCCCCC'
        read_aligned_seq = 'CCCCCC-------'
        allele_data = [['allele_A1','--CCC--------'],
                       ['allele_A2','--CC-C-------'],
                       ['allele_B1','-CCCCCC------'],
                       ['allele_B2','C-C-C-C-C-C-C'],
                       ['allele_C1','----------C--'],
                       ['allele_C2','--CCC--------']]
        Read_test = SelectHybridReads.Read(self.read_seq, read_aligned_seq, allele_data)
        self.assertDictEqual(Read_test.get_relative_position(), {'allele_A1': [0,1,2],
                                                                 'allele_A2': [0,1,2],
                                                                 'allele_B1': [0,1,2,3,4],
                                                                 'allele_B2': [0,1,2],
                                                                 'allele_C1': [],
                                                                 'allele_C2': [0,1,2]})

if __name__ == '__main__':
    unittest.main()
