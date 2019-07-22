"""
22-07-'19

A few unit test for SelectHybridReads.py


"""

import unittest
import SelectHybridReads

class TestRead(unittest.TestCase):

    def test_check_alignment(self):

        allele_data = {}
        # Test case 1
        read_seq = 'CCCCC'
        read_aligned_seq = '---CC--CCC---'
        Read_test = SelectHybridReads.Read(read_seq, read_aligned_seq, allele_data)
       # result = Read_test.check_alignment()
        self.assertEqual(Read_test.check_alignment(), True)

        # Test case 2
        read_seq = 'CCACC'
        read_aligned_seq = '---CC--CCC---'
        Read_test = SelectHybridReads.Read(read_seq, read_aligned_seq, allele_data)
        self.assertEqual(Read_test.check_alignment(), False)

        #Test case 3
        read_seq = ''
        read_aligned_seq = '---CC--CCC---'
        Read_test = SelectHybridReads.Read(read_seq, read_aligned_seq, allele_data)
        with self.assertRaises(ValueError):
            Read_test.check_alignment()


    
    def test_apply_qv(self):

        allele_data = {}

        #Test case 1, low quality values
        read_seq = 'CCCCC'
        read_aligned_seq = '---CC--CCC---'
        Read_test = SelectHybridReads.Read(read_seq, read_aligned_seq, allele_data)
        read1_qv = '!!!!!'
        self.assertEqual(Read_test.apply_qv(read1_qv), '---NN--NNN---')
        
        #Test case 1, high quality values
        read_seq = 'CCCCC'
        read_aligned_seq = '---CC--CCC---'
        Read_test = SelectHybridReads.Read(read_seq, read_aligned_seq, allele_data)
        read1_qv = 'IIIII'
        self.assertEqual(Read_test.apply_qv(read1_qv), '---CC--CCC---')

        #Test case 1, high and low quality values
        read_seq = 'CCCCC'
        read_aligned_seq = '---CC--CCC---'
        Read_test = SelectHybridReads.Read(read_seq, read_aligned_seq, allele_data)
        read1_qv = 'I!!I!'
        self.assertEqual(Read_test.apply_qv(read1_qv), '---CN--NCN---')

if __name__ == '__main__':
    unittest.main()
