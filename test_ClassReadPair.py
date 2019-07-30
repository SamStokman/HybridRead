"""
30-07-'19

This script contains unittests for the class ReadPair from SelectHybridReads.py.


"""

import unittest
import SelectHybridReads

class TestReadPair(unittest.TestCase):
    """
    This class contains unittests

    """

    def setUp(self):
        """

        """
        self.read_seq1 = 'CCCCC'
        self.read_aligned_seq1 = '---CN--NCC---'
        self.read_seq2 = 'CCCCC'
        self.read_aligned_seq2 = '---CC--CCC---'


    def test_check_read_pair(self):
        """
        The read length of read sequence 1 and 2 have to have a minimum length and a maximum number of N's

        """
        # Test case 1, read lengths and number of N's met requirement 
        Read_test = SelectHybridReads.ReadPair(self.read_aligned_seq1, self.read_aligned_seq2, self.read_seq1, self.read_seq2)
        min_read_length = 4
        N_quantity = 5
        self.assertEqual(Read_test.check_read_pair(min_read_length, N_quantity), True)

        # Test case 2, read lengths do not met requirement, number of N's does met requirement
        min_read_length = 6
        N_quantity = 1
        self.assertEqual(Read_test.check_read_pair(min_read_length, N_quantity), False)

        # Test case 3, read lengths do met requirement, number of N's does met requirement
        min_read_length = 4
        N_quantity = 1
        self.assertEqual(Read_test.check_read_pair(min_read_length, N_quantity), False)
        
        # Test case 4, read lengths and number of N's does met requirement
        min_read_length = 6
        N_quantity = 1
        self.read_aligned_seq1 = '---CN--NCC---'
        self.assertEqual(Read_test.check_read_pair(min_read_length, N_quantity), False)
    
    def test_create_read_consensus(self):

        # Test case 1, reads overlap completely, 1 read has N's, the other not
        Read_test = SelectHybridReads.ReadPair(self.read_aligned_seq1, self.read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(Read_test.create_read_consensus(), '---CC--CCC---')

        # Test case 2, reads have a gap
        read_aligned_seq1 = '---CCCC------------'
        read_aligned_seq2 = '----------TTTTTT---'
        Read_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(Read_test.create_read_consensus(), '---CCCC***TTTTTT---')

        # Test case 3, reads have a different length
        read_aligned_seq1 = '---CCCC------'
        read_aligned_seq2 = '----------TTTTTT---'
        Read_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        with self.assertRaises(ValueError):
            Read_test.create_read_consensus()

        # Test case 4, distinguish the gap between reads from deletions in the reads themselves.
        read_aligned_seq1 = '---CC-C------------'
        read_aligned_seq2 = '----------TTT--T---'
        Read_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(Read_test.create_read_consensus(), '---CC-C***TTT--T---')

        # Test case 5, same nucleotide position, different nucleotide 
        read_aligned_seq1 = '---CCCCCCC--------'
        read_aligned_seq2 = '-------TTTTTTT----'
        Read_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(Read_test.create_read_consensus(), '---CCCCNNNTTTT----')

        # Test case 6, special case, not a gap since they overlap
        read_aligned_seq1 = '---CCCCCCC-----C--'
        read_aligned_seq2 = '-------TTTTTTT----'
        Read_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(Read_test.create_read_consensus(), '---CCCCNNNNNNN-C--')

if __name__ == '__main__':
    unittest.main()

