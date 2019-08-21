"""
30-07-'19

This script contains 2 unittests for the class ReadPair from SelectHybridReads.py.
The test can be ran with the bash command line: python3 test_ClassReadPair.py
"""

import unittest
import SelectHybridReads

class TestReadPair(unittest.TestCase):
    """
    This class contains unittests for the methods check_read_pair() and create_read_consensus(). 
    """

    def setUp(self):
        self.read_seq1 = 'CCCCC'
        self.read_aligned_seq1 = '---CN--NCC---'
        self.read_seq2 = 'CCCCC'
        self.read_aligned_seq2 = '---CC--CCC---'

    def test_check_read_pair(self):
        """
        The paired-end reads are checked for a minimum length and a maximum number of N's.
        The length default is 80 and number of N's is 5, both requirements can be set by the user.
        These requirements should be met for each individual read per read pair.
        If both individual reads of the paired-end read met the requirements, then the boolean must return True.
        In all other cases, the boolean must return False.
        Note that the input sequences are identical for each test case, the requirement values are varied.
        """

        #Test case 1: read lengths and number of N's met requirement 
        ReadPair_test = SelectHybridReads.ReadPair(self.read_aligned_seq1, self.read_aligned_seq2, self.read_seq1, self.read_seq2)
        min_read_length = 4
        N_quantity = 5
        self.assertEqual(ReadPair_test.check_read_pair(min_read_length, N_quantity), True)

        #Test case 2: read lengths do not met requirement, number of N's does met requirement
        min_read_length = 6
        N_quantity = 1
        self.assertEqual(ReadPair_test.check_read_pair(min_read_length, N_quantity), False)

        #Test case 3: read lengths do met requirement, number of N's does met requirement
        min_read_length = 4
        N_quantity = 1
        self.assertEqual(ReadPair_test.check_read_pair(min_read_length, N_quantity), False)
        
        #Test case 4: read lengths and number of N's does met requirement
        min_read_length = 6
        N_quantity = 1
        self.read_aligned_seq1 = '---CN--NCC---'
        self.assertEqual(ReadPair_test.check_read_pair(min_read_length, N_quantity), False)
    
    def test_create_read_consensus(self):
        """
        The read consensus must be created out of the paired-end read.
        If the reads overlap and 1 nucleotide is 'N' and the other an actual nucleotide, then the actual nucleotide should be used.
        If the reads overlap and both nucleotides differ ('N' excluded) then an 'N'should be used.
        If the reads do not overlap, the gap should be filled in with '*'.
        The gaps (deletions) in the reads should be '-'.
        """

        #Test case 1: reads overlap completely, 1 read has N's, the other not
        ReadPair_test = SelectHybridReads.ReadPair(self.read_aligned_seq1, self.read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(ReadPair_test.create_read_consensus(), '---CC--CCC---')

        #Test case 2: reads have a gap
        read_aligned_seq1 = '---CCCC------------'
        read_aligned_seq2 = '----------TTTTTT---'
        ReadPair_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(ReadPair_test.create_read_consensus(), '---CCCC***TTTTTT---')

        #Test case 3: reads have a different length
        read_aligned_seq1 = '---CCCC------'
        read_aligned_seq2 = '----------TTTTTT---'
        ReadPair_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        with self.assertRaises(ValueError):
            ReadPair_test.create_read_consensus()

        #Test case 4: distinguish the gap between reads from deletions in the reads themselves
        read_aligned_seq1 = '---CC-C------------'
        read_aligned_seq2 = '----------TTT--T---'
        ReadPair_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(ReadPair_test.create_read_consensus(), '---CC-C***TTT--T---')

        #Test case 5: same nucleotide position, different nucleotide 
        read_aligned_seq1 = '---CCCCCCC--------'
        read_aligned_seq2 = '-------TTTTTTT----'
        ReadPair_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(ReadPair_test.create_read_consensus(), '---CCCCNNNTTTT----')

        #Test case 6: complex situation, not a gap since they overlap
        read_aligned_seq1 = '---CCCCCCC-----C--'
        read_aligned_seq2 = '-------TTTTTTT----'
        ReadPair_test = SelectHybridReads.ReadPair(read_aligned_seq1, read_aligned_seq2, self.read_seq1, self.read_seq2)
        self.assertEqual(ReadPair_test.create_read_consensus(), '---CCCCNNNNNNN-C--')

if __name__ == '__main__':
    unittest.main()

