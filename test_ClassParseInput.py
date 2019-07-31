"""
30-07-'19

This script contains 2 unittests for the class ParseInput from SelectHybridReads.py.
The test can be ran with the bash command line: python3 test_ClassParseInput.py
"""

import unittest
import SelectHybridReads

class TestParseInput(unittest.TestCase):
    """
    This class contains unittests for the method collect_all_data() and get_allele_combinations().

    """

    def test_collect_all_data(self):
        """
        The input data is a string that contains all lines from the given input file.
        The data needs to be separated and the data of each paired-end read needs to be stored in a list of lists.
        The all_data list must collect all read data.
        Also the allele names are extracted and stored in a list.
        The data content should not be changed only parsed. The number of alleles does not matter.
        """

        # Test case 1, MSA data with 5 alleles
        Input_test = SelectHybridReads.ParseInput()
        input_file = 'Header\nExtraInformation\n$$$\nReadName1\tCCCCCCCC\tIIIIIIII\nReadName2\tTTTTTTTT\tIIIIIIII\nRead1\t-----CCCCCCCC------\nRead2\t--------TTTTTTTT---\nallele_a1\t-CCCCCCCCCCCCCCC---\nallele_a2\t---CCCCCCCCCCCCCC--\nallele_b1\t--TTTTTTTTTTTT-----\nallele_b2\t----TTTTTTTTTTTTTTT\nallele_c1\tCCCCCCCCCTTTTTTTT--\n$$$\n'
        all_data, allele_names = Input_test.collect_all_data(input_file)
        self.assertEqual(all_data, [[['ReadName1', 'CCCCCCCC', 'IIIIIIII'], 
                                     ['ReadName2', 'TTTTTTTT', 'IIIIIIII'], 
                                     ['Read1', '-----CCCCCCCC------'],
                                     ['Read2', '--------TTTTTTTT---'], 
                                     ['allele_a1', '-CCCCCCCCCCCCCCC---'], 
                                     ['allele_a2', '---CCCCCCCCCCCCCC--'], 
                                     ['allele_b1', '--TTTTTTTTTTTT-----'],
                                     ['allele_b2', '----TTTTTTTTTTTTTTT'],
                                     ['allele_c1', 'CCCCCCCCCTTTTTTTT--']]])
        self.assertEqual(allele_names, ['allele_a1', 'allele_a2', 'allele_b1', 'allele_b2', 'allele_c1'])
  
        # Test case 2, MSA with 6 alleles
        Input_test = SelectHybridReads.ParseInput()
        input_file = 'Header\nExtraInformation\n$$$\nReadName1\tCCCCCCCC\tIIIIIIII\nReadName2\tTTTTTTTT\tIIIIIIII\nRead1\t-----CCCCCCCC------\nRead2\t--------TTTTTTTT---\nallele_a1\t-CCCCCCCCCCCCCCC---\nallele_a2\t---CCCCCCCCCCCCCC--\nallele_b1\t--TTTTTTTTTTTT-----\nallele_b2\t----TTTTTTTTTTTTTTT\nallele_c1\tCCCCCCCCCTTTTTTTT--\nallele_c2\tCCCTT-CCCTTTTTTTT--\n$$$\n'
        all_data, allele_names = Input_test.collect_all_data(input_file)
        self.assertEqual(all_data, [[['ReadName1', 'CCCCCCCC', 'IIIIIIII'], 
                                     ['ReadName2', 'TTTTTTTT', 'IIIIIIII'], 
                                     ['Read1', '-----CCCCCCCC------'],
                                     ['Read2', '--------TTTTTTTT---'], 
                                     ['allele_a1', '-CCCCCCCCCCCCCCC---'], 
                                     ['allele_a2', '---CCCCCCCCCCCCCC--'], 
                                     ['allele_b1', '--TTTTTTTTTTTT-----'],
                                     ['allele_b2', '----TTTTTTTTTTTTTTT'],
                                     ['allele_c1', 'CCCCCCCCCTTTTTTTT--'],
                                     ['allele_c2', 'CCCTT-CCCTTTTTTTT--']]])
        self.assertEqual(allele_names, ['allele_a1', 'allele_a2', 'allele_b1', 'allele_b2', 'allele_c1', 'allele_c2'])


    def test_get_allele_combinations(self):
        """
        All combinations of the names of the alleles must be generated. Each combination should be present once, in which the
        order does not matter. The combination should not contain the same allele name twice.
        The method can handle input with 4, 5 or 6 alleles.
        If the input contains less than 4 alleles, then a error is raised.
        """

        # Test case 1: 5 allele names
        Input_test = SelectHybridReads.ParseInput()

        allele_names = ['allele_A1', 'allele_B1', 'allele_B2', 'allele_C1', 'allele_C2']
        all_allele_combinations = Input_test.get_allele_combinations(allele_names)
        self.assertEqual(all_allele_combinations, [['allele_A1', 'allele_B1'], 
                                                   ['allele_A1', 'allele_B2'], 
                                                   ['allele_A1', 'allele_C1'],
                                                   ['allele_A1', 'allele_C2'], 
                                                   ['allele_B1', 'allele_B2'], 
                                                   ['allele_B1', 'allele_C1'], 
                                                   ['allele_B1', 'allele_C2'], 
                                                   ['allele_B2', 'allele_C1'], 
                                                   ['allele_B2', 'allele_C2'], 
                                                   ['allele_C1', 'allele_C2']])

        # Test case 2: 6 allele names
        Input_test = SelectHybridReads.ParseInput()

        allele_names = ['allele_A1', 'allele_A2', 'allele_B1', 'allele_B2', 'allele_C1', 'allele_C2']
        all_allele_combinations = Input_test.get_allele_combinations(allele_names)
        self.assertEqual(all_allele_combinations, [['allele_A1', 'allele_A2'],
                                                   ['allele_A1', 'allele_B1'],
                                                   ['allele_A1', 'allele_B2'], 
                                                   ['allele_A1', 'allele_C1'],
                                                   ['allele_A1', 'allele_C2'],
                                                   ['allele_A2', 'allele_B1'], 
                                                   ['allele_A2', 'allele_B2'], 
                                                   ['allele_A2', 'allele_C1'],
                                                   ['allele_A2', 'allele_C2'], 
                                                   ['allele_B1', 'allele_B2'],
                                                   ['allele_B1', 'allele_C1'],
                                                   ['allele_B1', 'allele_C2'],
                                                   ['allele_B2', 'allele_C1'],
                                                   ['allele_B2', 'allele_C2'], 
                                                   ['allele_C1', 'allele_C2']])
        
        # Test case 3: 4 allele names
        allele_names = ['allele_A1', 'allele_B1', 'allele_C1', 'allele_C2']
        all_allele_combinations = Input_test.get_allele_combinations(allele_names)
        self.assertEqual(all_allele_combinations, [['allele_A1', 'allele_B1'], 
                                                   ['allele_A1', 'allele_C1'],
                                                   ['allele_A1', 'allele_C2'],
                                                   ['allele_B1', 'allele_C1'],
                                                   ['allele_B1', 'allele_C2'],
                                                   ['allele_C1', 'allele_C2']])
        # Test case 4: 3 allele names
        allele_names = ['allele_A1', 'allele_B1', 'allele_C1']
        with self.assertRaises(ValueError):
            Input_test.get_allele_combinations(allele_names)
      

if __name__ == '__main__':
    unittest.main()


