"""
30-07-'19

This script contains unittests for the class ParseInput from SelectHybridReads.py.


"""

import unittest
import SelectHybridReads

class TestParseInput(unittest.TestCase):
    """
    This class contains unittests

    """

    def test_collect_all_data(self):
        """

        """
        Input_test = SelectHybridReads.ParseInput()
        input_file = 'Header\nExtraInformation\n$$$\nReadName1\tCCCCCCCC\tIIIIIIII\nReadName2\tTTTTTTTT\tIIIIIIII\nRead1\t-----CCCCCCCC------\nRead2\t--------TTTTTTTT---\nallele_a1\t-CCCCCCCCCCCCCCC---\nallele_a2\t---CCCCCCCCCCCCCC--\nallele_b1\t--TTTTTTTTTTTT-----\nallele_b2\t----TTTTTTTTTTTTTTT\nallele_c1\tCCCCCCCCCTTTTTTTT--\n$$$'
        


        all_data, allele_names = Input_test.collect_all_data(input_file)

        self.assertEqual(all_data, [[['ReadName1', 'CCCCCCCC', 'IIIIIIII'], 
                                     ['ReadName2', 'TTTTTTTT', 'IIIIIIII'], 
                                     ['Read1', '-----CCCCCCCC------'],
                                     ['Read2', '--------TTTTTTTT---'], 
                                     ['allele_a1', '-CCCCCCCCCCCCCCC---'], 
                                     ['allele_a2', '---CCCCCCCCCCCCCC--'], 
                                     ['allele_b1', '--TTTTTTTTTTTT-----'],
                                     ['allele_b2', '----TTTTTTTTTTTTTTT'],
                                     ['allele_c1', 'CCCCCCCCCTTTTTTTT--']], []])

        self.assertEqual(allele_names, ['allele_a1', 'allele_a2', 'allele_b1', 'allele_b2', 'allele_c1'])

    def test_get_allele_combinations(self):
        """

        """

        # Test case 1: 5 alleles
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

        # Test case 2: 6 alleles
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
        
        # Test case 3: 4 alleles
        allele_names = ['allele_A1', 'allele_B1', 'allele_C1', 'allele_C2']
        all_allele_combinations = Input_test.get_allele_combinations(allele_names)
        self.assertEqual(all_allele_combinations, [['allele_A1', 'allele_B1'], 
                                                   ['allele_A1', 'allele_C1'],
                                                   ['allele_A1', 'allele_C2'],
                                                   ['allele_B1', 'allele_C1'],
                                                   ['allele_B1', 'allele_C2'],
                                                   ['allele_C1', 'allele_C2']])
        # Test case 4: less than 4 alleles
        allele_names = ['allele_A1', 'allele_B1', 'allele_C1']
        with self.assertRaises(ValueError):
            Input_test.get_allele_combinations(allele_names)
      

if __name__ == '__main__':
    unittest.main()


