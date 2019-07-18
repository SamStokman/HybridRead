"""

Input needed: output file from align_all_reads.py

Script that finds hybrid reads


"""

from sys import argv
import os

class ParseInput:
    """

    Args:
        -
    """
    
    @staticmethod
    def collect_all_data(input_file):
        """
        
        Args:
            input_file (str): aligned reads and alleles, output from AllignAllRead.py
        
        Returns:
            all_data (list): list of lists, each list contains the alignment with the read pair and the best matches for HLA-A, B and C.
            allele_names (list): contains all allele types (max. 6)
        """
        all_data = []
        input_file = input_file.split('$$$')
        for read_data in input_file[1:-1]:   # 1:-1
            temp_collect_list = []
            read_data = read_data.split('\n')
            for line in read_data:
                line = line.split('\t')
                if len(line) > 1:
                    temp_collect_list += [line]
            all_data += [temp_collect_list]
        allele_names = all_data[0][4:]
        
        return all_data, allele_names

    @staticmethod
    def get_allele_combinations(allele_names):
        """
        This function creates all allele combinations (the type names) and returns them in a list
        
        Args:
            allele_names (list): list with all allele types (max. 6)
        
        Returns:
            all_combinations_list (list): contains all possible allele combinations
        """
        # for samples with 5 different alleles
        if len(allele_names) == 5:
            first_nr_list = [0,0,0,0,1,1,1,2,2,3]
            second_nr_list = [1,2,3,4,2,3,4,3,4,4] 
        
        # for samples with 6 different alleles
        if len(allele_names) == 6:
            first_nr_list = [0,0,0,0,0,1,1,1,1,2,2,2,3,3,4]
            second_nr_list = [1,2,3,4,5,2,3,4,5,3,4,5,4,5,5]

        all_alleles_list = []
        for allele, seq in allele_names:
           all_alleles_list += [allele]

        all_combinations_list = []
        for i in range(len(first_nr_list)):
            all_combinations_list += [[all_alleles_list[int(first_nr_list[i])], all_alleles_list[int(second_nr_list[i])]]]
        
        return (all_combinations_list)
        
class Read:
    """
    
    Args:
        read_seq (str):
        read_qv (str):
        read_aligned_seq (str):
        allele_data (list):
    """
    def __init__(self, read_seq, read_qv, read_aligned_seq, allele_data):
        self.read_seq = read_seq
        self.read_length = len(read_seq)
        self.read_qv = read_qv
        self.read_aligned_seq = read_aligned_seq
        self.allele_data = allele_data
    
    
    def check_alignment(self): 
        """

        Args:
            -
        Returns:
            correct_alignment (bool):
        """
        # Check if alignment is correct
        
        correct_alignment = True
        check_alignment = self.read_aligned_seq.replace('-', '')
        if self.read_seq != check_alignment:
            correct_alignment = False
        
        return correct_alignment
    
    def apply_qv(self):
        """
        Checks nucleotide quality values, if lower than a given value, then the nuclotide is replaced by a 'N'
        
        Args:
            -
        Returns:
           read_checked_aligned_seq (str):
        """
        quality_dict = {'!': 0, '"': 1, '#':2, '$':3, '%':4, '&':5, "'":6, '(':7, ')':8,\
        '*':9, '+':10, ',':11,'-':12, '.':13, '/':14, '0':15, '1':16, '2':17, '3':18,\
        '4':19, '5':20, '6':21, '7':22, '8':23, '9':24, ':':25, ';':26, '<':27, '=':28,\
        '>':29, '?':30, '@':31, 'A':32, 'B':33, 'C':34, 'D':35, 'E':36, 'F':37, 'G':38, 'H':39, 'I':40}
    
        read_checked_seq = '' 
        for i, qual in enumerate(self.read_qv):
            q_score = quality_dict[qual]
            nucleotide = self.read_seq[i]
            if q_score < 18:
                read_checked_seq += 'N'
            if q_score >= 18:
                read_checked_seq += nucleotide

        # count the number of '-' in front of aligned read (left)
        count_read_start = self.read_aligned_seq.lstrip('-')
        count_read_start = len(self.read_aligned_seq) - len(count_read_start)

        # remove '-' right from aligned read and count the number of '-' after the aligned read (right)
        read_aligned_seq_wo_r = self.read_aligned_seq.rstrip('-')
        count_read_end = len(self.read_aligned_seq) - len(read_aligned_seq_wo_r)

        read_checked_aligned_seq = ''
        read_start_seen = False   
        gap_count = 0
        for i, char in enumerate(read_aligned_seq_wo_r):
            if char == '-':
                read_checked_aligned_seq += char
            if char != '-' and read_start_seen == False:
                read_start_seen = True
            if read_start_seen == True:
                if char == '-':
                    gap_count += 1
                if char != '-':
                    read_checked_aligned_seq += read_checked_seq[i-count_read_start-gap_count]

        # add '-' at the right
        read_checked_aligned_seq = read_checked_aligned_seq + '-' * count_read_end
        return read_checked_aligned_seq
    

    def check_read_artefacts(self, read_aligned_qv):
        """" 
        Check if read has artefect
        Artefact definition: if all alleles have have a mismatch at the same position 
        If artefact found, read nucleotide replaced by 'N'

        Args:
            read_aligned_qv (str):
        Returns:
            read_aligned_fully_checked (str):
        """

        read_seq = read_aligned_qv

        # get the start of the read position (relative to the alleles) 
        start_relative_read_position = []
        start_relative_read_nucleotide = []
        nr_of_switches = 0
        nuc = read_aligned_qv[0]
        for i, char in enumerate(read_aligned_qv):
            if nuc == '-' and nuc != char:
                nr_of_switches += 1
                start_relative_read_position += [i]
                start_relative_read_nucleotide += [char]
            nuc = char
    
        # get all read positions (relative to the alleles)
        relative_read_position = []
        relative_read_nucleotide = []
        for i, char in enumerate(read_seq):
            if i >= start_relative_read_position[0] and i <= start_relative_read_position[-1]:
                relative_read_position += [i]
                relative_read_nucleotide += [char]
            if i >= start_relative_read_position[-1]+1 and char != '-':
                relative_read_position += [i]
                relative_read_nucleotide += [char]
       
        mismatch_track = self.__create_mismatch_track(allele_data, read_seq, relative_read_position, relative_read_nucleotide)

        # if a nucleotide in the read is considered as artefact (5 in mismatch track) then it is replaced by 'N'
        read_aligned_fully_checked = ''
        for i, nuc in enumerate(read_seq):
            try:
                mismatch_char = int(mismatch_track[i])
            except:
                mismatch_char = mismatch_track[i]
            if mismatch_char == 5:
                read_aligned_fully_checked += 'N'
            if mismatch_char != 5:
                read_aligned_fully_checked += nuc
        return read_aligned_fully_checked

    @staticmethod
    def __create_mismatch_track(allele_data, read_seq, relative_read_position, relative_read_nucleotide):
        """
        
        
        Args:
            allele_data
            read_seq
            relative_read_position
            relative_read_nucleotide
        Returns:
            mismatch_track:
        """
        # create a mismatch track string, for each mismatch in the allele '1' is added, up to 5  (where all alleles have mismatches)
        mismatch_track = '-' * len(read_seq)

        for allele, seq in allele_data:
            seq_string = seq
            for i, chari in enumerate(seq_string):
                if i in relative_read_position:
                    read_nuc = relative_read_nucleotide[i-min(relative_read_position)] 
                    if chari != read_nuc and read_nuc != '-' and read_nuc != 'N' and chari != '-' and  mismatch_track[i] == '4':
                        left_mismatch_track = mismatch_track[:i]
                        right_mismatch_track = mismatch_track[i+1:]
                        mismatch_track = left_mismatch_track + '5' + right_mismatch_track
                    if chari != read_nuc and read_nuc != '-' and read_nuc != 'N' and chari != '-' and  mismatch_track[i] == '3':
                        left_mismatch_track = mismatch_track[:i]
                        right_mismatch_track = mismatch_track[i+1:]
                        mismatch_track = left_mismatch_track + '4' + right_mismatch_track
                    if chari != read_nuc and read_nuc != '-' and read_nuc != 'N' and chari != '-' and  mismatch_track[i] == '2':
                        left_mismatch_track = mismatch_track[:i]
                        right_mismatch_track = mismatch_track[i+1:]
                        mismatch_track = left_mismatch_track + '3' + right_mismatch_track
                    if chari != read_nuc and read_nuc != '-' and read_nuc != 'N' and chari != '-' and  mismatch_track[i] == '1':
                        left_mismatch_track = mismatch_track[:i]
                        right_mismatch_track = mismatch_track[i+1:]
                        mismatch_track = left_mismatch_track + '2' + right_mismatch_track
                    if chari != read_nuc and read_nuc != '-' and read_nuc != 'N' and chari != '-' and mismatch_track[i] == '-':
                        left_mismatch_track = mismatch_track[:i]
                        right_mismatch_track = mismatch_track[i+1:]
                        mismatch_track = left_mismatch_track + '1' + right_mismatch_track

        return mismatch_track

    def get_mismatches(self, read_aligned_fully_checked):
        """
        Check for mismatches read vs all alleles (substitutions, insertions and deletions)
        
        Args:
            allele_data: list, with all alleles and sequences
            read_seq: str, read sequence
        Returns:
           mismatch_dict (dict):
           extended_mismatch_dict (dict):
        """
        read_seq = read_aligned_fully_checked

        # get the read position (relative to the alleles) 
        start_relative_read_position = []
        start_relative_read_nucleotide = []
        nr_of_switches = 0
        nuc = read_aligned_fully_checked[0]
        for i, char in enumerate(read_aligned_fully_checked):
            if nuc == '-' and nuc != char:
                nr_of_switches += 1
                start_relative_read_position += [i]
                start_relative_read_nucleotide += [char]
            nuc = char

        # get the allele position (relative to the alleles), to check whether the read starts in front of the allele
        deletion_correction_dict = self.__check_read_start(start_relative_read_position, read_aligned_fully_checked)

        relative_read_position = []
        relative_read_nucleotide = []
        for i, char in enumerate(read_seq):
            if i >= start_relative_read_position[0] and i <= start_relative_read_position[-1]:
                relative_read_position += [i]
                relative_read_nucleotide += [char]
            if i >= start_relative_read_position[-1]+1 and char != '-':
                relative_read_position += [i]
                relative_read_nucleotide += [char]

        extended_mismatch_dict = {}
        mismatch_dict = {}
        for allele, seq_string in self.allele_data:
            substitutions = 0
            insertions = 0
            deletions = int(deletion_correction_dict[allele])
            mismatches = 0
            for i, chari in enumerate(seq_string):
                if i in relative_read_position:
                    read_nuc = relative_read_nucleotide[i-min(relative_read_position)]
                    if chari != read_nuc:
                        if read_nuc != '-' and read_nuc != '*' and read_nuc != 'N' and chari != '-': 
                            substitutions += 1    
                        if read_nuc == '-' and chari != '-':
                            insertions += 1
                        if read_nuc != '-' and read_nuc != '*' and chari == '-':
                            deletions += 1
                        mismatches = substitutions + insertions + deletions
            mismatch_dict[allele] = [mismatches]
            extended_mismatch_dict[allele] = [substitutions, insertions, deletions, mismatches]
        return (mismatch_dict, extended_mismatch_dict)
  


    def __check_read_start(self, start_relative_read_position, read_aligned_fully_checked):
        """
       
        Args:
            start_relative_read_position (list):
            read_aligned_fully_checked (str):
        Returns:
            deletion_correction_dict (dict):                    
        """
        deletion_correction_dict = {} # if reads starts in front of allele
        for allele, seq_string in self.allele_data: 
            start_allele_post = 0
            if seq_string.startswith('-'):
                start_relative_allele_position = []
                start_relative_allele_nucleotide = []
                nuc = seq_string[0]
                nr_of_switches = 0
                for i, char in enumerate(seq_string):
                    if nuc == '-' and nuc != char:
                        nr_of_switches += 1
                        start_relative_allele_position += [i]
                        start_relative_allele_nucleotide += [char]
                    nuc = char
                start_allele_post = start_relative_allele_position[0]
                if start_allele_post > start_relative_read_position[0]:
                    deletion_correction = start_relative_read_position[0] - start_allele_post
                    deletion_correction_dict[allele] = deletion_correction
                else: 
                    deletion_correction = 0
                    deletion_correction_dict[allele] = deletion_correction
            else: 
                deletion_correction = 0
                deletion_correction_dict[allele] = deletion_correction
            count_read_insertion_for_deletion_correction = read_aligned_fully_checked[start_relative_read_position[0]:start_allele_post]
            read_inserts = count_read_insertion_for_deletion_correction.count('-')
            deletion_correction_dict[allele] += read_inserts

        return deletion_correction_dict
   
    @classmethod
    def classmethod_for_non_read(cls, aligned_sequence, allele_data):
        """
        Classmethod that generates data for the read consensus and turnover region (both in alignment) as input for the constructor. 
        The sequence (read_seq) is created  without the alignment. The read_qv is an empty string since the read_consensus and turnover
        region do not have quality values.
        
        Args:
            aligned_sequence (str):
            allele_data (list):
        Returns:
            read_seq (str):
            read_qv (str): Empty string
            aligned_sequence (str):
            allele_data (list):
        """
        read_seq = aligned_sequence.lstrip('-').rstrip('-')
        read_qv = ''

        return cls(read_seq, read_qv, aligned_sequence, allele_data)
    
    def print_mismatches(self, read_type, extended_mismatch_dict):
        """
        Print type of read, read length and the number mismatches
        
        Args:
            read_type (str): Discribes read type; First read, Second read or Read consensus
            extended_mismatch_dict (dict): Contains the number of SNP substitutions, insertions, deletions and total number of mismatches
        Returns:
           -
                    
        """
        print ('\n\nRead info', read_type)
        print ('Length: \t\t', self.read_length)
        
        print ('\nAllele\t\t\tSubstitutions\tInsertions\tDeletions\tTotal nr. or mismatches')
        for allele, mismatches in extended_mismatch_dict.items():
            print (allele, '\t\t', mismatches[0], '\t\t', mismatches[1], '\t\t', mismatches[2], '\t\t', mismatches[3])

    def get_relative_position(self):
        """
        Determines the (consensus) read position for each allele and returns it in 
        a dict.        
        
        Args:
            -            
        Returns:
            read_pos_dict (dict):  
        """
        read_pos_dict = {}
        for allele, allele_seq in self.allele_data:
            read_position = []
            # remove '-' left and right from allele seq and give read the same length
            allele_seq_wo_left = allele_seq.lstrip('-')
            left_difference = len(allele_seq) - len(allele_seq_wo_left)
            allele_seq_wo_right = allele_seq.rstrip('-')
            sequence_seq_wo_rl = self.read_aligned_seq[left_difference:len(allele_seq_wo_right)]
            allele_seq_wo_rl = allele_seq.strip('-')

            # Remove '-' right from the read
            sequence_read = sequence_seq_wo_rl.rstrip('-')

            read_start_seen = False
            deletion_count_allele_to_read_start = 0

            for i, char in enumerate(sequence_read):
                if char == '-' and read_start_seen ==  False:
                    if allele_seq_wo_rl[i] == '-':
                        deletion_count_allele_to_read_start += 1
                if char != '-':
                    read_start_seen = True

            read_start_seen = False
            pos = 0
        
            for i, char in enumerate(sequence_read):
                if char != '-' and allele_seq_wo_rl[i] != '-' and read_start_seen == False:
                    read_start_seen = True
                    pos = i - deletion_count_allele_to_read_start
                if char != '-' and allele_seq_wo_rl[i] != '-' and read_start_seen == True:
                    read_position += [pos]
                if char == '-' and allele_seq_wo_rl[i] != '-' and read_start_seen == True: # if read has deletion, allele position is taken into account
                    read_position += [pos]
                if len(read_position) != 0:
                    pos = read_position[-1] + 1
            read_pos_dict[allele] = [read_position]

        # get position if turnover region has a length of 0 and the allele has '-' as nucleotide
        allele_name =  self.allele_data[0][0]

        if 'K' in self.read_aligned_seq and read_pos_dict[allele_name] == [[]]:
            read_pos_dict = self.__get_special_case_pos(read_pos_dict, allele_name)

        return read_pos_dict


    def __get_special_case_pos(self, read_pos_dict, allele_name):  #NOG IETS VERZINNEN VOOR LENGTH IS 1
        """
        If turnover region has a length of 0 and the allele has '-' as nucleotide, still get the correct position
        
        Args:
            read_pos_dict (dict):
            allele_name (str):
        Returns:
            read_pos_dict (dict): an updated version of original read_pos_dict
        """

        read_start_seen = False
        pos = 0
        read_position = []
        for allele, seq in  self.allele_data:
            allele_seq = seq
            
            # remove '-' left and right from allele seq and give read the same length
            allele_seq_wo_left = allele_seq.lstrip('-')
            left_difference = len(allele_seq) - len(allele_seq_wo_left)
            allele_seq_wo_right = allele_seq.rstrip('-')
            sequence_seq_wo_rl = self.read_aligned_seq[left_difference:len(allele_seq_wo_right)]
            allele_seq_wo_rl = allele_seq.strip('-')

            # Remove '-' right from the read
            sequence_read = sequence_seq_wo_rl.rstrip('-')

            read_start_seen = False
            deletion_count_allele_to_read_start = 0
            for i, char in enumerate(sequence_read):
                if char == '-' and read_start_seen ==  False:
                    if allele_seq_wo_rl[i] == '-':
                        deletion_count_allele_to_read_start += 1
                if char != '-':
                    read_start_seen = True

            read_start_seen = False
            pos = 0
            for i, char in enumerate(sequence_read):
                pos = i - deletion_count_allele_to_read_start
                if char != '-' and read_start_seen == False:
                    read_position = [pos]

            del read_pos_dict[allele_name]
            read_pos_dict[allele_name] = [read_position]

        return (read_pos_dict)

        
class ReadPair():
    """
    
    Args:
        read1_aligned_checked (str):
        read2_aligned_checked (str):
        read1_seq (str):
        read2_seq (str):
        min_read_length (int):
        N_quantity (int):
    """
    
    min_read_length = 0
    N_quantity = 50000

    def __init__(self, read1_aligned_checked, read2_aligned_checked, read1_seq, read2_seq):
        self.read1_aligned_checked = read1_aligned_checked
        self.read2_aligned_checked = read2_aligned_checked
        self.read1_seq = read1_seq
        self.read2_seq = read2_seq

        
    def check_read_pair(self):
        """
        The minimum read length of the original reads and the number of allowed N's per read can be determined
        
        Args:
            -
            
        Returns:
            approve_reads (bool):
        """

        approve_reads = True

        read1_length = len(self.read1_seq)
        read2_length = len(self.read2_seq)
        read1_N_count = self.read1_aligned_checked.count('N')
        read2_N_count = self.read2_aligned_checked.count('N')

        # check read length
        if read1_length < self.min_read_length or read2_length < self.min_read_length:    
            approve_reads = False

        # check nr of Ns
        if read1_N_count > self.N_quantity or read2_N_count > self.N_quantity:    
            approve_reads = False

        return approve_reads

   
    def create_read_consensus(self):
        """
        This function combines the two paired-end reads into a 'read consensus'
        if the reads do not overlap (and thus have a gap) the empty nucleotide between the 
        reads are replaced by '*'. This to distinguish between real read deletions and
        the empty space between the reads. For overlapping reads: tf the reads have different nucleotides 
        at the same position then it is replaced by a 'N'. If 1 read has a 'N' and the other one a 
        nucleotide, the nucleotide is used.
        
        Args:
            -
            
        Returns:
            read_consensus (str): 
        """
        read1_seq = self.read1_aligned_checked
        read2_seq = self.read2_aligned_checked

        
        # replace '-' in front of and after read 1 with '*', but not in the read itself 
        read1_seq_wo_l = read1_seq.lstrip('-')
        length_left = len(read1_seq) - len(read1_seq_wo_l)
        read1_seq_temp = length_left*'*' + read1_seq_wo_l
        read1_seq_temp_wo_r = read1_seq_temp.rstrip('-')
        length_right = len(read1_seq) - len(read1_seq_temp_wo_r)
        read1_seq_temp = read1_seq_temp_wo_r + length_right*'*'

        # replace '-' in front of and after read 2 with '*', but not in the read itself 
        read2_seq_wo_l = read2_seq.lstrip('-')
        length_left = len(read2_seq) - len(read2_seq_wo_l)
        read2_seq_temp = length_left*'*' + read2_seq_wo_l
        read2_seq_temp_wo_r = read2_seq_temp.rstrip('-')
        length_right = len(read2_seq) - len(read2_seq_temp_wo_r)
        read2_seq_temp = read2_seq_temp_wo_r + length_right*'*'

        # Combine sequences read 1 and 2 (consensus)
        read_consensus_temp = ''
        for i in range(len(read1_seq)):
            if read1_seq_temp[i] == '*' and read2_seq_temp[i] == '*':
                read_consensus_temp += read1_seq_temp[i]
            if read1_seq_temp[i] != '*' and read2_seq_temp[i] == '*':
                read_consensus_temp += read1_seq_temp[i]
            if read1_seq_temp[i] == '*' and read2_seq_temp[i] != '*':
                read_consensus_temp += read2_seq_temp[i]
            if read1_seq_temp[i] != '*' and read2_seq_temp[i] != '*':
               if read1_seq_temp[i] == read2_seq_temp[i]:
                   read_consensus_temp += read1_seq_temp[i]
               if read1_seq_temp[i] == 'N' and read2_seq_temp[i] != 'N':
                   read_consensus_temp += read2_seq_temp[i]
               if read2_seq_temp[i] == 'N' and read1_seq_temp[i] != 'N':
                   read_consensus_temp += read1_seq_temp[i]
               if read1_seq_temp[i] != 'N' and read2_seq_temp[i] != 'N' and read1_seq_temp[i] != read2_seq_temp[i]:
                   read_consensus_temp += 'N'

        # replace '*' in front of and after read consensus with '-', but not between the reads
        read_consensus = ''
        read_consensus_wo_l = read_consensus_temp.lstrip('*')
        length_left = len(read_consensus_temp) - len(read_consensus_wo_l)
        read_consensus_temp = length_left*'-' + read_consensus_wo_l
        read_consensus_temp_wo_r = read_consensus_temp.rstrip('*')
        length_right = len(read_consensus_temp) - len(read_consensus_temp_wo_r)
        read_consensus = read_consensus_temp_wo_r + length_right*'-'
        
        return (read_consensus)

        
class CheckAlleleCombination():
    """
    Args:
        read_consensus (str):
        allele_combo (list): contains the two allele names
        allele_data (list):
    """
    
    def __init__(self, read_consensus, allele_combo, allele_data):
        self.read_consensus = read_consensus
        self.allele_combo = allele_combo
        self.allele1 = allele_combo[0]
        self.allele2 = allele_combo[1]
        self.allele_data = allele_data
        self.indicator_string = ''
        self.number_of_artefacts = 0

    
    def create_indicator_string(self):
        """
        
        Args:
            -
        
        Returns:
            allele_seq_list (list):
        """

        read_consensus = self.read_consensus
        allele_combo = self.allele_combo
        allele_info = self.allele_data

        allele_match_list = allele_combo

        # Create a dict with only the allele matches with the lowest number of mismatches. The sequences now contain '-' for matches and 'X' or 'Y' for a mismatch (mismatch indicator sequences).

        turn_seq_list = []
        allele_seq_list = []
        mismatch_char = 'X'
   
        start_consensus = read_consensus.rstrip('-').count('-') # start read pos (absolute)
        end_consensus = len(read_consensus.rstrip('-'))  # end read pos (absolute)

        for allele, seq_string in allele_info:
            mismatches = 0
            seq_string_temp = ''
            turn_over_seq_temp = ''

            if allele in allele_match_list:
                allele_seq_list += [seq_string]
            if allele in allele_match_list:
                for i, chari in enumerate(seq_string):
                    if chari != read_consensus[i]:
                        if read_consensus[i] != '-' and chari != '-' and read_consensus[i] != '*' and read_consensus[i] != 'N':
                            mismatches += 1    
                            turn_over_seq_temp += mismatch_char
                        if read_consensus[i] == '-' and chari != '-':  
                            if i < start_consensus or i >= end_consensus: # nucleotide in front of and after reads
                                turn_over_seq_temp += '-'
                            if i >= start_consensus and i < end_consensus:  # allele insertion
                                turn_over_seq_temp += mismatch_char
                        if read_consensus[i] != '-' and chari == '-' and read_consensus[i] != '*' and read_consensus[i] != 'N':   # allele deletion
                            turn_over_seq_temp += mismatch_char
                        if read_consensus[i] == '*' or read_consensus[i] == 'N':
                            turn_over_seq_temp += '-'
                    if chari == read_consensus[i]: 
                        seq_string_temp += chari
                        turn_over_seq_temp += '-'
                mismatch_char = 'Y'
                turn_seq_list += [[allele,turn_over_seq_temp]]

        # Extract the mismatch indicator sequences for allele match 1 and 2
        seq_allele_1 = turn_seq_list[0][1]
        seq_allele_2 = turn_seq_list[1][1]

        # Str with mismatch indicator sequences for allele match 1 and 2 combined
        mismatch_indicator_combo_str = ''
        for i in range(len(seq_allele_1)):
            if seq_allele_1[i] == seq_allele_2[i] and seq_allele_1[i] != 'N':
                mismatch_indicator_combo_str += seq_allele_1[i]
            if seq_allele_1[i] != seq_allele_2[i]:
                if seq_allele_1[i] == '-' and seq_allele_2[i] != '-':
                    mismatch_indicator_combo_str += seq_allele_2[i]
                if seq_allele_1[i] != '-' and seq_allele_2[i] == '-':
                    mismatch_indicator_combo_str += seq_allele_1[i]
                if seq_allele_1[i] != '-' and seq_allele_2[i] != '-': 
                    mismatch_indicator_combo_str += 'M'
        
        self.indicator_string = mismatch_indicator_combo_str

        return allele_seq_list
    
    def check_indicative_SNPs(self):
        """
        check if alleles have enough indicative mismatches (at least 2 per allele)
        
        Args:
            -
        
        Returns:
            accept_combo (bool):
        """
        accept_combo = True
        
        mismatch_indicator_string = self.indicator_string
        if mismatch_indicator_string.count('X') < 2 or mismatch_indicator_string.count('Y') < 2:
            accept_combo = False

        return accept_combo

    def check_mutual_SNPs(self):
        """
        check if alleles do not have too many mutual mismatches (max. 2)
                
        Args:
            -
        
        Returns:
            accept_combo (bool):
        """
        accept_combo = True
        mismatch_indicator_string = self.indicator_string
        self.number_of_artefacts += mismatch_indicator_string.count('M') 

        if mismatch_indicator_string.count('M') > 2:
            accept_combo = False

        return accept_combo
        
    def check_alternately_SNPs(self):
        """
        check if alleles do not have too many alternately [XYX or YXY] mismatches (max. 2)
                
        Args:
            -
        
        Returns:
            count_indicator_list (list):
            number_of_artefacts (int):
        """
        mismatch_indicator_string = self.indicator_string

        # replace 'M' for '-'  if the alleles have a mutual mismatch, it is ignored and it is regarded as a read artefact
        mismatch_indicator_string = mismatch_indicator_string.replace('M', '-')

        # to check for XYX situations
        only_mismatch_indicator_chars = mismatch_indicator_string.replace('-', '')
        start_allele = only_mismatch_indicator_chars[0]
        count_indicator_length = 0
        count_indicator_list = []
        for char in mismatch_indicator_string:
            if char != '-' and char == start_allele:
                count_indicator_length += 1
            if char != '-' and char != start_allele:
                count_indicator_length = 1
            if char != '-':
                count_indicator_list += [count_indicator_length]
                start_allele = char

        # The number of XYX situations, only a single char between two others
        pcr_artefact = 0

        for i, number in enumerate(count_indicator_list):
            if i < len(count_indicator_list)-1:
                next_number = count_indicator_list[i+1]
                if number == 1 and next_number == 1:
                    pcr_artefact += 1
        self.number_of_artefacts += pcr_artefact
        # check for allele artefact, XYX or YXY (max. 2) 
        if pcr_artefact >= 3:
            count_indicator_list = None

        return count_indicator_list, self.number_of_artefacts
        
    def update_indicator_string(self, count_indicator_list):
        """
        Remove alternately artefacts
        
        Args:
            count_indicator_list (list): 
            
        Returns:
            final_indicator_string (str):          
        """
        mismatch_indicator_string = self.indicator_string
        mismatch_indicator_string = mismatch_indicator_string.replace('M', '-')
        

        indicator_combo_str_wo_artefacts = ''
        only_mismatch_indicator_chars = mismatch_indicator_string.replace('-', '')
        new_indicator_str = ''
        
        for i, number in enumerate(count_indicator_list):
            if i < len(count_indicator_list)-1:
                next_number = count_indicator_list[i+1]
                if number != 1 or next_number != 1:
                    new_indicator_str += only_mismatch_indicator_chars[i]
                if number == 1 and next_number == 1:
                    new_indicator_str += '-'
        
        new_indicator_str += only_mismatch_indicator_chars[-1]
        for i, char in enumerate(mismatch_indicator_string):
            if char == '-':
                indicator_combo_str_wo_artefacts += char
            if char != '-':
               if len(new_indicator_str) >= 1:
                    replace_char = new_indicator_str[0]
                    indicator_combo_str_wo_artefacts += replace_char
                    new_indicator_str = new_indicator_str[1:]

        final_indicator_string = indicator_combo_str_wo_artefacts

        return final_indicator_string
        
    def get_switches(self, final_indicator_string):
        """
        
        Args:
            final_indicator_string (str):
        
        Returns:
            nr_of_switches (int):
            start_turn_pos (int):
            end_turn_pos (int):
        """

        end_turn_pos = 0
        nr_of_switches = 0
        only_mismatch_indicator_chars = final_indicator_string.replace('-', '')
        switch_char = only_mismatch_indicator_chars[0]
        start_char = only_mismatch_indicator_chars[0]


        switch_seen = False
        for i, char in enumerate(final_indicator_string):
            if char != '-' and char != switch_char:
                nr_of_switches += 1
                switch_char = char
            if char == start_char:
                start_turn_pos = i + 1 # the first nucleotide after X or Y
            if char != '-' and char != start_char and switch_seen == False:
                end_turn_pos = i - 1  # the first nucleotide in front of X or Y
                switch_seen = True
                
        return (nr_of_switches, start_turn_pos, end_turn_pos)
    
    def print_1_switch_alleles(self):
        """
        Prints allele combination if they resulted in 1 switch hybrid read and the number of artefacts (max. 4)
        
        Args:
            -
        
        Returns:
            -
        
        """
        print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
        print ('Allele combination:\t\t ', self.allele1.strip(' '), 'and', self.allele2)
        print ('Number of artefacts:\t\t ', self.number_of_artefacts, '\n')

class GetOneSwitchData(): 
    """
    
    Args:
        allele1 (str): name allele 1
        allele2 (str): name allele 2
        
    """    

    def __init__(self, allele1, allele2):
        self.allele1 = allele1
        self.allele2 = allele2
         
    def get_read_position(self, read1_pos_dict, read2_pos_dict):
        """
        Extracts the first and last value from read 1 and read 2 for both alleles positions
        
        Args:
            read1_pos_dict (dict):
            read2_pos_dict (dict):
            
        Returns:
            position_read1_allele1 (str):
            position_read2_allele1 (str):
            position_read1_allele2 (str):
            position_read2_allele2 (str):       
        """

        # Extract turnover sequence, based on start and end position and the sequence of allele 1 and allele 2 
        allele1 = self.allele1
        allele2 = self.allele2
        
        # Get all start and stop positions from best allele matches (reads)
        positions_list_allele1_read1 = read1_pos_dict[allele1]
        start_pos_allele1_read1 = positions_list_allele1_read1[0][0]
        end_pos_allele1_read1 = positions_list_allele1_read1[0][-1]

        positions_list_allele2_read1  = read1_pos_dict[allele2]
        start_pos_allele2_read1 = positions_list_allele2_read1[0][0]
        end_pos_allele2_read1 = positions_list_allele2_read1[0][-1] 

        positions_list_allele1_read2 = read2_pos_dict[allele1]
        start_pos_allele1_read2 = positions_list_allele1_read2[0][0]
        end_pos_allele1_read2 = positions_list_allele1_read2[0][-1] 

        positions_list_allele2_read2 = read2_pos_dict[allele2]
        start_pos_allele2_read2 = positions_list_allele2_read2[0][0]
        end_pos_allele2_read2 = positions_list_allele2_read2[0][-1] 

        position_read1_allele1 = str(start_pos_allele1_read1) + '-' +  str(end_pos_allele1_read1)
        position_read2_allele1 = str(start_pos_allele1_read2) + '-' +  str(end_pos_allele1_read2)
        position_read1_allele2 = str(start_pos_allele2_read1) + '-' +  str(end_pos_allele2_read1)
        position_read2_allele2 = str(start_pos_allele2_read2) + '-' +  str(end_pos_allele2_read2)
        
        return (position_read1_allele1, position_read2_allele1, position_read1_allele2, position_read2_allele2)

    def prep_for_turnover_position(self, start_pos, end_pos, allele_seq_list, allele_data): 
        """
        
        Args:
            start_pos (int):
            end_pos (int):
            allele_seq_list (list):
            allele_data (list):
            
        Returns:
            turn_over_region1_for_pos(str): turnover sequence region in alignment
            seq_list_allele1 (list): allele name and aligned sequence for first allele in allele combination
            turn_over_region2_for_pos(str): turnover sequence region in alignment
            seq_list_allele2 (list): allele name and aligned sequence for second allele in allele combination
        
        """
    
        allele1 = self.allele1
        allele2 = self.allele2

        turn_over_region1_for_pos = self.__get_pos_TO_region(allele_seq_list[0], start_pos, end_pos)
        turn_over_region2_for_pos = self.__get_pos_TO_region(allele_seq_list[1], start_pos, end_pos)
 
        if start_pos == end_pos +1:
            print ('Turnover sequence contains 0 nucleotides')

        # Extract allele sequences
        seq_dict_allele1 = {}
        seq_dict_allele2 = {}
        for allele, seq in allele_data:
            if allele1 == allele:
                seq_dict_allele1[allele] = seq
            if allele2 == allele:
                seq_dict_allele2[allele] = seq
      
        seq_list_allele1 = sorted(seq_dict_allele1.items())
        seq_list_allele2 = sorted(seq_dict_allele2.items())

        return turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2

    @staticmethod
    def __get_pos_TO_region(aligned_allele, start_pos, end_pos):
        """
        
        Args:
            aligned_allele (str):
            start_pos (int):
            end_pos (int):
            
        Returns:    
            turn_over_region_for_pos (str):
            
        
        """
        turn_over_region_for_pos = ''
        for i, char in enumerate(aligned_allele):
            if start_pos == end_pos:    # for TO with length 1
                if i != start_pos:
                    turn_over_region_for_pos += '-'
                if i == start_pos:
                    if char == '-':   #nog niet tegen gekomen
                        print ('heeerreeee')
                        quit()
                        char = 'Z'
                    turn_over_region_for_pos += char
            if start_pos == end_pos + 1:    # for TO with length 0
                char = 'K'     # maken we gewoon altijd een K van, deze sequence bestaat toch niet
                if i != start_pos:
                    turn_over_region_for_pos += '-'
                if i == start_pos:
                    turn_over_region_for_pos += char
            if start_pos != end_pos + 1 and start_pos != end_pos:  # for TO > length 1
                if i < start_pos:
                    turn_over_region_for_pos += '-'
                if i > end_pos:
                    turn_over_region_for_pos += '-'
                if i >= start_pos and i <= end_pos:
                    turn_over_region_for_pos += char
        return (turn_over_region_for_pos)
       
    def get_TO_position(self, TO_allele1_dict, TO_allele2_dict, turn_over_region1_for_pos, turn_over_region2_for_pos):
        """
        
        Args:
            TO_allele1_dict (dict):
            TO_allele2_dict (dict):
            turn_over_region1_for_pos (str):
            turn_over_region2_for_pos (str):
       
        Returns:
            position_to_region1 (list?):
            position_to_region2
            turn_over_region1 (str):
            turn_over_region2 (str):
        """

        allele1 = self.allele1
        allele2 = self.allele2
        # Get all start and stop positions from best allele matches (TO regions)
        # if to region 1 has 0 nucleotides
        if TO_allele1_dict[allele1] == [[]]:
            start_pos_allele1_TO_region = 0
            end_pos_allele1_TO_region = 0
        if TO_allele1_dict[allele1] != [[]]:
            positions_list_allele1_TO_region = TO_allele1_dict[allele1]
            start_pos_allele1_TO_region = positions_list_allele1_TO_region[0][0]
            end_pos_allele1_TO_region = positions_list_allele1_TO_region[0][-1]
        # if to region 2 has 0 nucleotides
        if TO_allele2_dict[allele2] == [[]]:
            start_pos_allele2_TO_region = 0
            end_pos_allele2_TO_region = 0
        if TO_allele2_dict[allele2] != [[]]:
            positions_list_allele2_TO_region  = TO_allele2_dict[allele2]
            start_pos_allele2_TO_region = positions_list_allele2_TO_region[0][0]
            end_pos_allele2_TO_region = positions_list_allele2_TO_region[0][-1] 

        # delete deletions (just keep the sequence), and K, since it is a non existing TO sequence (lenght 0)
        turn_over_region1 = turn_over_region1_for_pos.replace('K','-').replace('-','')
        turn_over_region2 = turn_over_region2_for_pos.replace('K','-').replace('-','')

        # if length TO position is 0 or 1, just the start position is taken into account
        if start_pos_allele1_TO_region != end_pos_allele1_TO_region:
            position_to_region1 = str(start_pos_allele1_TO_region) + '-'+ str(end_pos_allele1_TO_region)
        if start_pos_allele2_TO_region != end_pos_allele2_TO_region:
            position_to_region2 = str(start_pos_allele2_TO_region) + '-'+ str(end_pos_allele2_TO_region)
        if start_pos_allele1_TO_region == end_pos_allele1_TO_region:
            position_to_region1 = str(start_pos_allele1_TO_region)
        if start_pos_allele2_TO_region == end_pos_allele2_TO_region:
            position_to_region2 = str(start_pos_allele2_TO_region)
       
        return (position_to_region1, position_to_region2, turn_over_region1, turn_over_region2)
    
    @staticmethod
    def print_TO_output(allele_name, pos_read1_allele, pos_read2_allele, turn_over_region, pos_to_region):
        """
        Prints all 1 switch hybrid read data (per allele).
        
        Args:
            allele_name (str): Allele name
            pos_read1_allele (str): Start and stop position for read 1 relative to the allele
            pos_read2_allele (str): Start and stop position for read 2 relative to the allele
            turn_over_region (str): Sequence of turnover region
            pos_to_region (str): Start and stop position of turnover region relative to the allele
        
        Returns:
            -        
        """

        print ('Allele match:\t\t\t ', allele_name)
        print ('Read 1 position:\t\t ', pos_read1_allele)
        print ('Read 2 position:\t\t ', pos_read2_allele)
        print ('Turnover sequence length:\t ', len(turn_over_region))
        if turn_over_region == '':
            turn_over_region = '-'
        print ('Turnover sequence:\t\t ', turn_over_region)
        print ('Turnover region position:\t ', pos_to_region)
        print ('\n')
         
class CreateOutput():
    """
    
    Args:
        read_name (str) = Name of read
    
    """


    def __init__(self, read_name):
        self.read_name = read_name
    
    @staticmethod
    def prep_output_files(input_file_name):
        """
        Creates all output files names and creates the files themselves including the headers.
        
        Args:
            input_file_name (str): Name of input file
            
        Returns:
            -
        
        """
        data_type = input_file_name[-9:-4]
        CreateOutput.output_file_non_hybrids = 'non_hybrid_reads_{0}.txt'.format(data_type)
        CreateOutput.output_file_zero_reads = 'zero_reads_{0}.txt'.format(data_type)
        CreateOutput.output_file_more_switches = 'hybrid_reads_more_switches_{0}.txt'.format(data_type)
        CreateOutput.output_file_1_switch = 'hybrid_reads_1_switch_{0}.txt'.format(data_type)
        CreateOutput.output_file_overall = 'metadata_{0}.txt'.format(data_type)

        #create output file for non hybrid reads
        with open(CreateOutput.output_file_non_hybrids, 'w') as db_file:
            db_file.write('Read name\tAllele match\n') 
        db_file.close()
    
        #create output file for zero reads
        with open(CreateOutput.output_file_zero_reads, 'w') as db_file:
            db_file.write('Read name\tNote\n') 
        db_file.close()

        #create output file for hyrbid reads with more than 1 switch
        with open(CreateOutput.output_file_more_switches, 'w') as db_file:
            db_file.write('Read name\n') 
        db_file.close()

        #create output file for hyrbid reads with 1 switch
        with open(CreateOutput.output_file_1_switch, 'w') as db_file:
           db_file.write('Read name\tAllele match\tRead1 pos\tRead2 pos\tRead1 mis\tRead2 mis\tRead con mis\tArtefacts\tTurnover region pos\tTurnover sequence\n') 
        db_file.close()      
    
    def non_hybrid_read(self, allele_match, note):
        """
        Adds the non hybrid reads to the output file.
        Zero mismatches between both reads and for same allele
        
        Args:
            allele_match (str): Name of allele (best match)
            note (str): if read consensus has 1 mismatch, then it is stored as a note here 
            
        Returns:
            -    
        """
        print ('Non hybrid read: ', self.read_name)
        
        # Write read data into outfile
        if note == '':
            with open(CreateOutput.output_file_non_hybrids, 'a') as db_file:
                db_file.write(self.read_name + '\t' + allele_match + '\n') 
   
        if note != '':
            with open(CreateOutput.output_file_non_hybrids, 'a') as db_file:
                db_file.write(self.read_name + '\t' + allele_match + '\t' + note + '\n') 
 
    def zero_reads(self, note):
        """
        Adds the zero reads (reads that have 0 mismatches for multiple alleles) to the output file
        
        Args:
            note (str): The alleles with 0 mismatches are noted here
            
        Returns:
            -
        """
        print ('Read with 0 mismatches for multiple alleles:', read_name)

        # add reads that have 0 mismatches for multiple alleles
        with open(CreateOutput.output_file_zero_reads, 'a') as db_file:
            db_file.write(self.read_name + '\t' + note + '\n') 
    
    def hybrid_read_more_switches(self):
        """
        Adds the hybrid reads with more than 1 switch to the output file, these read pairs gave too many
        mismatches which were not indicative
        
        Args:
            -
            
        Returns:
            -    
        """

        print ('Read with more switches: ', read_name)
        
        # add hybrid reads with more switches  
        with open(CreateOutput.output_file_more_switches, 'a') as db_file:
            db_file.write(read_name + '\n') 
   
    def hybrid_read_1_switch(self, allele_name, pos_read1_allele, pos_read2_allele, allele_read1_mismatches, allele_read2_mismatches, allele_consensus_mismatches, turn_over_region, pos_to_region, read_artefacts):
        """
        Adds all hybrid reads with 1 switch data to the output file
        
        Args:
            allele_name (str):
            pos_read1_allele (str):
            pos_read2_allele (srt):
            allele_read1_mismatches (str):
            allele_read2_mismatches (str):
            allele_consensus_mismatches (str):
            turn_over_region (str):
            pos_to_region (str):
            read_artefacts (int):
        
        Return:
            -      
        """

        if turn_over_region == '':
            turn_over_region = '-'

        # add hybrid reads with one switch
        with open(CreateOutput.output_file_1_switch, 'a') as db_file:
            db_file.write(str(self.read_name) + '\t' + str(allele_name) + '\t' + str(pos_read1_allele) + '\t' + str(pos_read2_allele) + '\t' + str(allele_read1_mismatches)  + '\t' + str(allele_read2_mismatches)\
              + '\t' + str(allele_consensus_mismatches) + '\t' + str(read_artefacts) + '\t' + str(pos_to_region) + '\t' + str(turn_over_region) + '\n') 

    @staticmethod
    def metadata(incorrect_aligned_reads, rejected_read_count, non_hybrid_count, zero_count, more_switches_count, one_switch_hybrid_count, total_nr_of_reads):
        """
        Creates metadata output file and adds all read counts
        
        Args:
            incorrect_aligned_reads (int): Number of reads which were incorrect aligned by Clustal Omega
            rejected_read_count (int): Number of reads which did not met the set requirements (read length/ number of N's)
            non_hybrid_count (int): Number of non hybrid reads 
            zero_count (int): Number of zero reads
            more_switches_count (int): Number of reads with more switches
            one_switch_hybrid_count (int): Number of hybrid reads with 1 switch
            total_nr_of_reads (int): The number of reads in total

        Return:
            -
        """
        with open(CreateOutput.output_file_overall, 'w') as db_file:
            db_file.write('Incorrect aligned reads\t' + str(incorrect_aligned_reads) + '\n')       
            db_file.write('Rejected reads\t' + str(rejected_read_count) + '\n')  
            db_file.write('Non hybrid reads\t' + str(non_hybrid_count) + '\n')
            db_file.write('Hybrid reads with more switches\t' + str(more_switches_count) + '\n')
            db_file.write('Hybrid reads with 1 switch\t' + str(one_switch_hybrid_count) + '\n')
            db_file.write('Read with 0 mismatches for multiple alleles\t' + str(zero_count) + '\n')
            db_file.write('Total nr. of reads\t' + str(total_nr_of_reads) + '\n')






if __name__ == "__main__":

    input_file = argv[1]


    CreateOutput.prep_output_files(input_file)
  
    # Class ParseInput (Parse input file and get all allele combinations)

    with open (input_file) as file_object:
            input_file = file_object.read()
            
    all_data, allele_data = ParseInput.collect_all_data(input_file)
    all_combinations_list = ParseInput.get_allele_combinations(allele_data)

    incorrect_aligned_reads = 0
    rejected_read_count = 0
    non_hybrid_count = 0
    zero_count = 0
    more_switches_count = 0
    one_switch_hybrid_count = 0
    read_counter = 1
    for read_info in all_data:
        print ('Number of analyzed reads :', read_counter, '\n')
        read_counter += 1
        print ('~~~~~~~~~~~~~~~~ Read analysis started ~~~~~~~~~~~~~~~~')
        read_name = read_info[0][0]
        print ('Read name :', read_name)
        note = ''
        read1_seq = read_info[0][1]
        read1_qv = read_info[0][2]
        read2_seq = read_info[1][1]
        read2_qv = read_info[1][2]
        read1_aligned_seq = read_info[2][1]
        read2_aligned_seq = read_info[3][1]
        
        allele_data = read_info[4:]
        
        # Perform checks for read 1
        R1_read = Read(read1_seq, read1_qv, read1_aligned_seq, allele_data)
        check_alignment = R1_read.check_alignment()
        if check_alignment == False:
            continue        
        R1_alignment_after_first_check = R1_read.apply_qv()
        R1_alignment_after_second_check = R1_read.check_read_artefacts(R1_alignment_after_first_check)

        # Perform checks for read 2
        R2_read = Read(read2_seq, read2_qv, read2_aligned_seq, allele_data)
        check_alignment = R2_read.check_alignment()
        if check_alignment == False:
            incorrect_aligned_reads += 1
            continue  
        R2_alignment_after_first_check = R2_read.apply_qv()
        R2_alignment_after_second_check = R2_read.check_read_artefacts(R2_alignment_after_first_check)

        # Check if read pair met the requirements
        R1_and_R2 = ReadPair(R1_alignment_after_second_check, R2_alignment_after_second_check, read1_seq, read2_seq)
        approve_reads = R1_and_R2.check_read_pair()
        if approve_reads == True:
            print ('Paired-end read is accepted')
        if approve_reads == False:
            rejected_read_count += 1
            print ('Paired-end read is rejected')
            continue
        R1_mismatch_dict, R1_mismatch_dict_ex = R1_read.get_mismatches(R1_alignment_after_second_check)
        R2_mismatch_dict, R2_mismatch_dict_ex = R2_read.get_mismatches(R2_alignment_after_second_check)

        # Sort alleles, alleles with lowest nr of mismatches first
        mismatch_dict_read1_sorted = sorted(R1_mismatch_dict.items(), key=lambda kv: kv[1])     
        mismatch_dict_read2_sorted = sorted(R2_mismatch_dict.items(), key=lambda kv: kv[1])

        # Count number of alleles with 0 mismatches for read 1
        zero_mismatch_count_read1 = 0
        zero_mismatch_allele_read1 = []
        for allele, mismatches in R1_mismatch_dict.items():
            if mismatches[0] == 0:
                zero_mismatch_count_read1 += 1
                zero_mismatch_allele_read1 += [allele]

        # Count number of alleles with 0 mismatches for read 2
        zero_mismatch_count_read2 = 0
        zero_mismatch_allele_read2 = []
        for allele, mismatches in R2_mismatch_dict.items():
            if mismatches[0] == 0:
                zero_mismatch_count_read2 += 1
                zero_mismatch_allele_read2 += [allele]


        ###########
        ###########  Mismatches per read
        ###########
        ### For non hybrid reads (perfect non hybrid)
        non_hybrid = False
        if zero_mismatch_count_read1 == zero_mismatch_count_read2 == 1:
            # First check if reads are hybrid, if same allele has mismatches for both reads, then it is a non hybrid
            if mismatch_dict_read1_sorted[0][1][0] == 0:
                check_allel = mismatch_dict_read1_sorted[0][0]
                if R2_mismatch_dict[check_allel][0] == 0:
                    allele_match =  mismatch_dict_read1_sorted[0][0]
                    non_hybrid_count += 1
                    read_output = CreateOutput(read_name)
                    read_output.non_hybrid_read(allele_match, note)
                    non_hybrid = True
                    continue

        ### Multiple alleles with 0 mismatches (zero_mismatch_multiple_alleles_count)
        if zero_mismatch_count_read1 != 0 and zero_mismatch_count_read2 != 0:
            if zero_mismatch_count_read1 > 1 or zero_mismatch_count_read2 > 1:
                note = 'Note: Allele(s) {0} has/have 0 mismatches with read 1\tAllele(s) {1} has/have 0 mismatches with read 2'.format(zero_mismatch_allele_read1, zero_mismatch_allele_read2)
                zero_count += 1
                read_output = CreateOutput(read_name)
                read_output.zero_reads(note)
                non_hybrid = True
                continue

        ###########
        ###########  Mismatches read consensus
        ###########
        ### For non hybrid reads with 0 or 1 mismatches in read consensus

        
        # Create read consensus
        alignment_read_consensus = R1_and_R2.create_read_consensus()
        consensus_read = Read.classmethod_for_non_read(alignment_read_consensus, allele_data)
        mismatch_dict_read_con, mismatch_dict_read_con_ex = consensus_read.get_mismatches(alignment_read_consensus)
        mismatch_dict_read_con_sorted = sorted(mismatch_dict_read_con.items(), key=lambda kv: kv[1])
            
        ### if read consensus has 0 or 1 mismatches
        if mismatch_dict_read_con_sorted[0][1][0] == 0 or mismatch_dict_read_con_sorted[0][1][0] == 1:
            allele_match = mismatch_dict_read_con_sorted[0][0]
            if mismatch_dict_read_con_sorted[0][1][0] == 1:
                note = 'Read consensus has 1 mismatch'
            non_hybrid_count += 1
            read_output = CreateOutput(read_name)
            read_output.non_hybrid_read(allele_match, note)
            continue
        
        # print all mismatch info
        R1_read.print_mismatches('First read', R1_mismatch_dict_ex)
        R2_read.print_mismatches('Second read', R2_mismatch_dict_ex)
        consensus_read.print_mismatches('Read consensus', mismatch_dict_read_con_ex)


        ###########
        ###########  Determine number of switches for all allele combinations
        ###########
        more_switches = True
        for allele_combo in all_combinations_list:
            allele1 = allele_combo[0]
            allele2 = allele_combo[1]
         
            per_allele_info = CheckAlleleCombination(alignment_read_consensus, allele_combo, allele_data)
            allele_seq_list = per_allele_info.create_indicator_string()

            check1 = per_allele_info.check_indicative_SNPs()
            if check1 == True:
                check2 = per_allele_info.check_mutual_SNPs()
            else:
                continue
                
            if check2 == True:
                count_indicator_list, number_of_artefacts = per_allele_info.check_alternately_SNPs()
            else:
                continue
                
            if count_indicator_list != None:
                final_indicator_string = per_allele_info.update_indicator_string(count_indicator_list)
            else:
                continue
            # Check if updated indicator string has enough indicative SNPs
            repeat_check1 = per_allele_info.check_indicative_SNPs()
            
            if repeat_check1 == True:
                nr_of_switches, start_turn_pos, end_turn_pos = per_allele_info.get_switches(final_indicator_string)
            
            else:
                continue

            if nr_of_switches == 1:
                more_switches = False
                per_allele_info.print_1_switch_alleles()
                read1_pos_dict = R1_read.get_relative_position()
                read2_pos_dict = R2_read.get_relative_position()
                
                final_to_region = GetOneSwitchData(allele1, allele2)

                pos_read1_allele1, pos_read2_allele1, pos_read1_allele2, pos_read2_allele2 = final_to_region.get_read_position(read1_pos_dict, read2_pos_dict)
                turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = final_to_region.prep_for_turnover_position(start_turn_pos, end_turn_pos, allele_seq_list, allele_data)

                TO1_seq = Read.classmethod_for_non_read(turn_over_region1_for_pos, seq_list_allele1)
                TO2_seq = Read.classmethod_for_non_read(turn_over_region2_for_pos, seq_list_allele2)

                TO_allele1_dict = TO1_seq.get_relative_position()
                TO_allele2_dict = TO2_seq.get_relative_position()

                pos_to_region1, pos_to_region2, turn_over_region1, turn_over_region2 = final_to_region.get_TO_position(TO_allele1_dict, TO_allele2_dict, turn_over_region1_for_pos, turn_over_region2_for_pos)
                GetOneSwitchData.print_TO_output(allele1, pos_read1_allele1, pos_read2_allele1, turn_over_region1, pos_to_region1)
                GetOneSwitchData.print_TO_output(allele2, pos_read1_allele2, pos_read2_allele2, turn_over_region2, pos_to_region2)
                
                allele1_read1_mismatches = str(R1_mismatch_dict[allele1][0])
                allele1_read2_mismatches = str(R2_mismatch_dict[allele1][0])
                allele2_read1_mismatches = str(R1_mismatch_dict[allele2][0])
                allele2_read2_mismatches = str(R2_mismatch_dict[allele2][0])
                allele1_consensus_mismatches = str(mismatch_dict_read_con[allele1][0])
                allele2_consensus_mismatches = str(mismatch_dict_read_con[allele2][0])

                read_output = CreateOutput(read_name)
                read_output.hybrid_read_1_switch(allele1, pos_read1_allele1, pos_read2_allele1, allele1_read1_mismatches, allele1_read2_mismatches, allele1_consensus_mismatches, turn_over_region1, pos_to_region1, number_of_artefacts)
                read_output.hybrid_read_1_switch(allele2, pos_read1_allele2, pos_read2_allele2, allele2_read1_mismatches, allele2_read2_mismatches, allele2_consensus_mismatches, turn_over_region2, pos_to_region2, number_of_artefacts)

        if more_switches == False:
            print ('Hybrid read with 1 switch: ', read_name)
            one_switch_hybrid_count += 1
        
        # For hybrid reads with more switches (too many mismatches)        
        if more_switches == True:
            more_switches_count += 1
            read_output = CreateOutput(read_name)
            read_output.hybrid_read_more_switches()
       
            
    total_nr_of_reads = read_counter-1 + incorrect_aligned_reads + rejected_read_count
    CreateOutput.metadata(incorrect_aligned_reads, rejected_read_count, non_hybrid_count, zero_count, more_switches_count, one_switch_hybrid_count, total_nr_of_reads)
   