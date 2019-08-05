"""
10-07-'19

This script categorizes reads in non hybrid reads, zero reads, hybrid reads with 1 switch and hybrid reads with more switches.
The reads are categorized based on mismatches (SNPs) and the number of switches. 
Also metadata is generated and the file with 1 switch data contains the most extended information.

Input required: output file with alignments from AlignReads.py

"""

from sys import argv
import os

class ParseInput:
    """
    This class prepares the input data for further processing.

    Args:
        -
    """
    
    @staticmethod
    def collect_all_data(input_file):
        """
        Separetes the input data (aligned read pairs and the best matching alleles) per read pair.
        
        Args:
            input_file (str): aligned read pairs and the best matching alleles (all lines from output file AlignRead.py).
        Returns:
            all_data (list): list of lists, each list contains the read information (read name and quality values), allele
            names and all alignments for the reads and best matches (HLA-A, B and C). 
            allele_names (list): contains all allele names (max. 6)
        """

        all_data = []
        input_file = input_file.split('$$$')
        for per_read_pair_data in input_file[1:-2]:   # 1:-1
            temp_collect_list = []
            read_data = per_read_pair_data.split('\n')
            for line in read_data:
                line = line.split('\t')
                if len(line) > 1:
                    temp_collect_list += [line]
            all_data += [temp_collect_list]
        allele_names_list = all_data[0][4:]

        allele_names = []
        for allele_name_line in allele_names_list:
            allele_names += [allele_name_line[0]]

        return all_data, allele_names

    @staticmethod
    def get_allele_combinations(allele_names):
        """
        Generates all possible allele combinations (allele names) within and between loci.
        
        Args:
            allele_names (list): contains all allele names (max. 6)
        Returns:
            all_allele_combinations (list): contains all possible allele name combinations
        """

        if len(allele_names) < 4:
            raise ValueError ('The sample contains less than 4 alleles!')

        # for samples with 4 different alleles
        if len(allele_names) == 4:
            first_nr_list = [0,0,0,1,1,2]
            second_nr_list = [1,2,3,2,3,3] 

        # for samples with 5 different alleles
        if len(allele_names) == 5:
            first_nr_list = [0,0,0,0,1,1,1,2,2,3]
            second_nr_list = [1,2,3,4,2,3,4,3,4,4] 
        
        # for samples with 6 different alleles
        if len(allele_names) == 6:
            first_nr_list = [0,0,0,0,0,1,1,1,1,2,2,2,3,3,4]
            second_nr_list = [1,2,3,4,5,2,3,4,5,3,4,5,4,5,5]

        all_allele_combinations = []
        for i in range(len(first_nr_list)):
            all_allele_combinations += [[allele_names[int(first_nr_list[i])], allele_names[int(second_nr_list[i])]]]
        
        return (all_allele_combinations)
        
class Read:
    """
    This class processes single read data but some functions can also process single sequence data other than reads.
    Several checks are included (quality and artefacts), also mismatch and position information can be generated.

    Args:
        read_seq (str): sequence, not in alignment, from read, read consensus or turnover region
        read_aligned_seq (str): sequence in alignment, from read, read consensus or turnover region
        allele_data (list): list of lists with all allele names and aligned sequences
    """
    def __init__(self, read_seq, read_aligned_seq, allele_data):
        self.read_seq = read_seq
        self.read_length = len(read_seq)
        self.read_aligned_seq = read_aligned_seq
        self.allele_data = allele_data

    def check_alignment(self): 
        """
        Checks if the read pair is correctly aligned. If clustal omega was not able to do so, then the alignment of the
        read prior to it was copied. 

        Args:
            -
        Returns:
            correct_alignment (bool): True if sequences are correctly aligned, False if sequences are incorrecly aligned
        """
        # Check if alignment is correct
        
        correct_alignment = True
        check_alignment = self.read_aligned_seq.replace('-', '')

        if self.read_seq != check_alignment:
            correct_alignment = False

        if self.read_seq == '':
            raise ValueError ('Read sequence is empty!')

        return correct_alignment
    
    def apply_qv(self, read_qv):
        """
        Checks nucleotide quality values, if lower than a given value (minimum_q_score), then the nucleotide is replaced by a 'N'
        
        Args:
            read_qv (str): read quality values
        Returns:
            read_checked_aligned_seq (str): the updated aligned read sequence
        """

        # Adjust here the minimum quality score 
        minimum_q_score = 18

        quality_dict = {'!': 0, '"': 1, '#':2, '$':3, '%':4, '&':5, "'":6, '(':7, ')':8,\
        '*':9, '+':10, ',':11,'-':12, '.':13, '/':14, '0':15, '1':16, '2':17, '3':18,\
        '4':19, '5':20, '6':21, '7':22, '8':23, '9':24, ':':25, ';':26, '<':27, '=':28,\
        '>':29, '?':30, '@':31, 'A':32, 'B':33, 'C':34, 'D':35, 'E':36, 'F':37, 'G':38, 'H':39, 'I':40}
    
        read_checked_seq = '' 
        for i, qual in enumerate(read_qv):
            q_score = quality_dict[qual]
            nucleotide = self.read_seq[i]
            if q_score < minimum_q_score:
                read_checked_seq += 'N'
            if q_score >= minimum_q_score:
                read_checked_seq += nucleotide

        # count the number of '-' in front of aligned read (left)
        count_read_start = self.read_aligned_seq.lstrip('-')
        count_read_start = len(self.read_aligned_seq) - len(count_read_start)

        # remove '-' right from aligned read and count the number of '-' after the aligned read (right)
        read_aligned_seq_wo_r = self.read_aligned_seq.rstrip('-')
        count_read_end = len(self.read_aligned_seq) - len(read_aligned_seq_wo_r)

        read_aligned_seq_checked = ''
        read_start_seen = False   
        gap_count = 0
        for i, char in enumerate(read_aligned_seq_wo_r):
            if char == '-':
                read_aligned_seq_checked += char
            if char != '-' and read_start_seen == False:
                read_start_seen = True
            if read_start_seen == True:
                if char == '-':
                    gap_count += 1
                if char != '-':
                    read_aligned_seq_checked += read_checked_seq[i-count_read_start-gap_count]

        # add '-' at the right
        read_aligned_seq_checked = read_aligned_seq_checked + '-' * count_read_end

        return read_aligned_seq_checked
    

    def check_read_artefacts(self, read_aligned_seq_checked):
        """" 
        Check if read has artefect type 1. If artefact is found then the read nucleotide is replaced by a 'N'
        Artefact type  1 definition: if all alleles have have a mismatch at the same position, we assume that
        the read nucleotide is incorrect (e.g. PCR artefact) not the allele nucleotide.
        
        Args:
            read_aligned_seq_checked (str): the updated aligned read sequence
        Returns:
            read_aligned_seq_fully_checked (str): the updated aligned read sequence after both checks
        """

        # get the start of the read position (based on the read in its alignment) 
        start_absolute_read_position = []
        start_absolute_read_nucleotide = []
        nr_of_switches = 0
        nuc = read_aligned_seq_checked[0]
        for i, char in enumerate(read_aligned_seq_checked):
            if nuc == '-' and nuc != char:
                nr_of_switches += 1
                start_absolute_read_position += [i]
                start_absolute_read_nucleotide += [char]
            nuc = char
    
        # get all read positions (based on the read in its alignment)
        absolute_read_position = []
        absolute_read_nucleotide = []
        for i, char in enumerate(read_aligned_seq_checked):
            if i >= start_absolute_read_position[0] and i <= start_absolute_read_position[-1]:
                absolute_read_position += [i]
                absolute_read_nucleotide += [char]
            if i >= start_absolute_read_position[-1]+1 and char != '-':
                absolute_read_position += [i]
                absolute_read_nucleotide += [char]
       
        mismatch_track = self.__create_mismatch_track(read_aligned_seq_checked, absolute_read_position, absolute_read_nucleotide)

        # if a nucleotide in the read is considered as artefact (5 in mismatch track) then it is replaced by 'N'
        read_aligned_seq_fully_checked = ''
        for i, nuc in enumerate(read_aligned_seq_checked):
            try:
                mismatch_char = int(mismatch_track[i])
            except:
                mismatch_char = mismatch_track[i]
            number_of_alleles = len(self.allele_data)
            if number_of_alleles == 5:
                if mismatch_char == 5:
                    read_aligned_seq_fully_checked += 'N'
                if mismatch_char != 5:
                    read_aligned_seq_fully_checked += nuc
            elif number_of_alleles == 6:
                if mismatch_char == 6:
                    read_aligned_seq_fully_checked += 'N'
                if mismatch_char != 6:
                    read_aligned_seq_fully_checked += nuc
            else:
                raise ValueError ('The number of alleles is incorrect!')

        return read_aligned_seq_fully_checked

    def __create_mismatch_track(self, read_aligned_seq_checked, absolute_read_position, absolute_read_nucleotide):
        """
        Checks for each read position if alleles have mismatches (substitution). If a mismatch is found, value 1 is added.
        If the input data contains 6 alleles, then this number can be max. 6 for each read position.
        
        Args:
            read_aligned_seq_checked (str): the updated aligned read sequence
            absolute_read_position (list): contains all positions (int) of the read (based on its alignment)
            absolute_read_nucleotide (list): contains all read nucleotides which corresponds to the absolute read positions
        Returns:
            mismatch_track (str): contains the number of mismatches for all alleles per position
        """

        # create a mismatch track string, for each mismatch in the allele '1' is added, up to 5  (where all alleles have mismatches)
        mismatch_track = '-' * len(read_aligned_seq_checked)

        for allele, seq in self.allele_data:
            seq_string = seq
            for i, chari in enumerate(seq_string):
                if i in absolute_read_position:
                    read_nuc = absolute_read_nucleotide[i-min(absolute_read_position)] 
                    if chari != read_nuc and read_nuc != '-' and read_nuc != 'N' and chari != '-' and  mismatch_track[i] == '5':
                        left_mismatch_track = mismatch_track[:i]
                        right_mismatch_track = mismatch_track[i+1:]
                        mismatch_track = left_mismatch_track + '6' + right_mismatch_track
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
        Gets all allele mismatches (substitutions, insertions and deletions) based on the read.
        
        Args:
            read_aligned_seq_fully_checked (str): the updated aligned read sequence after both checks
        Returns:
           mismatch_dict (dict): contains allele names and number of total mismatches
           extended_mismatch_dict (dict): contains allele names and number of substitutions, insertions and deletions
        """

        # get the start of the read position (based on the read in its alignment) 
        start_absolute_read_position = []
        nr_of_switches = 0
        nuc = read_aligned_fully_checked[0]
        for i, char in enumerate(read_aligned_fully_checked):
            if nuc == '-' and nuc != char:
                nr_of_switches += 1
                start_absolute_read_position += [i]
            nuc = char

        # If read is aligned in front of the allele, the mismatches are ignored (they will be extracted from deletions). 
        deletion_correction_dict = self.__check_read_start(start_absolute_read_position, read_aligned_fully_checked)

        absolute_read_position = []
        absolute_read_nucleotide = []
        for i, char in enumerate(read_aligned_fully_checked):
            if i >= start_absolute_read_position[0] and i <= start_absolute_read_position[-1]:
                absolute_read_position += [i]
                absolute_read_nucleotide += [char]
            if i >= start_absolute_read_position[-1]+1 and char != '-':
                absolute_read_position += [i]
                absolute_read_nucleotide += [char]

        extended_mismatch_dict = {}
        mismatch_dict = {}
        for allele, seq_string in self.allele_data:
            substitutions = 0
            insertions = 0
            deletions = int(deletion_correction_dict[allele])
            mismatches = 0
            for i, chari in enumerate(seq_string):
                if i in absolute_read_position:
                    read_nuc = absolute_read_nucleotide[i-min(absolute_read_position)]
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
  


    def __check_read_start(self, start_absolute_read_position, read_aligned_fully_checked):
        """
        Checks if aligned read starts in front of allele (does not occur often). If so, then the number of nucleotides
        which are aligned in front of the allele are counted and stored in a dictionary.  
       
        Args:
            start_absolute_read_position (list): contains absolute start position of the read and the start positions after
            insertions (if the read has any).
            read_aligned_seq_fully_checked (str): the updated aligned read sequence after both checks
        Returns:
            deletion_correction_dict (dict): contains allele names and number of nucleotides which are aligned in front 
            of the allele
        """

        deletion_correction_dict = {} # if reads starts in front of allele
        for allele, seq_string in self.allele_data: 
            start_allele_post = 0
            if seq_string.startswith('-'):
                start_absolute_allele_position = []
                nuc = seq_string[0]
                nr_of_switches = 0
                for i, char in enumerate(seq_string):
                    if nuc == '-' and nuc != char:
                        nr_of_switches += 1
                        start_absolute_allele_position += [i]
                    nuc = char
                start_allele_post = start_absolute_allele_position[0]
                if start_allele_post > start_absolute_read_position[0]:
                    deletion_correction = start_absolute_read_position[0] - start_allele_post
                    deletion_correction_dict[allele] = deletion_correction
                else: 
                    deletion_correction = 0
                    deletion_correction_dict[allele] = deletion_correction
            else: 
                deletion_correction = 0
                deletion_correction_dict[allele] = deletion_correction
            count_read_insertion_for_deletion_correction = read_aligned_fully_checked[start_absolute_read_position[0]:start_allele_post]
            read_inserts = count_read_insertion_for_deletion_correction.count('-')
            deletion_correction_dict[allele] += read_inserts

        return deletion_correction_dict
   
    @classmethod
    def classmethod_for_non_read(cls, aligned_sequence, allele_data):
        """
        Classmethod that generates data for the read consensus and turnover region (both in alignment) as input for the constructor. 
        The sequence (read_seq) is created  without its alignment. 
        
        Args:
            aligned_sequence (str): read consensus or turnover region sequence, in its alignment
            allele_data (list): list of lists with all allele names and aligned sequences
        Returns:
            read_seq (str): sequence, not in alignment, read consensus or turnover region
            read_aligned_seq (str): sequence in alignment, read consensus or turnover region
            allele_data (list): list of lists with all allele names and aligned sequences
        """
        read_seq = aligned_sequence.lstrip('-').rstrip('-')

        return cls(read_seq, aligned_sequence, allele_data)
    
    def print_mismatches(self, read_type, extended_mismatch_dict):
        """
        Print type of read, read length and the number mismatches
        
        Args:
            read_type (str): discribes read type; First read, Second read or Read consensus
            extended_mismatch_dict (dict): contains the number of SNP substitutions, insertions, deletions and total number of mismatches
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
        Determines all positions of a given sequence (read, read consensus or turnover region) per nucleotide relative to 
        the allele (the positions of the allele nucleotides that are covered by the given sequence). All alleles are included. If 
        the sequence starts in front of the allele (does not occur often), then those positions are ignored, they do not exist. 
        Same goes for gaps in alleles (sequence has nucleotide and allele does not).
        
        Args:
            -            
        Returns:
            read_pos_dict (dict): contains allele names and all positions of the given sequence relative to the alleles
        """
        read_pos_dict = {}
        for allele, allele_seq in self.allele_data:
            read_position = []
            # remove '-' left and right from allele seq and give read the same length
            allele_seq_wo_left = allele_seq.lstrip('-')
            left_difference = len(allele_seq) - len(allele_seq_wo_left)
            allele_seq_wo_right = allele_seq.rstrip('-')

            read_seq_wo_rl = self.read_aligned_seq[left_difference:len(allele_seq_wo_right)]
            allele_seq_wo_rl = allele_seq.strip('-')

            # Remove '-' right from the read, example sequence read: '--------------CCCCC'
            sequence_read = read_seq_wo_rl.rstrip('-')

            read_start_seen = False
            deletion_count_allele_to_read_start = 0
            deletion_count_activated = False
            #TODO: now the allele can have max. 2 gaps in front of read start

            for i, char in enumerate(sequence_read):
                if deletion_count_activated == True and char == '-' and read_start_seen ==  False:
                    if allele_seq_wo_rl[i] == '-':  # if allele has gap in front of read start, for second gap
                        deletion_count_allele_to_read_start += 1
                if char == '-' and read_start_seen ==  False and deletion_count_activated == False:
                    if allele_seq_wo_rl[i] == '-':  # if allele has gap in front of read start, for first gap
                        deletion_count_allele_to_read_start += 1
                        deletion_count_activated = True
                if char != '-' and deletion_count_activated == False:
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
            read_pos_dict[allele] = read_position

        # get position if turnover region has a length of 0 and the allele has '-' as nucleotide
        allele_name =  self.allele_data[0][0]

        if 'K' in self.read_aligned_seq or 'Z' in self.read_aligned_seq and read_pos_dict[allele_name] == []:
            read_pos_dict = self.__get_special_case_pos(read_pos_dict, allele_name)

        return read_pos_dict


    def __get_special_case_pos(self, read_pos_dict, allele_name):
        """
        If turnover region has a length of 0 (indicated by a 'K') or 1  (indicated by a 'Z') and the allele 
        has '-' as nucleotide. This function gets the correct position.
        
        Args:
            read_pos_dict (dict): contains allele names and all positions of the given sequence relative to the alleles
            allele_name (str): name of allele
        Returns:
            read_pos_dict (dict): an updated version of original read_pos_dict
        """

        read_start_seen = False
        pos = 0
        read_position = []
        for allele, allele_seq in  self.allele_data:
            
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
            read_pos_dict[allele_name] = read_position

        return (read_pos_dict)

        
class ReadPair():
    """
     
    Args:
        read1_aligned_checked (str): sequence read 1, in alignement
        read2_aligned_checked (str): sequence read 2, in alignement
        read1_seq (str): sequence read 1
        read2_seq (str): sequence read 2
        min_read_length (int): the minimum read length allowed
        N_quantity (int): the maximum number of N's allowed per read
    """
    
    def __init__(self, read1_aligned_checked, read2_aligned_checked, read1_seq, read2_seq):
        self.read1_aligned_checked = read1_aligned_checked
        self.read2_aligned_checked = read2_aligned_checked
        self.read1_seq = read1_seq
        self.read2_seq = read2_seq

    def check_read_pair(self, min_read_length, N_quantity):
        """
        Checks the minimum read length of the original reads and the number N's per aligned and checked read. Returns a 
        bool that indicates or the reads of a read pair met the requirements or not.
        
        Args:
            -
            
        Returns:
            approve_reads (bool): True if reads are longer than set value and the number of N's is lower than set value.
            False if one or both reads do not met the set values.
        """

        approve_reads = True

        read1_length = len(self.read1_seq)
        read2_length = len(self.read2_seq)
        read1_N_count = self.read1_aligned_checked.count('N')
        read2_N_count = self.read2_aligned_checked.count('N')

        # check read length
        if read1_length < min_read_length or read2_length < min_read_length:    
            approve_reads = False

        # check nr of Ns
        if read1_N_count > N_quantity or read2_N_count > N_quantity:    
            approve_reads = False

        return approve_reads

   
    def create_read_consensus(self):
        """
        This function combines the two paired-end reads into a 'read consensus' if the reads do not 
        overlap (and thus have a gap) the empty nucleotide between the reads are replaced by '*'. 
        This to distinguish between real read deletions and the empty space between the reads. 
        For overlapping reads: if the reads have different nucleotides at the same position then it 
        is replaced by a 'N'. If 1 read has a 'N' and the other one a nucleotide, the nucleotide is used.
        
        Args:
            -
        Returns:
            read_consensus (str): contains read pair sequences combined, '*' indicates the gap between the reads
        """
        read1_seq = self.read1_aligned_checked
        read2_seq = self.read2_aligned_checked
        
        if len(read1_seq) != len(read2_seq):
            raise ValueError ('Aligned reads have a different length!')

        
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
    This class processes each given allele combination. First an indicator string is created with indicative
    mismatches, then this string is checked and updated. Based on this string, the number of swicthes are
    determined. Allele combinations that result 1 switch (perfect hybrid reads) are the main focus.

    Args:
        read_consensus (str): contains read pair sequences combined, '*' indicates the gap between the reads
        allele_combo (list): allele names of given combination
        allele_data (list): list of lists with all allele names and aligned sequences
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
        Here, the indicator string is created based on the alleles of the given allele combination. The read 
        consensus is used as reference, mismatches for the first allele are indicated by a 'X' and mismatches for
        the second allele are indicated by a 'Y'. The order does not matter. The indicator strings are first 
        created separately and then they are combined. The allele_seq_list is needed in a later stage.
        
        Args:
            -
        Returns:
            allele_seq_list (list): contains allele sequences in alignment for given allele combination
        """

        # Create a dict with only the given allele combination. The mismatch indicator string contains '-' for matches and 'X' or 'Y' for a mismatches.
        turn_seq_list = []
        allele_seq_list = []
        mismatch_char = 'X'
        start_consensus = self.read_consensus.rstrip('-').count('-') # start read pos (absolute)
        end_consensus = len(self.read_consensus.rstrip('-'))  # end read pos (absolute)

        for allele, seq_string in self.allele_data:
            mismatches = 0
            seq_string_temp = ''
            turn_over_seq_temp = ''

            if allele in self.allele_combo:
                allele_seq_list += [seq_string]
            if allele in self.allele_combo:
                for i, chari in enumerate(seq_string):
                    if chari != self.read_consensus[i]:
                        if self.read_consensus[i] != '-' and chari != '-' and self.read_consensus[i] != '*' and self.read_consensus[i] != 'N':
                            mismatches += 1    
                            turn_over_seq_temp += mismatch_char
                        if self.read_consensus[i] == '-' and chari != '-':  
                            if i < start_consensus or i >= end_consensus: # nucleotide in front of and after reads
                                turn_over_seq_temp += '-'
                            if i >= start_consensus and i < end_consensus:  # allele insertion
                                turn_over_seq_temp += mismatch_char
                        if self.read_consensus[i] != '-' and chari == '-' and self.read_consensus[i] != '*' and self.read_consensus[i] != 'N':   # allele deletion
                            turn_over_seq_temp += mismatch_char
                        if self.read_consensus[i] == '*' or self.read_consensus[i] == 'N':
                            turn_over_seq_temp += '-'
                    if chari == self.read_consensus[i]: 
                        seq_string_temp += chari
                        turn_over_seq_temp += '-'
                mismatch_char = 'Y'
                turn_seq_list += [[allele,turn_over_seq_temp]]

        # Extract the mismatch indicator sequences for allele match 1 and 2
        seq_allele_1 = turn_seq_list[0][1]
        seq_allele_2 = turn_seq_list[1][1]
        if allele_seq_list[0] == allele_seq_list[1]:
            raise ValueError ('Aligned allele sequences (from allele combo) are identical!')

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
        Checks if alleles have enough indicative mismatches, based on the mismatch indicator string. 
        At least 2 mismatches per allele are required.
        
        Args:
            -
        Returns:
            accept_combo (bool): True if the mismatch indicator string contains at least 2 'X' and at least 2 
            'Y' mismatches. False if not, the allele combination has not enough indicative mismatches.
        """

        accept_combo = True
        
        if self.indicator_string.count('X') < 2 or self.indicator_string.count('Y') < 2:
            accept_combo = False

        return accept_combo

    def check_mutual_SNPs(self):
        """
        Checks if alleles do not have too many mutual mismatches (max. 2). If both alleles have a mismatch at
        the same position then we assume that the read has an artefact. But this is only allowed twice.
                
        Args:
            -
        Returns:
            accept_combo (bool): True if alleles have maximum 2 mutual mismatches, False if alleles
            contain more mutual mismatches.
        """
        accept_combo = True

        self.number_of_artefacts += self.indicator_string.count('M') 
        if self.indicator_string.count('M') > 2:
            accept_combo = False

        return accept_combo
        
    def check_alternately_SNPs(self):
        """
        Checks if alleles do not have too many alternately mismatches (max. 2). An alternately mismatch is a single
        mismatch of one allele between 2 mismatches of the other allele. 'XYX' or 'YXY' in the indicator string.
                
        Args:
            -
        Returns:
            count_indicator_list (list): contains ascending values (int), starting from 1, each new start
            indicates a switch for 'X' to 'Y' or vice versa. Two 1's next to each other indicate an
            alternately mismatch.  
            number_of_artefacts (int): total number of artefacts (caused by mutual or alternately SNPs)
        """

        # replace 'M' for '-'  if the alleles have a mutual mismatch, it is ignored and it is regarded as a read artefact
        mismatch_indicator_string = self.indicator_string.replace('M', '-')

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
        
        if count_indicator_list[1] == 1 and count_indicator_list[2] == 1 :  # if the indiactor string starts with a single X or Y
            pcr_artefact -= 1

        for i, number in enumerate(count_indicator_list):
            if i < len(count_indicator_list)-1:
                next_number = count_indicator_list[i+1]
                if number == 1 and next_number == 1:
                    pcr_artefact += 1

        self.number_of_artefacts += pcr_artefact
        number_of_artefacts = self.number_of_artefacts

        # check for allele artefact, XYX or YXY (max. 2) 
        if pcr_artefact > 2:
            count_indicator_list = None

        return count_indicator_list, number_of_artefacts
        
    def update_indicator_string(self, count_indicator_list):
        """
        Here the mutual and alternately mismatches are removed (if there are less then 5) and the mismatch indicator string
        is updated.
        
        Args:
            count_indicator_list (list): contains ascending values (int), starting from 1, each new start
            indicates a switch for 'X' to 'Y' or vice versa. Two 1's next to each other indicate an
            alternately mismatch. 
        Returns:
            final_indicator_string (str): an updated version of the indicator string without mutual and alternately
            mismatches.
        """

        mismatch_indicator_string = self.indicator_string.replace('M', '-')

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
        Determines the number of switches based on the indicator string. A switch is a shift from an 'X' to an 'Y' or
        vice versa. The ideal situation results in 1 switch.
        
        Args:
            final_indicator_string (str): an updated version of the indicator string without mutual and alternately
            mismatches.
        Returns:
            nr_of_switches (int): The number of switches; from 'X' to 'Y' or vice versa
            start_turn_pos (int): absolute start position, in alignment, first nucleotide in turnover region
            end_turn_pos (int): absolute end position, in alignment, last nucleotide in turnover region
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

        if nr_of_switches != 1:
            start_turn_pos = None
            end_turn_pos = None

        return (nr_of_switches, start_turn_pos, end_turn_pos)
    
    def print_1_switch_alleles(self):
        """
        Prints allele combination names if they resulted in an 1 switch hybrid read, also the
        number of artefacts is printed (max. 4).
        
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
    This class processes reads and allele combinations which resulted 1 switch. The read positions relative
    to both alleles are parsed. Also the turnover region sequence and its relative position are determined.
    
    Args:
        allele1 (str): name allele 1
        allele2 (str): name allele 2
    """    

    def __init__(self, allele1, allele2):
        self.allele1 = allele1
        self.allele2 = allele2
         
    def get_read_position(self, read1_pos_dict, read2_pos_dict):
        """
        Parses and extracts the first and last value from the relative read 1 and read 2 positions for both alleles.
        
        Args:
            read1_pos_dict (dict): contains allele names and all positions of read 1 relative to the alleles
            read2_pos_dict (dict): contains allele names and all positions of read 2 relative to the alleles
            
        Returns:
            position_read1_allele1 (str): start and end position of read 1 relative to allele 1
            position_read2_allele1 (str): start and end position of read 2 relative to allele 1
            position_read1_allele2 (str): start and end position of read 1 relative to allele 2
            position_read2_allele2 (str): start and end position of read 2 relative to allele 2     
        """

        # Get all start and stop positions from best allele matches
        positions_list_allele1_read1 = read1_pos_dict[self.allele1]
        start_pos_allele1_read1 = positions_list_allele1_read1[0]
        end_pos_allele1_read1 = positions_list_allele1_read1[-1]

        positions_list_allele2_read1  = read1_pos_dict[self.allele2]
        start_pos_allele2_read1 = positions_list_allele2_read1[0]
        end_pos_allele2_read1 = positions_list_allele2_read1[-1] 

        positions_list_allele1_read2 = read2_pos_dict[self.allele1]
        start_pos_allele1_read2 = positions_list_allele1_read2[0]
        end_pos_allele1_read2 = positions_list_allele1_read2[-1] 

        positions_list_allele2_read2 = read2_pos_dict[self.allele2]
        start_pos_allele2_read2 = positions_list_allele2_read2[0]
        end_pos_allele2_read2 = positions_list_allele2_read2[-1] 

        position_read1_allele1 = str(start_pos_allele1_read1) + '-' +  str(end_pos_allele1_read1)
        position_read2_allele1 = str(start_pos_allele1_read2) + '-' +  str(end_pos_allele1_read2)
        position_read1_allele2 = str(start_pos_allele2_read1) + '-' +  str(end_pos_allele2_read1)
        position_read2_allele2 = str(start_pos_allele2_read2) + '-' +  str(end_pos_allele2_read2)
        
        return (position_read1_allele1, position_read2_allele1, position_read1_allele2, position_read2_allele2)

    def prep_for_turnover_position(self, start_pos, end_pos, allele_seq_list, allele_data): 
        """
        Prepares turnover region data. This function gets the turnover region in alignment for both alleles
        and also the allele sequences in separeted lists, which can be used in order to determine the postions
        of the nucleotides of the turnover region sequence later on.

        Args:
            start_pos (int): absolute start position, in alignment, first nucleotide in turnover region
            end_pos (int): absolute end position, in alignment, last nucleotide in turnover region
            allele_seq_list (list): contains allele sequences in alignment for given allele combination
            allele_data (list): list of lists with all allele names and aligned sequences
        Returns:
            turn_over_region1_for_pos(str): turnover sequence region in alignment for first allele
            seq_list_allele1 (list): allele name and aligned sequence for first allele in allele combination
            turn_over_region2_for_pos(str): turnover sequence region in alignment for second allele
            seq_list_allele2 (list): allele name and aligned sequence for second allele in allele combination
        """

        # First get the absolute turnover sequence positions for both alleles 
        turn_over_region1_for_pos = self.__get_pos_TO_region(allele_seq_list[0], start_pos, end_pos)
        turn_over_region2_for_pos = self.__get_pos_TO_region(allele_seq_list[1], start_pos, end_pos)
 
        if start_pos == end_pos +1:
            print ('Turnover sequence contains 0 nucleotides')

        # Extract allele sequences
        seq_dict_allele1 = {}
        seq_dict_allele2 = {}
        for allele, seq in allele_data:
            if self.allele1 == allele:
                seq_dict_allele1[allele] = seq
            if self.allele2 == allele:
                seq_dict_allele2[allele] = seq
      
        seq_list_allele1 = sorted(seq_dict_allele1.items())
        seq_list_allele2 = sorted(seq_dict_allele2.items())

        return turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2

    @staticmethod
    def __get_pos_TO_region(aligned_allele, start_pos, end_pos):
        """
        Generates the turnover region in alignment for the given allele. This need to be done for
        both alleles independently as the alleles can differ from each other at those positions.
        If the turnover region consist of 0 nucleotides, then the '-' of the second (imaginary)
        nucleotide in the switch (e.g. XY, position Y) is replaces by a 'K', which is the indicator.      
      
        Args:
            aligned_allele (str): the allele sequence in alignment
            start_pos (int): absolute start position, in alignment, first nucleotide in turnover region
            end_pos (int): absolute end position, in alignment, last nucleotide in turnover region
        Returns:    
            turn_over_region_for_pos (str): the turnover region sequence in alignment from given allele
        """

        turn_over_region_for_pos = ''
        check_for_empty_char_string = ''
        empty_string =  True
        for i, char in enumerate(aligned_allele):
            if start_pos == end_pos:    # for TO with length 1
                empty_string = False
                if i != start_pos:
                    turn_over_region_for_pos += '-'
                if i == start_pos:
                    if char == '-':     # if allele has '-' as nucleotide
                        char = 'Z'
                    turn_over_region_for_pos += char
            if start_pos == end_pos + 1:    # for TO with length 0
                empty_string = False
                char = 'K'     
                if i != start_pos:
                    turn_over_region_for_pos += '-'
                if i == start_pos:
                    turn_over_region_for_pos += char
            if start_pos != end_pos + 1 and start_pos != end_pos:  # for TO > length 1
                if i < start_pos:
                    turn_over_region_for_pos += '-'
                    check_for_empty_char_string += '-'
                if i > end_pos:
                    turn_over_region_for_pos += '-'
                    check_for_empty_char_string += '-'
                if i >= start_pos and i <= end_pos:
                    if char != '-':
                        turn_over_region_for_pos += char
                        empty_string = False
                    if char == '-':
                        turn_over_region_for_pos += char
                        check_for_empty_char_string += 'Z'
       
        if empty_string == True:
            turn_over_region_for_pos = check_for_empty_char_string

        return (turn_over_region_for_pos)
       
    def get_TO_position(self, TO_allele1_dict, TO_allele2_dict, turn_over_region1_for_pos, turn_over_region2_for_pos):
        """
        Gets start and end positions for turnover region for both alleles. If the turnover region consist
        of 0 or 1 nucleotide, then the start and end position are equal, and just the start position is selected. 
        
        Args:
            TO_allele1_dict (dict): contains allele name and turnover positions relative to allele 1
            TO_allele2_dict (dict): contains allele name and turnover positions relative to allele 2
            turn_over_region1_for_pos (str): turnover region sequence in aligenment for allele 1
            turn_over_region2_for_pos (str): turnover region sequence in aligenment for allele 2
        Returns:
            position_to_region1 (str): start and end position of turnover region for allele 1
            position_to_region2 (str): start and end position of turnover region for allele 2
            turn_over_region1 (str): turnover region sequence for allele 1
            turn_over_region2 (str): turnover region sequence for allele 2
        """



        # Get all turnover start and end positions for allele 1
        positions_list_allele1_TO_region = TO_allele1_dict[self.allele1]
        start_pos_allele1_TO_region = positions_list_allele1_TO_region[0]
        end_pos_allele1_TO_region = positions_list_allele1_TO_region[-1]
        
        # Get all turnover start and end positions for allele 2
        positions_list_allele2_TO_region  = TO_allele2_dict[self.allele2]
        start_pos_allele2_TO_region = positions_list_allele2_TO_region[0]
        end_pos_allele2_TO_region = positions_list_allele2_TO_region[-1] 

        # delete deletions (just keep the sequence), and K, since it is a non existing TO sequence (lenght 0)
        turn_over_region1 = turn_over_region1_for_pos.replace('K','-').replace('-','')
        turn_over_region2 = turn_over_region2_for_pos.replace('K','-').replace('-','')

        # Parse start and end position of turnover regions
        if start_pos_allele1_TO_region != end_pos_allele1_TO_region:
            position_to_region1 = str(start_pos_allele1_TO_region) + '-'+ str(end_pos_allele1_TO_region)
        if start_pos_allele2_TO_region != end_pos_allele2_TO_region:
            position_to_region2 = str(start_pos_allele2_TO_region) + '-'+ str(end_pos_allele2_TO_region)
        
        # if length TO position is 0 or 1, just the start position is taken into account  
        if start_pos_allele1_TO_region == end_pos_allele1_TO_region or start_pos_allele1_TO_region-1 == end_pos_allele1_TO_region:
            position_to_region1 = str(start_pos_allele1_TO_region)
        if start_pos_allele2_TO_region == end_pos_allele2_TO_region or start_pos_allele2_TO_region-1 == end_pos_allele2_TO_region:
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
    First, all output files including the headers are created, 5 in total. After the 
    analysis of each read, which is categorized, the data (at least the read name) is
    written into the correct output file. Ultimately, the metadata consisting of 
    category quantities is added to its output file.

    Args:
        read_name (str) = name of read
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
        Adds the non hybrid reads to the output file, if the read consensus has 1 mismatch, then
        it is added in a note.
        
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
        Adds the zero reads (reads that have 0 mismatches for multiple alleles) to the output file.
        The read name and a note with the alleles with zero mismatches are added.
        
        Args:
            note (str): The alleles with 0 mismatches are noted here
        Returns:
            -
        """
        print ('Read with 0 mismatches for multiple alleles:', self.read_name)

        # add reads that have 0 mismatches for multiple alleles
        with open(CreateOutput.output_file_zero_reads, 'a') as db_file:
            db_file.write(self.read_name + '\t' + note + '\n') 
    
    def hybrid_read_more_switches(self):
        """
        Adds the hybrid reads with more than 1 switch to the output file, these read pairs gave too many
        mismatches which were not indicative. Only the read name is added.
        
        Args:
            -
        Returns:
            -    
        """

        print ('Read with more switches: ', self.read_name)
        
        # add hybrid reads with more switches  
        with open(CreateOutput.output_file_more_switches, 'a') as db_file:
            db_file.write(self.read_name + '\n') 
   
    def hybrid_read_1_switch(self, allele_name, pos_read1_allele, pos_read2_allele, allele_read1_mismatches, allele_read2_mismatches, allele_consensus_mismatches, turn_over_region, pos_to_region, read_artefacts):
        """
        Adds all hybrid reads with 1 switch data to the output file
        
        Args:
            allele_name (str): Allele name
            pos_read1_allele (str): Start and stop position for read 1 relative to the allele
            pos_read2_allele (srt): Start and stop position for read 2 relative to the allele
            allele_read1_mismatches (str): Total number of mismatches in read 1
            allele_read2_mismatches (str): Total number of mismatches in read 2
            allele_consensus_mismatches (str): Total number of mismatches in read consensus
            turn_over_region (str): Sequence of turnover region
            pos_to_region (str): Start and stop position of turnover region relative to the allele
            read_artefacts (int): Total number of read artefacts, mutual and alternately mismatches (max. 4)
        Returns:
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

        Returns:
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



def main():
    """
    This is the main function of the script and calls all methods according to the sequence diagram. All reads are monitored and counted. 
    If one does not met the set requirements then it is skipped by using 'continue' (these reads are also monitored). 
        
    Args:
        -
    Returns:
        -
    """

    input_file = argv[1]

    # Create all output files
    CreateOutput.prep_output_files(input_file)
  
    # Parse input file
    with open (input_file) as file_object:
            input_file = file_object.read()
    all_data, allele_data = ParseInput.collect_all_data(input_file)
    all_allele_combinations = ParseInput.get_allele_combinations(allele_data)

    # Track all reads
    incorrect_aligned_reads = 0
    rejected_read_count = 0
    non_hybrid_count = 0
    zero_count = 0
    more_switches_count = 0
    one_switch_hybrid_count = 0
    read_counter = 1

    # Loop through each read pair
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
        
        # Get allele data
        allele_data = read_info[4:]
       
        # Perform checks for read 1
        R1_read = Read(read1_seq, read1_aligned_seq, allele_data)
        check_alignment = R1_read.check_alignment()
        if check_alignment == False:  # Check if alignement correct
            incorrect_aligned_reads += 1
            continue        
        R1_alignment_after_first_check = R1_read.apply_qv(read1_qv)
        R1_alignment_after_second_check = R1_read.check_read_artefacts(R1_alignment_after_first_check)

        # Perform checks for read 2
        R2_read = Read(read2_seq, read2_aligned_seq, allele_data)
        R2_alignment_after_first_check = R2_read.apply_qv(read2_qv)
        R2_alignment_after_second_check = R2_read.check_read_artefacts(R2_alignment_after_first_check)

        # Check if read pair met the requirements
        R1_and_R2 = ReadPair(R1_alignment_after_second_check, R2_alignment_after_second_check, read1_seq, read2_seq)

        # Adjust requirement values here:
        min_read_length = 80
        N_quantity = 5

        approve_reads = R1_and_R2.check_read_pair(min_read_length, N_quantity)
        if approve_reads == True:
            print ('Paired-end read is accepted')
        if approve_reads == False:
            rejected_read_count += 1
            print ('Paired-end read is rejected')
            continue

        ###########
        ###########  Mismatches per read
        ###########
        # Get mismatches per read for each allele
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

        ### For zero reads (multiple alleles with 0 mismatches)
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
        
        # Create read consensus
        alignment_read_consensus = R1_and_R2.create_read_consensus()
        consensus_read = Read.classmethod_for_non_read(alignment_read_consensus, allele_data)
        mismatch_dict_read_con, mismatch_dict_read_con_ex = consensus_read.get_mismatches(alignment_read_consensus)
        mismatch_dict_read_con_sorted = sorted(mismatch_dict_read_con.items(), key=lambda kv: kv[1])
            
        ### For non hybrid reads, if read consensus has 0 or 1 mismatches
        if mismatch_dict_read_con_sorted[0][1][0] == 0 or mismatch_dict_read_con_sorted[0][1][0] == 1:
            allele_match = mismatch_dict_read_con_sorted[0][0]
            if mismatch_dict_read_con_sorted[0][1][0] == 1:
                note = 'Read consensus has 1 mismatch'
            non_hybrid_count += 1
            read_output = CreateOutput(read_name)
            read_output.non_hybrid_read(allele_match, note)
            continue
        
        # Print all mismatch information
        R1_read.print_mismatches('First read', R1_mismatch_dict_ex)
        R2_read.print_mismatches('Second read', R2_mismatch_dict_ex)
        consensus_read.print_mismatches('Read consensus', mismatch_dict_read_con_ex)


        ###########
        ###########  Determine number of switches for all allele combinations
        ###########
        more_switches = True
        # Loop through all allele combinations
        for allele_combo in all_allele_combinations:
            allele1 = allele_combo[0]
            allele2 = allele_combo[1]
         
            per_allele_info = CheckAlleleCombination(alignment_read_consensus, allele_combo, allele_data)

            # Create indicator string
            allele_seq_list = per_allele_info.create_indicator_string()

            # Apply first and second check, and update indicator string
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

            # Generate all data if allele combo resulted in a 1 switch indicator string
            if nr_of_switches == 1:
                more_switches = False

                # Print 1 switch pre data
                per_allele_info.print_1_switch_alleles()
                read1_pos_dict = R1_read.get_relative_position()
                read2_pos_dict = R2_read.get_relative_position()
                
                # Get read positions
                final_to_region = GetOneSwitchData(allele1, allele2)
                pos_read1_allele1, pos_read2_allele1, pos_read1_allele2, pos_read2_allele2 = final_to_region.get_read_position(read1_pos_dict, read2_pos_dict)
               
                # Get turnover region sequence and positions
                turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = final_to_region.prep_for_turnover_position(start_turn_pos, end_turn_pos, allele_seq_list, allele_data)
                TO1_seq = Read.classmethod_for_non_read(turn_over_region1_for_pos, seq_list_allele1)
                TO2_seq = Read.classmethod_for_non_read(turn_over_region2_for_pos, seq_list_allele2)

                TO_allele1_dict = TO1_seq.get_relative_position()
                TO_allele2_dict = TO2_seq.get_relative_position()

                pos_to_region1, pos_to_region2, turn_over_region1, turn_over_region2 = final_to_region.get_TO_position(TO_allele1_dict, TO_allele2_dict, turn_over_region1_for_pos, turn_over_region2_for_pos)
                
                # Print 1 switch extended data              
                GetOneSwitchData.print_TO_output(allele1, pos_read1_allele1, pos_read2_allele1, turn_over_region1, pos_to_region1)
                GetOneSwitchData.print_TO_output(allele2, pos_read1_allele2, pos_read2_allele2, turn_over_region2, pos_to_region2)
                
                # Parse 1 switch extended data and add it to output file
                allele1_read1_mismatches = str(R1_mismatch_dict[allele1][0])
                allele1_read2_mismatches = str(R2_mismatch_dict[allele1][0])
                allele2_read1_mismatches = str(R1_mismatch_dict[allele2][0])
                allele2_read2_mismatches = str(R2_mismatch_dict[allele2][0])
                allele1_consensus_mismatches = str(mismatch_dict_read_con[allele1][0])
                allele2_consensus_mismatches = str(mismatch_dict_read_con[allele2][0])

                read_output = CreateOutput(read_name)
                read_output.hybrid_read_1_switch(allele1, pos_read1_allele1, pos_read2_allele1, allele1_read1_mismatches, allele1_read2_mismatches, allele1_consensus_mismatches, turn_over_region1, pos_to_region1, number_of_artefacts)
                read_output.hybrid_read_1_switch(allele2, pos_read1_allele2, pos_read2_allele2, allele2_read1_mismatches, allele2_read2_mismatches, allele2_consensus_mismatches, turn_over_region2, pos_to_region2, number_of_artefacts)

        # If at least one the allele combinations resulted in indicator string with 1 switch
        if more_switches == False:
            print ('Hybrid read with 1 switch: ', read_name)
            one_switch_hybrid_count += 1

        # If non of the allele combinations resulted in indicator string with 1 switch, then we found a
        # hybrid read with more switches (too many mismatches).  
        if more_switches == True:
            more_switches_count += 1
            read_output = CreateOutput(read_name)
            read_output.hybrid_read_more_switches()
       
    # Output metadata
    total_nr_of_reads = incorrect_aligned_reads + rejected_read_count + non_hybrid_count + zero_count + more_switches_count + one_switch_hybrid_count
    CreateOutput.metadata(incorrect_aligned_reads, rejected_read_count, non_hybrid_count, zero_count, more_switches_count, one_switch_hybrid_count, total_nr_of_reads)

if __name__ == "__main__":
    main()