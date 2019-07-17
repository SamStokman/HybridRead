"""

Input needed: output file from align_all_reads.py

Script that finds hybrid reads


"""

from sys import argv

class ParseInput:

    
    def __init__(self, file_name):
        with open (file_name) as file_object:
            self.input = file_object.read()
        self.data_type = file_name[-9:-4]
    
    def collect_all_data(self):
        collect_all_data = []
        input_file = self.input.split('$$$')
        for read_data in input_file[1:2]:   # 1:-1
            temp_collect_list = []
            read_data = read_data.split('\n')
            for line in read_data:
                line = line.split('\t')
                if len(line) > 1:
                    temp_collect_list += [line]
            collect_all_data += [temp_collect_list]
        allele_data = collect_all_data[0][4:]
        
        return collect_all_data, allele_data

    @staticmethod
    def get_allele_combinations(allele_data):
        """
        This function creates all allele combinations (the type names) and returns them in a list

        allele_data: list, with all alleles and sequences
        """
        # for samples with 5 different alleles
        if len(allele_data) == 5:
            first_nr_list = [0,0,0,0,1,1,1,2,2,3]
            second_nr_list = [1,2,3,4,2,3,4,3,4,4] 
        
        # for samples with 6 different alleles
        if len(allele_data) == 6:
            first_nr_list = [0,0,0,0,0,1,1,1,1,2,2,2,3,3,4]
            second_nr_list = [1,2,3,4,5,2,3,4,5,3,4,5,4,5,5]

        all_alleles_list = []
        for allele, seq in allele_data:
           all_alleles_list += [allele]

        all_combinations_list = []
        for i in range(len(first_nr_list)):
            all_combinations_list += [[all_alleles_list[int(first_nr_list[i])], all_alleles_list[int(second_nr_list[i])]]]
        
        return (all_combinations_list)
        
class Read:
    def __init__(self, read_seq, read_qv, read_aligned_seq, allele_data):
        self.read_seq = read_seq
        self.read_qv = read_qv
        self.read_aligned_seq = read_aligned_seq
        self.allele_data = allele_data
    
    
    def check_alignment(self): 
        # Check if alignment is correct
        
        correct_alignment = True
        check_alignment = self.read_aligned_seq.replace('-', '')
        if self.read_seq != check_alignment:
            correct_alignment = False
        
        return correct_alignment
    
    def apply_qv(self):
        """
        Checks nucleotide quality values, if lower than a given value, then the nuclotide is replaced by a 'N'
    
        Returns the updated aligned read

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

    def get_mismatches(self, which_read, read_aligned_fully_checked):
        """
        Check for mismatches read vs all alleles (substitutions, insertions and deletions)

        allele_data: list, with all alleles and sequences
        read_seq: str, read sequence

        """
        read_seq = read_aligned_fully_checked
        allele_data = self.allele_data
    
        # get the read position (relative to the alleles) 
        start_relative_read_position = []
        start_relative_read_nucleotide = []
        nr_of_switches = 0
        nuc = read_seq[0]
        for i, char in enumerate(read_seq):
            if nuc == '-' and nuc != char:
                nr_of_switches += 1
                start_relative_read_position += [i]
                start_relative_read_nucleotide += [char]
            nuc = char

        # get the allele position (relative to the alleles), to check whether the read starts in front of the allele
    
        deletion_correction_dict = {} # if reads starts in front of allele
        for allele, seq in allele_data:  
            seq_string = seq
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
            count_read_insertion_for_deletion_correction = read_seq[start_relative_read_position[0]:start_allele_post]
            read_inserts = count_read_insertion_for_deletion_correction.count('-')
            deletion_correction_dict[allele] += read_inserts

    
        relative_read_position = []
        relative_read_nucleotide = []
        for i, char in enumerate(read_seq):
            if i >= start_relative_read_position[0] and i <= start_relative_read_position[-1]:
                relative_read_position += [i]
                relative_read_nucleotide += [char]
            if i >= start_relative_read_position[-1]+1 and char != '-':
                relative_read_position += [i]
                relative_read_nucleotide += [char]
   

        print ('Read info', which_read)
        print ('Length: \t\t', len(relative_read_nucleotide))
        print ('\nAllele\t\t\tSubstitutions\tInsertions\tDeletions\tTotal nr. or mismatches')

        mismatch_dict = {}
        for allele, seq in allele_data:
            substitutions = 0
            insertions = 0
            deletions = int(deletion_correction_dict[allele])
            mismatches = 0
            seq_string = seq
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
        

            print (allele, '\t\t', substitutions, '\t\t', insertions, '\t\t',deletions, '\t\t', mismatches) 
        print ('\n')
        return (mismatch_dict)

    def get_relative_position(self):
        """
        Determines the (consensus) read position for each allele and returns it in 
        a dict.
    
        """

        read_pos_dict = {}
        for allele, seq in  self.allele_data:
            allele_seq = seq
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
      
        # if turnover region has a length of 1 and the allele has '-' as nucleotide, select next nucleotide
        allele_name =  self.allele_data[0][0]
        if 'K' in self.read_aligned_seq and read_pos_dict[allele_name] == [[]]:
            read_pos_dict = self.__get_special_case_pos(read_pos_dict, allele_name)
        return read_pos_dict

    @staticmethod
    def __get_special_case_pos(read_pos_dict, allele_name):  #DEZE MOET NOG GETEST WORDEN
        # if turnover region has a length of 1 and the allele has '-' as nucleotide, select next nucleotide
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
    def __init__(self, read_consensus, allele_combo, allele_data):
        self.read_consensus = read_consensus
        self.allele_combo = allele_combo
        self.allele1 = allele_combo[0]
        self.allele2 = allele_combo[1]
        self.allele_data = allele_data
        self.indicator_string = ''
        self.number_of_artefacts = 0

    
    def create_indicator_string(self):

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
        # check if alleles have enough indicative mismatches (at least 2 per allele)
        accept_combo = True
        
        mismatch_indicator_string = self.indicator_string
        if mismatch_indicator_string.count('X') < 2 or mismatch_indicator_string.count('Y') < 2:
            accept_combo = False

        return accept_combo

    def check_mutual_SNPs(self):
        # check if alleles do not have too many mutual mismatches (max. 2)
        accept_combo = True
        mismatch_indicator_string = self.indicator_string
        self.number_of_artefacts += mismatch_indicator_string.count('M') 

        if mismatch_indicator_string.count('M') > 2:
            accept_combo = False

        return accept_combo
        
    def check_alternately_SNPs(self):
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

        return count_indicator_list
        
    def update_indicator_string(self, count_indicator_list):
        """
        Remove alternately artefacts
        
        
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
        print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
        print ('Allele combination:\t\t ', self.allele1.strip(' '), 'and', self.allele2)
        print ('Number of artefacts:\t\t ', self.number_of_artefacts, '\n')

class GetOneSwitchData():      

    def __init__(self, allele1, allele2):
             
        self.allele1 = allele1
        self.allele2 = allele2
        
        self.turn_over_region1 = ''
        self.turn_over_region2 = ''

    #def select_turnover_region(read_consensus, allele_combo, allele_info, orientation, read1_pos_dict, read2_pos_dict, position_consensus_dict):
    
    def get_read_position(self, read1_pos_dict, read2_pos_dict): 

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
        allele1 = self.allele1
        allele2 = self.allele2

        turn_over_region1_for_pos = ''
        turn_over_region2_for_pos = ''       
 
        if start_pos == end_pos +1:
            print ('Turnover sequence contains 0 nucleotides')

        for i, char in enumerate(allele_seq_list[0]):
            if start_pos == end_pos:    # for TO with length 1
                if i != start_pos:
                    turn_over_region1_for_pos += '-'
                if i == start_pos:
                    if char == '-':
                        char = 'Z'
                    turn_over_region1_for_pos += char
                    self.turn_over_region1 += char
            if start_pos == end_pos + 1:    # for TO with length 0
                char = 'K'     # maken we gewoon altijd een K van, deze sequence bestaat toch niet
                if i != start_pos:
                    turn_over_region1_for_pos += '-'
                if i == start_pos + 1:
                    turn_over_region1_for_pos += char
            if start_pos != end_pos + 1 and start_pos != end_pos:
                if i < start_pos:
                    turn_over_region1_for_pos += '-'
                if i > end_pos:
                    turn_over_region1_for_pos += '-'
                if i >= start_pos and i <= end_pos:
                    turn_over_region1_for_pos += char
                    self.turn_over_region1 += char
     
        for i, char in enumerate(allele_seq_list[1]):
            if start_pos == end_pos:    # for TO with length 1
                if i != start_pos:
                    turn_over_region2_for_pos += '-'
                if i == start_pos:
                    if char == '-':
                        char = 'Z'
                    turn_over_region2_for_pos += char
                    self.turn_over_region2 += char
            if start_pos == end_pos + 1:   # for TO with length 0
                char = 'K'
                if i != start_pos:
                    turn_over_region2_for_pos += '-'
                if i == start_pos + 1:
                    turn_over_region2_for_pos += char
            if start_pos != end_pos + 1 and start_pos != end_pos:
                if i == start_pos == end_pos + 1:
                    turn_over_region2_for_pos += char
                if i < start_pos:
                    turn_over_region2_for_pos += '-'
                if i > end_pos:
                    turn_over_region2_for_pos += '-'
                if i >= start_pos and i <= end_pos:
                    turn_over_region2_for_pos += char
                    self.turn_over_region2 += char
        
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
        
    def get_TO_position(self, TO_allele1_dict, TO_allele2_dict):
    
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

        # delete deletions
        turn_over_region1 = self.turn_over_region1.replace('-','')
        turn_over_region2 = self.turn_over_region2.replace('-','')

        position_to_region1 = str(start_pos_allele1_TO_region) + '-'+ str(end_pos_allele1_TO_region)
        position_to_region2 = str(start_pos_allele2_TO_region) + '-'+ str(end_pos_allele2_TO_region)

        return (position_to_region1, position_to_region2, turn_over_region1, turn_over_region2)
    
    @staticmethod
    def print_TO_output(allele, pos_read1_allele, pos_read2_allele, turn_over_region, pos_to_region):

        print ('Allele match:\t\t\t ', allele)
        print ('Read 1 position:\t\t ', pos_read1_allele)
        print ('Read 2 position:\t\t ', pos_read2_allele)
        print ('Turnover sequence length:\t ', len(turn_over_region))
        if turn_over_region == '':
            turn_over_region = '-'
        print ('Turnover sequence:\t\t ', turn_over_region)
        print ('Turnover region position:\t ', pos_to_region)
        print ('\n')
         

        
if __name__ == "__main__":
    
    input_file = argv[1]
        

    # Class ParseInput (Parse input file and get all allele combinations)
    msa_test_output = ParseInput(input_file)

    collect_all_data, allele_data  = msa_test_output.collect_all_data()
    all_combinations_list = msa_test_output.get_allele_combinations(allele_data)

    for read_info in collect_all_data:

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

        
        R1_mismatch_dict = R1_read.get_mismatches('First read', R1_alignment_after_second_check)
        
        # Perform checks for read 2
        R2_read = Read(read2_seq, read2_qv, read2_aligned_seq, allele_data)
        check_alignment = R2_read.check_alignment()
        if check_alignment == False:
            continue  
        R2_alignment_after_first_check = R2_read.apply_qv()
        R2_alignment_after_second_check = R2_read.check_read_artefacts(R2_alignment_after_first_check)

        R2_mismatch_dict = R2_read.get_mismatches('Second read', R2_alignment_after_second_check)


        # Check if read pair met the requirements
        R1_and_R2 = ReadPair(R1_alignment_after_second_check, R2_alignment_after_second_check, read1_seq, read2_seq)
        approve_reads = R1_and_R2.check_read_pair()
        if approve_reads == True:
            print ('Paired-end read is accepted')
        if approve_reads == False:
            print ('Paired-end read is rejected')

            continue
        
        # Create read consensus
        alignment_read_consensus = R1_and_R2.create_read_consensus()
        consensus_read = Read('', '', alignment_read_consensus, allele_data)
        mismatch_dict_read_con = consensus_read.get_mismatches('Read consensus', alignment_read_consensus)

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
                count_indicator_list = per_allele_info.check_alternately_SNPs()
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
                per_allele_info.print_1_switch_alleles()
            
            else:
                continue

            if nr_of_switches == 1:
                read1_pos_dict = R1_read.get_relative_position()
                read2_pos_dict = R2_read.get_relative_position()
                
                final_to_region = GetOneSwitchData(allele1, allele2)
                pos_read1_allele1, pos_read2_allele1, pos_read1_allele2, pos_read2_allele2 = final_to_region.get_read_position(read1_pos_dict, read2_pos_dict)
                turn_over_region1_for_pos, seq_list_allele1, turn_over_region2_for_pos, seq_list_allele2 = final_to_region.prep_for_turnover_position(start_turn_pos, end_turn_pos, allele_seq_list, allele_data)
                
                # Get TO position dicts
                TO1_seq = Read(seq_list_allele1, '', turn_over_region1_for_pos, allele_data)     
                TO2_seq = Read(seq_list_allele2, '', turn_over_region2_for_pos, allele_data)    
                TO_allele1_dict = TO1_seq.get_relative_position()
                TO_allele2_dict = TO2_seq.get_relative_position()

                pos_to_region1, pos_to_region2, turn_over_region1, turn_over_region2 = final_to_region.get_TO_position(TO_allele1_dict, TO_allele2_dict)
                GetOneSwitchData.print_TO_output(allele1, pos_read1_allele1, pos_read2_allele1, turn_over_region1, pos_to_region1)
                GetOneSwitchData.print_TO_output(allele2, pos_read1_allele2, pos_read2_allele2, turn_over_region2, pos_to_region2)
