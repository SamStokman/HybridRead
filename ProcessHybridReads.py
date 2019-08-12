"""
10-07-'19

Script that creates the tables with  number of hybrid reads (1 switch) for each allele combo. 
Use output file hybrid_read_1_switch.txt from SelectHybridReads.py

Command line: python3 ProcessHybridReads.py hybrid_1_switch_data_hla_a.txt hybrid_1_switch_data_hla_b.txt hybrid_1_switch_data_hla_b.txt 

"""

from sys import argv
 
 
def create_combo_output(input_file1, input_file2, input_file3):
    """
    Function that processes three input files (from SelectHybridReads.py) with hybrid read 1 switch data. Per allele
    combination the number of reads is determined, per file but also for the three file together. If the read is present
    more than once, then it is ignored. This data is written into a output file.

    Args:
        input_file1 (str): contains all hybrid reads, and their allele data, that were initially assigned to HLA-A.
        input_file2 (str): contains all hybrid reads, and their allele data, that were initially assigned to HLA-B.
        input_file3 (str): contains all hybrid reads, and their allele data, that were initially assigned to HLA-C.
    Returns:
        -
    """
   
    collect_all_data = []
    approved_read_list = []
    input_file1 = input_file1.split('\n')
    for line in input_file1[1:-1]:
        line = line.split('\t')
        read_name = line[0]
        artefacts = int(line[7])
        if artefacts == 0:   # change here the number of allowed artefacts for HLA-A
            approved_read_list += [read_name]
        collect_all_data += [['A'] + line]

    input_file2 = input_file2.split('\n')
    for line in input_file2[1:-1]:
        line = line.split('\t')
        read_name = line[0]
        artefacts = int(line[7])
        if artefacts == 0:    # change here the number of allowed artefacts for HLA-B
            approved_read_list += [read_name]
        collect_all_data += [['B'] + line]
  
    input_file3 = input_file3.split('\n')
    for line in input_file3[1:-1]:
        line = line.split('\t')
        read_name = line[0]
        artefacts = int(line[7])
        if artefacts == 0:   # change here the number of allowed artefacts for HLA-C
            approved_read_list += [read_name]
        collect_all_data += [['C'] + line]

    read_dict = {}
    paired_end_read = False
    for line in collect_all_data:
        read_name = line[1]
        if read_name in approved_read_list and approved_read_list.count(read_name) == 2:    # only use reads if they only occur once (no double reads)
            hla_type = line[0]
            allele = str(line[2])
            if paired_end_read == False:
                allele1 = allele
                paired_end_read = True
                continue
            if paired_end_read == True:
                paired_end_read = False
                allele_combo = str(allele1) + ', ' + str(allele)
                allele_combo_v = str(allele) + ', ' + str(allele1)
           
                if allele_combo in read_dict.keys() or allele_combo_v in read_dict.keys():
                    read_dict[allele_combo] += [[hla_type, read_name]]
                
                if allele_combo not in read_dict.keys() and allele_combo_v not in read_dict.keys():
                    read_dict[allele_combo] = [[hla_type, read_name]]
                continue
  
    read_list_sorted = sorted(read_dict.items(), key=lambda x: x[1])
    total_quantity_dict = {}
    for allele, reads in read_list_sorted:
        nr_of_reads = len(reads)
        allele = allele.split(',')
        allele1 = allele[0].strip(' ')
        allele2 = allele[1].strip(' ')
        if allele1 in total_quantity_dict.keys():
            total_quantity_dict[allele1] += nr_of_reads
        if allele2 in total_quantity_dict.keys():
            total_quantity_dict[allele2] += nr_of_reads
        if allele1 not in total_quantity_dict.keys():
            total_quantity_dict[allele1] = nr_of_reads
        if allele2 not in total_quantity_dict.keys():
            total_quantity_dict[allele2] = nr_of_reads

        total_quantity_dict[allele1]
        print (allele, '$', len(reads)) 

    # Create output file with read counts, ratios and total read count per allele.
    output_file_name = 'hybrid_read_summary.txt'
    with open(output_file_name, 'w') as db_file: 
        db_file.write('Allele combination\tHLA-A\tHLA-B\tHLA-C\tTotal\n')
        
        total_A_count = 0
        total_B_count = 0
        total_C_count = 0
        total_ABC_count = 0
        for read_info in read_list_sorted:
            alleles = read_info[0]
            A_count = 0
            B_count = 0
            C_count = 0
            ABC_count = 0
            all_reads = read_info[1]
            for read in all_reads:
                hla_type = read[0]
                if hla_type == 'A': 
                    A_count += 1
                if hla_type == 'B': 
                    B_count += 1
                if hla_type == 'C': 
                    C_count += 1
                ABC_count = A_count + B_count + C_count
            total_A_count += A_count
            total_B_count += B_count
            total_C_count += C_count
            
            db_file.write(str(alleles) +'\t'+ str(A_count) + '\t' + str(B_count) + '\t' + str(C_count) +'\t' + str(ABC_count) + '\n')
        total_ABC_count += total_A_count + total_B_count + total_C_count
        db_file.write('Total\t'+ str(total_A_count) + '\t' + str(total_B_count) + '\t' + str(total_C_count) +'\t' + str(total_ABC_count) + '\n')
        db_file.write('$$$\n')
        db_file.write('Allele combination\tHLA-A\tHLA-B\tHLA-C\n')
        for read_info in read_list_sorted:
            alleles = read_info[0]
            A_count = 0
            B_count = 0
            C_count = 0
            ABC_count = 0
            all_reads = read_info[1]
            for read in all_reads:
                hla_type = read[0]
                if hla_type == 'A': 
                    A_count += 1
                if hla_type == 'B': 
                    B_count += 1
                if hla_type == 'C': 
                    C_count += 1
                ABC_count = A_count + B_count + C_count
            db_file.write(str(alleles) +'\t'+ str(round(A_count/ABC_count, 2)) + '\t' + str(round(B_count/ABC_count, 2)) + '\t' + str(round(C_count/ABC_count, 2)) + '\n')
        
        db_file.write('$$$\n')
        db_file.write('Hybrid read quantities per allele\n')
        for key, value in total_quantity_dict.items():
            db_file.write(str(key) +'\t'+ str(value) + '\n')        


if __name__ == "__main__":

    
    with open(argv[1]) as file_object1:
        input_file1 = file_object1.read()
        
        
    with open(argv[2]) as file_object2:
        input_file2 = file_object2.read()
    
    with open(argv[3]) as file_object3:
        input_file3 = file_object3.read()
      
      
    create_combo_output(input_file1, input_file2, input_file3)

