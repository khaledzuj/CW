#this program search for a desierd gene from an input file (single file) the user can select a desierd parameter such as protien sequince length and specific amino acid numper inside the sequince
# for runing this code use the formula  python 3 "file_location" minimum_length maximum_lingth minimum_cys_count max_cys_count

import sys
file= (sys.argv[1])
def translate_seq(data):
    """this function translates the nucleotides sequence in data file into amino acid"""
    # create a dictionary containing the nucleotides with their translated amino acid
    genetic_code={'ATT':'I','ATC':'I','ATA':'I','CTT':'L','CTC':'L','CTA':'L','CTG':'L','TTA':'L','TTG':'L','GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A','GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V','TTC':'F', 'TTT':'F','ATG':'M','TGC':'C', 'TGT':'C','GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G','CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P','ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T','AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S','TAC':'Y', 'TAT':'Y','TGG':'W','CAA':'Q', 'CAG':'Q','AAC':'N', 'AAT':'N','CAC':'H', 'CAT':'H','GAA':'E', 'GAG':'E','GAC':'D', 'GAT':'D','AAA':'K', 'AAG':'K','AGA':'R', 'AGG':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R','TAA':'_','TAG':'_', 'TGA':'_'}
    protien_seq="" #empty string to add the amino acid to it
    for i in data: # translation process
        for x in range(0, len(data)-(3+len(data)%3),3):# divide the sequence into 3 nucleotides pairs for the translation prosses
            protien_seq+=genetic_code[data[x:x+3]] # add the translated amino acid into the empty string
    return protien_seq


def get_protien_seq(file,len_range_start,len_range_end,cys_count_range_start,cys_count_range_end):
    """return the function to list of possible sequences that match the desired gene  criteria """
    with open(file, 'r') as data: # to parse the nucleotide sequence file to translate it
        for lines in data:
            if not lines.startswith('>'):
                protien_seq=translate_seq(lines[lines.find("ATG"):])# translate the amino acid and start reading from the first M

    #the second step is to search for the desired gene
    k= 1000 #the size of the protein sequence imported each time to study
    possiple_gene=[] #empty list to add on
    for i in range(0,len(protien_seq),900):# select the methods in which the data introduced to the searching prosses each time we shift to take 900 new amino acid sequence
        gene=protien_seq[i:i+k].split("_")# split the sequence in the ending codon places
        for x in gene: # look for our desierd gene
            if 'M' in x and len(x)>=len_range_start and len(x)<=len_range_end: # primary selection, select the fragments that have starting codon with the desired length
                target_gene=x
                if target_gene.count("C")>=cys_count_range_start and target_gene.count("C") <=cys_count_range_end: # secondary selection, select the gene fragment that contains the specified Cys range.
                    z= target_gene
                    if z not in possiple_gene: # add the unique potential genes into the empty list
                        possiple_gene.append(z)
    return possiple_gene # return the function into the possible gene list



def blast_protien_seq(data,database_type,base):
    """blast the genes sequince debend on (database_type,base) ,and produce an xML file"""
    import Bio
    from Bio.Blast import NCBIWWW # import nesseary pakeges
    result_handle = NCBIWWW.qblast(database_type, base, str(data)) # this code used to blast the sequences
    blast_result = result_handle.read()
    print (blast_result)
    return blast_result

blast_protien_seq(get_protien_seq(file,int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])),"blastp","nr")
