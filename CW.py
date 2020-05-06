fasta= ('/home/khaled/Downloads/KJ_problem/Cgr-B-H4-1G22_corrected_assembled_S-locus.fas')

def translate(data):
    """this function translates the nucleotides sequence in the data file into amino acid"""
    # create a dictionary containing the nucleotides with their translated amino acid
    genetic_code={'ATT':'I','ATC':'I','ATA':'I','CTT':'L','CTC':'L','CTA':'L','CTG':'L','TTA':'L','TTG':'L','GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A','GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V','TTC':'F', 'TTT':'F','ATG':'M','TGC':'C', 'TGT':'C','GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G','CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P','ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T','AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S','TAC':'Y', 'TAT':'Y','TGG':'W','CAA':'Q', 'CAG':'Q','AAC':'N', 'AAT':'N','CAC':'H', 'CAT':'H','GAA':'E', 'GAG':'E','GAC':'D', 'GAT':'D','AAA':'K', 'AAG':'K','AGA':'R', 'AGG':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R','TAA':'_','TAG':'_', 'TGA':'_'}
    protien_seq="" #empty string to add the amino acid to it
    for i in data: # translation process
        for x in range(0, len(data)-(3+len(data)%3),3):# divide the sequence into 3 nucleotides pairs for the translation prosses
            protien_seq+=genetic_code[data[x:x+3]] # add the translated amino acid into the empty string
    return protien_seq


def get_gene(file,len_range_start,len_range_end,cys_count_range_start,cys_count_range_end):
    """return the function to list of possible sequences that match the desired gene  criteria """
    with open(file, 'r') as data: # to parse the nucleotide sequence file to translate it
        for lines in data:
            if not lines.startswith('>'):
                protien_seq=translate(lines[lines.find("ATG"):])# translate the amino acid and start reading from the first M

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



def blast_gene(data,database_type,base):
    """blast the genes sequince debend on (database_type,base) ,and produce an xML file"""
    import Bio
    from Bio.Blast import NCBIWWW # import nesseary pakeges
    result_handle = NCBIWWW.qblast(database_type, base, str(data)) # this code used to blast the sequences
    blast_result = result_handle.read()
    return (blast_result)


gene_blast= blast_gene(get_gene(file,60,150,6,9),"blastp","nr")
