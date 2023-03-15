from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from collections import defaultdict
from collections import Counter
import argparse

def get_prots_and_strains(in_file):
    all_prots = []
    all_strains = []
    with open(in_file, 'r') as infile:
        for line in infile:
            prot = line.split()[2]
            strain = line.split()[0].split('_')[0]
            all_strains.append(strain)
            all_prots.append(prot)
    prot_counts = Counter(all_prots)
    all_strains = set(all_strains)
    all_prots = [prot for prot, count in prot_counts.items() if count > 0.3 * len(all_strains)]
    return all_prots

def main():
    parser = argparse.ArgumentParser(
        description=('Concatenate a series of individually aligned genes into one alignment'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-s', '--singletons', help='file containing singletons')
    parser.add_argument('-m', '--min_frac', help='Mininum fraction of column to be ungapped')
    parser.add_argument('-o', '--out_name', help='Output concatenated alignment filename')
    args = parser.parse_args()

    singletons = get_prots_and_strains(args.singletons)
    filelist = [os.path.join("output/03.phylo_tree/ribo_aligned/", "{}_aligned.fasta".format(prot)) for prot in singletons]
    min_frac = float(args.min_frac)
    out_name = args.out_name
    sequence_dict = defaultdict(str)
    final_seqs = []
    names1 = []
    all_strains = []
    num = 0

    for files in filelist:
        strains = [seqs.id for seqs in SeqIO.parse(files, 'fasta')]
        all_strains += strains
    all_strains = set(all_strains)

    for files in filelist:
        prot_strains = []
        for seqs in SeqIO.parse(files, 'fasta'):
            taxon = seqs.id
            prot_strains.append(taxon)
            sequence_dict[taxon] += str(seqs.seq)
            ali_length = len(seqs)
        gap_seq = "-" * ali_length
        for taxa in (all_strains - set(prot_strains)):
            sequence_dict[taxa] += gap_seq

    for name, sequence in sequence_dict.items():
        seq = Seq(sequence)
        temp_record = SeqRecord(seq, id = name, description='')
        final_seqs.append(temp_record)

    SeqIO.write(final_seqs, out_name, 'fasta')

    columns_to_filter = []

    aln = AlignIO.read(out_name, 'fasta')
    new_align = aln
        
    AlignIO.write(new_align, out_name, 'fasta')
    
if __name__ == '__main__':
    main()