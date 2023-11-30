import requests
import pandas as pd
from tqdm import tqdm

def read_hg38alt():
    alt_contigs_list = list()
    with open("Homo_sapiens_assembly38.fasta.64.alt") as f:
        for line in f.readlines():
            cols = line.split("\t")
            if "_alt" in cols[0]:
                alt_contigs_list.append(cols[0])

    return alt_contigs_list


def get_contig_seq(chr):
    url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=gold;chrom={chr}"
    return requests.get(url).json()


def parse_info(index, contig_json):
    seqs = list()
    name = contig_json["chrom"]
    length = contig_json["end"]
    for seq in contig_json["gold"]:
        seqs.append([index+1, name, length, seq["chromStart"], seq["chromEnd"], # ALT CONTIGS
                seq["frag"], seq["fragStart"], seq["fragEnd"], seq['strand'], seq["fragEnd"]-seq["fragStart"], # Accessions coordinator
                ]) 

    return seqs


def write_list(list):
    df = pd.DataFrame(list)
    df.to_csv("sequence.csv", index=False, header=False, sep="\t")


def main():
    alt_contigs_list = read_hg38alt()
    all_seq_list = list()
    for index, chr in tqdm(enumerate(alt_contigs_list)):
        contig_json = get_contig_seq(chr)
        seqs = parse_info(index, contig_json)
        all_seq_list.extend(seqs)

    write_list(all_seq_list)


if __name__=="__main__":
    main()