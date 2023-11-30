import pandas as pd
import requests
# from tqdm import tqdm

def read_sequences():
    # header = ["index", "name", "chrLength", "chromStart", "chromEnd", "frag", "fragStart", "fragEnd", "fragLength", "isAnchor"]
    df = pd.read_csv("sequence.isAnchor.csv", sep="\t")
    return df


def get_contig_seq(chr):
    url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=gold;chrom={chr}"
    return requests.get(url).json()


def update_mainChr_assemblyList(seq, main_chr, assembly_df):
    update_main_chr = main_chr
    update_assembly_df = assembly_df
    if update_main_chr not in seq.name:
        update_main_chr = seq.name.split("_")[0]
        update_assembly_list = get_contig_seq(update_main_chr)["gold"]
        update_assembly_df = pd.DataFrame(update_assembly_list)

    return update_main_chr, update_assembly_df


def calculate_region_on_hg38(seq, selected_seq_hg38):
    # if seq.fragStrand != selected_seq_hg38.strand:
    #     print("Strand different")

    if (seq.fragStart - selected_seq_hg38.fragStart) >= 0:
        anchor_hg38_start = selected_seq_hg38.chromStart + (seq.fragStart - selected_seq_hg38.fragStart)
    else:           
        anchor_hg38_start = selected_seq_hg38.chromStart

    if selected_seq_hg38.fragEnd - seq.fragEnd >= 0:
        anchor_hg38_end = selected_seq_hg38.chromEnd - (selected_seq_hg38.fragEnd - seq.fragEnd)
    else:
        anchor_hg38_end = selected_seq_hg38.chromEnd

    return anchor_hg38_start, anchor_hg38_end


def write_list(lst, file):
    header = ["pdindex", "index", "name", "chrLength", "chromStart", "chromEnd", "frag", "fragStart", "fragEnd", "strand", "fragLength", "isAnchor",
            "hg38_chrom", "hg38_chromStart", "hg38_chromEnd", "hg38_fragStart", "hg38_fragEnd", "hg38_fragLength",
            "anchor_hg38_start", "anchor_hg38_end"]
    df = pd.DataFrame(lst, columns=header)
    df.drop(["pdindex"], axis=1).to_csv(file, sep="\t", index=False)


def main():
    df = read_sequences()
    df = df[df["anchor"]==1]
    seq_full_info_list = list()
    
    main_chr = "chr0" # init main_chr not in alt contig name and it will update and get info from api
    assembly_df = pd.DataFrame()
    for seq in df.itertuples():
        main_chr, assembly_df = update_mainChr_assemblyList(seq, main_chr, assembly_df)

        seq_hg38 = assembly_df[(assembly_df["frag"] == seq.frag) & 
                                (assembly_df["fragEnd"] > int(seq.fragStart)) & 
                                (assembly_df["fragStart"] < int(seq.fragEnd)) & 
                                (assembly_df["fragStart"] <= int(seq.fragStart)) &
                                (assembly_df["fragEnd"] >= int(seq.fragEnd))]
        if len(seq_hg38) != 1:
            seq_hg38 = assembly_df[(assembly_df["frag"] == seq.frag) & 
                                    (assembly_df["fragEnd"] > int(seq.fragStart)) &
                                    (assembly_df["fragStart"] < int(seq.fragEnd)) & 
                                    ((assembly_df["fragStart"] <= int(seq.fragStart)) | (assembly_df["fragEnd"] >= int(seq.fragEnd)))]

        if len(seq_hg38) == 0:
            print("Frag is anchor but dont exist on main chr assemble:", end=" ")
            print(seq.name, seq.frag, seq.fragStart, seq.fragEnd)
        # elif len(seq_hg38) == 1:
        #     selected_seq_hg38 = seq_hg38.iloc[0]
        #     # get pos if seq38 not contain all seq
        #     anchor_hg38_start, anchor_hg38_end = calculate_region_on_hg38(seq, selected_seq_hg38)

        #     update_seq = list(seq)
        #     update_seq.extend([selected_seq_hg38.chrom, selected_seq_hg38.chromStart, selected_seq_hg38.chromEnd, selected_seq_hg38.fragStart, selected_seq_hg38.fragEnd, selected_seq_hg38.fragEnd-selected_seq_hg38.fragStart])
        #     update_seq.extend([anchor_hg38_start, anchor_hg38_end])
        #     seq_full_info_list.append(update_seq)
        else:
            # print(seq.name, seq.frag, seq.fragStart, seq.fragEnd)
            # print(seq_hg38)
            max_isec_base = 0
            selected_seq_hg38 = seq_hg38.iloc[0]
            for s in seq_hg38.itertuples():
                anchor_hg38_start, anchor_hg38_end = calculate_region_on_hg38(seq, s)
                if anchor_hg38_end - anchor_hg38_start > max_isec_base:
                    selected_seq_hg38 = s
            anchor_hg38_start, anchor_hg38_end = calculate_region_on_hg38(seq, selected_seq_hg38)

            update_seq = list(seq)
            update_seq.extend([selected_seq_hg38.chrom, selected_seq_hg38.chromStart, selected_seq_hg38.chromEnd, selected_seq_hg38.fragStart, selected_seq_hg38.fragEnd, selected_seq_hg38.fragEnd-selected_seq_hg38.fragStart])
            update_seq.extend([anchor_hg38_start, anchor_hg38_end])
            seq_full_info_list.append(update_seq)
            
    
    write_list(seq_full_info_list, "sequence.isAnchor.hg38pos.csv")


if __name__=="__main__":
    main()