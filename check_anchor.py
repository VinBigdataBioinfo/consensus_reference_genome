from get_sequences import main
import pandas as pd
import os

def read_sequences():
    header = ["index", "name", "chrLength", "chromStart", "chromEnd", "frag", "fragStart", "fragEnd", "fragStrand", "fragLength"]
    df = pd.read_csv("sequence.csv", names=header, sep="\t")
    return df


def get_anchors(file_name):
    anchors_list = list()
    with open(file_name) as f:
        lines = f.readlines()

    number_of_anchors = int((len(lines) - 1) / 6) -1
    chr = lines[0].replace("\n","")

    for i in range(number_of_anchors):
        frag = lines[(i+2)*6].replace("\n","").split(": ")[-1].replace("\"","")
        chromStart, chromEnd = lines[(i+1)*6 + 1].replace("\n","").split()[-1].split("..")
        anchors_list.append([chr, frag, int(chromStart), int(chromEnd)])
    
    return anchors_list


def get_all_anchors():
    all_anchors_list = list()
    dir = "crawl_data/"
    for file in os.listdir(dir):
        all_anchors_list.extend(get_anchors(dir + file))
    
    return all_anchors_list


def main():
    df = read_sequences()
    all_anchors_list = get_all_anchors()
    anchor_list_df = pd.DataFrame(all_anchors_list, columns=["chrom", "frag", "chromStart", "chromEnd"])

    is_anchor = list()
    for seq in df.itertuples():
        df_check = anchor_list_df[(anchor_list_df['chrom']==seq.name) & 
                        (anchor_list_df["frag"]==seq.frag) & 
                        (anchor_list_df['chromStart']==seq.chromStart+1) &
                        (anchor_list_df['chromEnd']==seq.chromEnd)]
        is_anchor.append(len(df_check))
        
    df["anchor"] = is_anchor
    df.to_csv("sequence.isAnchor.csv", index=False, sep="\t")


if __name__=="__main__":
    main()