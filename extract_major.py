from pysam import VariantFile
import pandas as pd

def main():
    anchor_df = pd.read_csv("sequence.isAnchor.hg38pos.csv", sep="\t")
    bcf_in = VariantFile("majorSet.fix_delins.vcf.gz")
    bcf_out = VariantFile("majorSet.fix_delins.extract_anchor.vcf", "w", header=bcf_in.header)

    for anchor in anchor_df.itertuples():
        print(anchor)
        for rec in bcf_in.fetch(anchor.hg38_chrom, anchor.anchor_hg38_start, anchor.anchor_hg38_end):
            changed_rec = rec.copy()
            changed_rec.contig = anchor.name
            if (anchor.hg38_fragLength < anchor.fragLength) and (anchor.fragStart < anchor.hg38_fragStart):
                changed_rec.pos = rec.pos - anchor.anchor_hg38_start + anchor.chromStart + (anchor.hg38_fragStart - anchor.fragStart)
            else:
                changed_rec.pos = rec.pos - anchor.anchor_hg38_start + anchor.chromStart
            
            bcf_out.write(changed_rec)


if __name__=="__main__":
    main()