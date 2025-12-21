import pandas as pd
import argparse

def prep_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["SE"] = pd.to_numeric(df["SE"], errors="coerce")
    df["Z"] = df["BETA"] / df["SE"]
    df.rename(columns={"FRQ": "MAF"}, inplace=True)
    return df

def main():
    parser = argparse.ArgumentParser(description="Add Z-scores to GWAS summary stats (TSV).")
    parser.add_argument("-i", "--input", required=True, help="Input GWAS summary stats TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()
    df = pd.read_csv(args.input, sep="\t")
    df = prep_data(df)
    df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()