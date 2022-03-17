import pandas as pd

VERSION = "9.1"
PATH_TO_IOC_NAMES = f"./data/external/IOC_Names_File_Plus-{VERSION}.xlsx"
PATH_TO_OUT = f"./data/interim/taxonomy-{VERSION}.csv"


def main():
    df = pd.read_excel(PATH_TO_IOC_NAMES)
    cols = ["Rank", "Scientific Name"]

    required_ranks = ['TAXON', 'ORDER', 'Family', 'Genus', 'Species']
    current_taxonomy = {r: None for r in required_ranks}

    species = []
    for _, row in df[cols].iterrows():
        rank = row["Rank"]
        sciname = row["Scientific Name"]

        if rank in current_taxonomy:
            current_taxonomy[rank] = sciname
            if rank == "Species":
                species.append(current_taxonomy.copy())

    taxa_df = pd.DataFrame(species)
    taxa_df["TAXON"] = taxa_df.TAXON.str.capitalize()
    taxa_df["ORDER"] = taxa_df.ORDER.str.removeprefix("ORDER ").str.capitalize()
    taxa_df["Family"] = taxa_df.Family.str.removeprefix("Family ")
    taxa_df.columns = taxa_df.columns.str.capitalize()

    taxa_df.to_csv(PATH_TO_OUT, index=None)


if __name__ == "__main__":
    main()
