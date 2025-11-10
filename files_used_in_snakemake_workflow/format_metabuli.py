from pathlib import Path

import typer


def main(filepath_clas: Path, filepath_report: Path):
    # Create mapping of tax_id to taxonomy based on the metabuli_report
    curr_taxonomy = []
    tax_id2taxonomy = {}
    with open(filepath_report) as f:
        for line in f:
            id_line = line.split("\t")[4].strip()
            line = line.split("\t")[5].rstrip()
            curr_space = len(line.split(" ")) - 1
            while len(curr_taxonomy) - 1 < curr_space:
                curr_taxonomy.append("")
            curr_taxonomy[curr_space] = line

            to_add = []
            for i, entry in enumerate(curr_taxonomy):
                if i > curr_space:
                    break
                to_add.append(entry)

            tax_id2taxonomy[id_line] = ";".join([x.strip() for x in to_add if x != ""])

    # Convert the metabuli_classifications to taxonomy using the tax_id2taxonomy mapping
    print("contigs\tpredictions")
    with open(filepath_clas, "r") as f:
        for line in f:
            sp_line = line.split("\t")
            search_id = sp_line[2]
            tax = tax_id2taxonomy[search_id]
            if tax == "unclassified":
                tax = ""
            print(f"{sp_line[1]}\t{tax[5:]}")


if __name__ == "__main__":
    typer.run(main)
