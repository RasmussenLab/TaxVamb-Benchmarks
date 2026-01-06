# Script to ensure that the format of the kalmari and trembl database outputs are consistent with the fomrat taxvamb expect

import enum
from pathlib import Path

import typer


def fix_tax(string):
    when_to_start = {
        "k_": 0,
        "p_": 1,
        "c_": 2,
        "o_": 3,
        "f_": 4,
        "g_": 5,
        "s_": 6,
    }
    levels = ["" for x in range(len(when_to_start.keys()))]

    s_found = None
    split_string = string.split(";")

    if len(split_string) in [1, 2]:
        return string

    for i, tax in enumerate(split_string[2:]):
        if tax[:2] in when_to_start.keys():
            levels[when_to_start[tax[:2]]] = tax
        if tax[:2] == "s_":
            s_found = i + 2

    prev_non_empty = ""
    for i in reversed(range(len(levels))):
        tax = levels[i]
        if tax == "" and prev_non_empty != "":
            levels[i] = "LEVEL_" + str(i) + "_ADDED_FROM_" + prev_non_empty
        if tax != "":
            prev_non_empty = tax

    levels = split_string[:2] + levels
    if s_found is not None:
        levels = levels + split_string[s_found + 1 :]

    tax_list = ";".join([x for x in levels if x != ""])
    return tax_list


def format_tax_kalmari(string):
    formatted_tax = []
    last_tax = ""
    for i, tax in enumerate(string.split(";")):
        if (
            not tax.startswith(("-")) or i in [0, 1] or last_tax.startswith("s_")
        ):  # Filter out all taxnomy not all entries have. These start with -. But keep subspecies (which is last) and bacteria (which start with -)
            formatted_tax.append(tax)
        last_tax = tax
    formatted_tax = ";".join(formatted_tax)
    return formatted_tax


def test_all_levels_present(string):
    formated = format_tax_kalmari(string)

    how_many_levels_should_be = {
        "s": 7,
        "g": 6,
        "f": 5,
        "o": 4,
        "c": 3,
        "p": 2,
        "k": 1,
    }
    start = False
    count = 0
    first_tax = None
    for tax in reversed(formated.split(";")):
        if tax.startswith("-_"):
            start = False
        if not tax.startswith("-_"):
            start = True
            if first_tax is None:
                first_tax = tax
        if start:
            count += 1
    if (
        string == "-_cellular organisms;-_Bacteria"
        or string == "-_cellular organisms"
        or string == "-_root"
        or string == "-_Viruses"
        or string == "-_Viruses;-_unclassified viruses"
        or string == "-_cellular organisms;-_Archaea"
        or string == "-_cellular organisms;-_Eukaryota"
        or string == "-_Viruses;-_unclassified bacterial viruses"
    ):
        pass
    elif how_many_levels_should_be[first_tax[0]] != count:
        raise ValueError(f"""
                {how_many_levels_should_be}
                {count}
                {formated}
        """)
        # {how_many_levels_should_be}
        # {count}
        # print(f"""
        #         {string}
        # """)


file = "./head_cut.tsv"


def main(file: Path):
    levels = [
        "s",
        "g",
        "f",
        "o",
        "c",
        "p",
        "k",
    ]
    with open(file, "r") as f:
        for line in f:
            name, tax = line.split("\t")
            tax = tax.strip()
            if tax == "":
                tax = ""
            else:
                tax = fix_tax(tax)
                test_all_levels_present(tax)

                # print(tax)
            print(f"{name}\t{tax}")
            # print(name)

string = "-_cellular organisms;-_Bacteria;k_Pseudomonadati;p_Acidobacteriota;-_unclassified Acidobacteriota;s_Acidobacteriota bacterium"
fix_tax(string)

if __name__ == "__main__":
    typer.run(main)


