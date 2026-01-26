"""
Format MMseqs2 Taxonomy Output (TrEMBL/Kalmari) for Taxvamb

This script standardizes and validates taxonomy lineages from MMseqs2 classifications using
TrEMBL or Kalmari databases, converting them into a format compatible with Taxvamb.

Purpose:
--------
MMseqs2 can classify contigs against multiple databases (TrEMBL, Kalmari, GTDB). The TrEMBL
and Kalmari databases produce taxonomy strings that may:
- Have incomplete taxonomic levels (missing ranks in the hierarchy)
- Include non-standard prefixes (e.g., "-_" for intermediate taxa)
- Lack proper formatting for downstream binning tools

This script performs three main operations:
1. fix_tax(): Reorganizes taxonomy strings to ensure all 7 standard ranks are present
   (k=kingdom, p=phylum, c=class, o=order, f=family, g=genus, s=species).
   Missing levels are filled with placeholder names derived from the nearest known rank.

2. format_tax_kalmari(): Filters out non-standard taxonomy entries (those starting with "-_")
   while preserving domain-level classifications and subspecies annotations.

3. test_all_levels_present(): Validates that the formatted taxonomy has the correct number
   of levels based on the deepest rank present, raising an error if validation fails.

Pipeline Integration:
--------------------
This script is called in the Snakemake workflow in two rules:
- run_taxvamb_kalmari: Processes MMseqs2 output from Kalmari database (column 5)
- run_taxvamb_trembl: Processes MMseqs2 output from TrEMBL database (column 9)

The workflow extracts the relevant columns using 'cut', then pipes through this script
to produce the final taxonomy file for Taxvamb.

Input:
------
- file: Tab-separated file with two columns (contig_name\ttaxonomy_string)
  The taxonomy_string should be semicolon-separated taxonomic lineage from MMseqs2

Output:
-------
Prints to stdout a tab-separated file with:
- Column 1: Contig name
- Column 2: Standardized and validated taxonomy lineage (semicolon-separated)

Special Cases:
-------------
- Empty taxonomy strings are preserved as empty
- Unclassified contigs at high levels (e.g., "-_cellular organisms;-_Bacteria") are
  allowed without validation
- Missing intermediate ranks are filled with "LEVEL_N_ADDED_FROM_<parent_taxon>"
"""

# Script to ensure that the format of the kalmari and trembl database outputs are consistent with the fomrat taxvamb expect

import enum
from pathlib import Path

import typer


def fix_tax(string):
    """
    Reorganize taxonomy string to ensure all 7 standard taxonomic ranks are present.

    This function fills in missing taxonomic levels with placeholder names derived from
    the nearest known parent taxon. For example, if genus is present but family is missing,
    it creates a placeholder family name.
    """
    # Map taxonomic rank prefixes to their position in the hierarchy (0=kingdom to 6=species)
    when_to_start = {
        "k_": 0,  # Kingdom
        "p_": 1,  # Phylum
        "c_": 2,  # Class
        "o_": 3,  # Order
        "f_": 4,  # Family
        "g_": 5,  # Genus
        "s_": 6,  # Species
    }
    # Initialize empty array for the 7 taxonomic levels
    levels = ["" for x in range(len(when_to_start.keys()))]

    s_found = None
    split_string = string.split(";")

    # If taxonomy is very short (1-2 entries), return as-is (likely root/domain only)
    if len(split_string) in [1, 2]:
        return string

    # Parse taxonomy string starting from position 2 (skip first 2 entries: root and domain)
    # Place each taxonomic rank into its correct position in the levels array
    for i, tax in enumerate(split_string[2:]):
        if tax[:2] in when_to_start.keys():
            levels[when_to_start[tax[:2]]] = tax
        if tax[:2] == "s_":
            s_found = i + 2  # Track position of species for later use

    # Fill in missing intermediate ranks by working backwards from deepest to shallowest
    # For missing ranks, create placeholder names based on the next available lower rank
    prev_non_empty = ""
    for i in reversed(range(len(levels))):
        tax = levels[i]
        if tax == "" and prev_non_empty != "":
            # Create placeholder: e.g., "LEVEL_4_ADDED_FROM_g_Escherichia" for missing family
            levels[i] = "LEVEL_" + str(i) + "_ADDED_FROM_" + prev_non_empty
        if tax != "":
            prev_non_empty = tax

    # Reconstruct final taxonomy: [root, domain] + [7 standard ranks] + [any subspecies]
    levels = split_string[:2] + levels
    if s_found is not None:
        # Append any subspecies or strain-level annotations that come after species
        levels = levels + split_string[s_found + 1 :]

    # Join all non-empty levels with semicolons
    tax_list = ";".join([x for x in levels if x != ""])
    return tax_list


def format_tax_kalmari(string):
    """
    Filter out non-standard intermediate taxonomy entries from Kalmari database output.

    Removes entries prefixed with "-_" (which indicate intermediate or non-standard ranks)
    while preserving:
    - The first two entries (root and domain, even if they start with "-_")
    - Any entries after species rank (subspecies/strains)
    - All standard rank entries (k_, p_, c_, o_, f_, g_, s_)
    """
    formatted_tax = []
    last_tax = ""
    for i, tax in enumerate(string.split(";")):
        # Keep taxonomy entries that:
        # 1. Don't start with "-_" (standard ranks like k_, p_, etc.)
        # 2. Are in positions 0 or 1 (root/domain level, always keep)
        # 3. Come after species (subspecies annotations)
        if (
            not tax.startswith(("-")) or i in [0, 1] or last_tax.startswith("s_")
        ):  # Filter out all taxnomy not all entries have. These start with -. But keep subspecies (which is last) and bacteria (which start with -)
            formatted_tax.append(tax)
        last_tax = tax
    formatted_tax = ";".join(formatted_tax)
    return formatted_tax


def test_all_levels_present(string):
    """
    Validate that the formatted taxonomy string has the correct number of levels.

    After formatting, the taxonomy should have a specific number of levels depending on
    the deepest rank present:
    - If deepest is species (s_): should have 7 levels (k,p,c,o,f,g,s)
    - If deepest is genus (g_): should have 6 levels (k,p,c,o,f,g)
    - And so on up the hierarchy

    Raises ValueError if validation fails (incorrect number of levels for the deepest rank).
    """
    formated = format_tax_kalmari(string)

    # Define expected number of levels for each deepest rank
    how_many_levels_should_be = {
        "s": 7,  # Species: k,p,c,o,f,g,s
        "g": 6,  # Genus: k,p,c,o,f,g
        "f": 5,  # Family: k,p,c,o,f
        "o": 4,  # Order: k,p,c,o
        "c": 3,  # Class: k,p,c
        "p": 2,  # Phylum: k,p
        "k": 1,  # Kingdom: k
    }

    # Count taxonomy levels by iterating backwards from deepest to shallowest
    start = False
    count = 0
    first_tax = None  # Will hold the deepest (most specific) taxonomic rank found
    for tax in reversed(formated.split(";")):
        if tax.startswith("-_"):
            start = False  # Skip non-standard entries (shouldn't happen after format_tax_kalmari)
        if not tax.startswith("-_"):
            start = True
            if first_tax is None:
                first_tax = tax  # Capture the deepest rank (found first when going backwards)
        if start:
            count += 1

    # Allow special cases: high-level classifications that are intentionally incomplete
    # (e.g., only classified to domain level, not to species)
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
        pass  # These are valid incomplete taxonomies
    # Validate: does the number of levels match what we expect for this deepest rank?
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
    """
    Main function: Process a TSV file of contig names and taxonomies.

    For each line:
    1. Parse contig name and taxonomy string (tab-separated)
    2. If taxonomy is empty, pass through as-is
    3. Otherwise, fix missing taxonomic levels and validate the result
    4. Output formatted tab-separated contig\ttaxonomy to stdout
    """
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
            # Parse tab-separated contig name and taxonomy
            name, tax = line.split("\t")
            tax = tax.strip()
            if tax == "":
                # Empty taxonomy: pass through as empty
                tax = ""
            else:
                # Fix and validate the taxonomy
                tax = fix_tax(tax)
                test_all_levels_present(tax)

                # print(tax)
            # Output formatted result to stdout
            print(f"{name}\t{tax}")
            # print(name)

## Test example based on a real taxonomy string from the pipeline
def test_fix_tax():
    """
    Test the fix_tax function with a realistic example.

    Input taxonomy has:
    - Root and domain entries: "-_cellular organisms;-_Bacteria"
    - Kingdom: "k_Pseudomonadati"
    - Phylum: "p_Acidobacteriota"
    - Intermediate non-standard entry: "-_unclassified Acidobacteriota"
    - Species: "s_Acidobacteriota bacterium"
    - Missing: class, order, family, genus

    Expected behavior:
    - Preserve root, domain, kingdom, phylum, species
    - Fill in missing class, order, family, genus with placeholder names
    - Remove or handle the "-_unclassified Acidobacteriota" entry appropriately
    """
    test_string = "-_cellular organisms;-_Bacteria;k_Pseudomonadati;p_Acidobacteriota;-_unclassified Acidobacteriota;s_Acidobacteriota bacterium"

    result = fix_tax(test_string)

    print("Test fix_tax():")
    print(f"Input:  {test_string}")
    print(f"Output: {result}")
    print()

    # Verify the result has the expected structure
    parts = result.split(";")

    # Should have root, domain, k, p, c (filled), o (filled), f (filled), g (filled), s
    assert len(parts) >= 9, f"Expected at least 9 parts, got {len(parts)}: {parts}"

    # Check that standard ranks are present
    has_k = any("k_" in p for p in parts)
    has_p = any("p_" in p for p in parts)
    has_s = any("s_" in p for p in parts)

    assert has_k, f"Expected kingdom (k_) in result: {result}"
    assert has_p, f"Expected phylum (p_) in result: {result}"
    assert has_s, f"Expected species (s_) in result: {result}"

    # Check that missing levels were filled with placeholders
    has_placeholder = any("LEVEL_" in p for p in parts)
    assert has_placeholder, f"Expected placeholder levels to be added: {result}"

    print("Test passed: fix_tax() correctly fills missing taxonomic levels\n")

    return result


if __name__ == "__main__":
    # Run test if script is executed directly (not when called from Snakemake)
    import sys
    if len(sys.argv) == 1:
        # No arguments: run test
        test_fix_tax()
    else:
        # Arguments provided: run main function (normal pipeline usage)
        typer.run(main)


