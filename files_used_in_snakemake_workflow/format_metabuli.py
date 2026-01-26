"""
Format Metabuli Taxonomy Output for Taxvamb

This script converts Metabuli taxonomy classification output into a format compatible with Taxvamb.

Purpose:
--------
Metabuli produces two output files during taxonomic classification:
1. A classification file (filepath_clas) containing contig-to-taxID assignments
2. A report file (filepath_report) containing the hierarchical taxonomy structure with taxIDs

This script bridges these files by:
- Parsing the hierarchical taxonomy structure from the report file
- Creating a mapping from taxonomic IDs to full taxonomy lineages
- Converting the classification file to use full lineage strings instead of taxIDs
- Formatting the output to match Taxvamb's expected input format (tab-separated: contigs\tpredictions)

Pipeline Integration:
--------------------
This script is called in the Snakemake workflow during the Metabuli-based Taxvamb runs.
It processes the raw Metabuli output before passing it to the Taxvamb binning algorithm.

Input:
------
- filepath_clas: Metabuli classification file (tab-separated, with taxID in column 3, contig name in column 2)
- filepath_report: Metabuli report file (tab-separated, with taxID in column 5, taxonomy name in column 6)

Output:
-------
Prints to stdout a tab-separated file with header:
- Column 1: Contig name
- Column 2: Full taxonomy lineage (semicolon-separated), or empty string for unclassified contigs

Note: The first 5 characters of each taxonomy string are stripped (removes prefix characters).
"""

# Script to ensure that the format of the metabulti database outputs are consistent with the fomrat taxvamb expect

from pathlib import Path

import typer


def main(filepath_clas: Path, filepath_report: Path):
    """
    Convert Metabuli classification output to Taxvamb-compatible format.

    Args:
        filepath_clas: Path to Metabuli classification file (contig assignments)
        filepath_report: Path to Metabuli report file (taxonomy hierarchy)
    """
    # ===== PHASE 1: Build tax_id -> taxonomy lineage mapping from report file =====
    # Create mapping of tax_id to taxonomy based on the metabuli_report
    curr_taxonomy = []  # Tracks current taxonomy path as we traverse the hierarchy
    tax_id2taxonomy = {}  # Final mapping: tax_id -> full semicolon-separated lineage

    with open(filepath_report) as f:
        for line in f:
            # Extract taxonomic ID (column 5, index 4)
            id_line = line.split("\t")[4].strip()

            # Extract taxonomy name (column 6, index 5)
            # The name is indented with spaces to show hierarchy depth
            line = line.split("\t")[5].rstrip()

            # Determine hierarchy level by counting leading spaces
            # More spaces = deeper in the taxonomy tree
            curr_space = len(line.split(" ")) - 1

            # Expand curr_taxonomy array if we're going deeper than before
            while len(curr_taxonomy) - 1 < curr_space:
                curr_taxonomy.append("")

            # Store this taxon at its hierarchy level (overwrites previous at same level)
            curr_taxonomy[curr_space] = line

            # Build the full lineage from root to current level
            to_add = []
            for i, entry in enumerate(curr_taxonomy):
                if i > curr_space:
                    break  # Don't include levels deeper than current
                to_add.append(entry)

            # Create semicolon-separated lineage string (e.g., "Bacteria;Proteobacteria;...")
            tax_id2taxonomy[id_line] = ";".join([x.strip() for x in to_add if x != ""])

    # ===== PHASE 2: Convert classification file using the taxonomy mapping =====
    # Convert the metabuli_classifications to taxonomy using the tax_id2taxonomy mapping
    print("contigs\tpredictions")  # Header line for Taxvamb input format

    with open(filepath_clas, "r") as f:
        for line in f:
            sp_line = line.split("\t")

            # Column 2 (index 1): contig name
            # Column 3 (index 2): taxonomic ID assigned to this contig
            search_id = sp_line[2]

            # Look up the full taxonomy lineage for this tax_id
            tax = tax_id2taxonomy[search_id]

            # Handle unclassified contigs: convert to empty string
            if tax == "unclassified":
                tax = ""

            # Output: contig_name \t taxonomy_lineage
            # Note: tax[5:] strips the first 5 characters (likely removes a prefix)
            print(f"{sp_line[1]}\t{tax[5:]}")


## Test example based on realistic Metabuli output
def test_format_metabuli():
    """
    Test the Metabuli formatting with synthetic but realistic data.

    Creates temporary test files simulating Metabuli output and verifies
    that the script correctly maps tax_ids to full lineages.
    """
    import tempfile
    import os

    # Create temporary test files
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a sample Metabuli report file
        # Format: tab-separated with tax_id in column 5 (index 4) and taxonomy in column 6 (index 5)
        # Indentation (spaces) in column 6 indicates hierarchy depth
        report_path = os.path.join(tmpdir, "test_report.txt")
        with open(report_path, "w") as f:
            # Root level (0 spaces)
            f.write("col1\tcol2\tcol3\tcol4\t1\troot\n")
            # Domain level (1 space)
            f.write("col1\tcol2\tcol3\tcol4\t2\t Bacteria\n")
            # Phylum level (2 spaces)
            f.write("col1\tcol2\tcol3\tcol4\t3\t  Proteobacteria\n")
            # Class level (3 spaces)
            f.write("col1\tcol2\tcol3\tcol4\t4\t   Gammaproteobacteria\n")
            # Order level (4 spaces)
            f.write("col1\tcol2\tcol3\tcol4\t5\t    Enterobacterales\n")
            # Family level (5 spaces)
            f.write("col1\tcol2\tcol3\tcol4\t6\t     Enterobacteriaceae\n")
            # Genus level (6 spaces)
            f.write("col1\tcol2\tcol3\tcol4\t7\t      Escherichia\n")
            # Species level (7 spaces)
            f.write("col1\tcol2\tcol3\tcol4\t8\t       Escherichia coli\n")
            # Unclassified entry
            f.write("col1\tcol2\tcol3\tcol4\t0\tunclassified\n")

        # Create a sample Metabuli classification file
        # Format: tab-separated with contig_name in column 2 (index 1) and tax_id in column 3 (index 2)
        clas_path = os.path.join(tmpdir, "test_classification.txt")
        with open(clas_path, "w") as f:
            # Contig assigned to E. coli (tax_id 8)
            f.write("col1\tcontig_001\t8\tother_cols\n")
            # Contig assigned to Enterobacteriaceae family (tax_id 6)
            f.write("col1\tcontig_002\t6\tother_cols\n")
            # Contig assigned to Bacteria domain (tax_id 2)
            f.write("col1\tcontig_003\t2\tother_cols\n")
            # Unclassified contig (tax_id 0)
            f.write("col1\tcontig_004\t0\tother_cols\n")

        # Run the main function and capture output
        print("Test format_metabuli():")
        print(f"Report file: {report_path}")
        print(f"Classification file: {clas_path}")
        print("\nOutput:")
        print("-" * 60)

        # Redirect stdout to capture output
        import sys
        from io import StringIO
        old_stdout = sys.stdout
        sys.stdout = captured_output = StringIO()

        # Run the function
        main(Path(clas_path), Path(report_path))

        # Restore stdout
        sys.stdout = old_stdout
        output = captured_output.getvalue()

        print(output)
        print("-" * 60)

        # Validate output
        lines = output.strip().split("\n")

        # Check header
        assert lines[0] == "contigs\tpredictions", f"Expected header 'contigs\\tpredictions', got '{lines[0]}'"

        # Check contig_001 -> E. coli (full lineage)
        # Expected: root;Bacteria;Proteobacteria;...;Escherichia coli (with first 5 chars removed)
        assert "contig_001" in lines[1], f"Expected contig_001 in line 1, got '{lines[1]}'"
        assert "Escherichia coli" in lines[1], f"Expected 'Escherichia coli' in line 1, got '{lines[1]}'"

        # Check contig_002 -> Enterobacteriaceae (partial lineage)
        assert "contig_002" in lines[2], f"Expected contig_002 in line 2, got '{lines[2]}'"
        assert "Enterobacteriaceae" in lines[2], f"Expected 'Enterobacteriaceae' in line 2, got '{lines[2]}'"

        # Check contig_003 -> Bacteria (domain only)
        assert "contig_003" in lines[3], f"Expected contig_003 in line 3, got '{lines[3]}'"

        # Check contig_004 -> unclassified (should be empty)
        assert "contig_004" in lines[4], f"Expected contig_004 in line 4, got '{lines[4]}'"
        # After tab, should be empty or very short (just the stripped portion)
        parts = lines[4].split("\t")
        assert len(parts) == 2, f"Expected 2 columns, got {len(parts)}"
        # The taxonomy part should be empty for unclassified
        # (Note: might be empty string or might have been set to empty in the function)

        print("\n✓ All tests passed: format_metabuli() correctly converts tax_ids to lineages\n")


if __name__ == "__main__":
    # Run test if script is executed directly (not when called from Snakemake)
    import sys
    if len(sys.argv) == 1:
        # No arguments: run test
        test_format_metabuli()
    else:
        # Arguments provided: run main function (normal pipeline usage)
        typer.run(main)
