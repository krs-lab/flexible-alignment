import argparse
import os
from Bio import PDB


# ==============================================================================
# SECTION 1: CUSTOM PDB RESIDUE SELECTOR CLASS
# ==============================================================================
class ResidueSelect(PDB.Select):
  """Custom selection rule passed to Bio.PDB.PDBIO to filter out residues outside

  a specific numeric range (start_res to end_res).
  """

  def __init__(self, start_res, end_res):
    self.start_res = start_res
    self.end_res = end_res

  def accept_residue(self, residue):
    # Retrieve residue sequence identifier (ignores insertion codes)
    res_num = residue.get_id()[1]
    # Return True only if residue falls within target window
    return self.start_res <= res_num <= self.end_res


# ==============================================================================
# SECTION 2: CORE PDB SLICING FUNCTION
# ==============================================================================
def slice_pdb(pdb_path, output_dir, segment_len=200, min_tail_len=87):
  """Reads a single PDB file, extracts residue ranges in steps of segment_len,

  and writes out individual sliced PDB files.
  """
  # Initialize BioPython PDB parser and file writer
  parser = PDB.PDBParser(QUIET=True)
  structure = parser.get_structure("protein", pdb_path)
  writer = PDB.PDBIO()

  # Extract baseline filename without extension (e.g., 'AF-P12345-F1' from 'AF-P12345-F1.pdb')
  base_name = os.path.splitext(os.path.basename(pdb_path))[0]

  # Iterate over models and polypeptide chains in the structure
  for model in structure:
    for chain in model:
      # Filter for standard amino acid residues to avoid empty/non-protein records
      residues = [r for r in chain if PDB.is_aa(r)]
      if not residues:
        continue

      # Determine start and end residue sequence numbers
      start_id = residues[0].get_id()[1]
      end_id = residues[-1].get_id()[1]

      seg_num = 1
      curr_start = start_id

      # Loop through the amino acid sequence step-by-step
      while curr_start <= end_id:
        curr_end = curr_start + segment_len - 1
        remaining = end_id - curr_start + 1

        # ----------------------------------------------------------------------
        # LOGIC: TERMINAL TAIL DISCARD THRESHOLD
        # If remaining length < 87aa AND this is not the first segment (seg_num > 1),
        # skip processing and exit loop. Single short proteins (< 200aa total, seg_num == 1)
        # bypass this check and are saved as-is.
        # ----------------------------------------------------------------------
        if remaining < min_tail_len and seg_num > 1:
          break  # Skip leftover tails smaller than 87aa

        # Build output path (e.g., output_dir/P12345_seg1_1_200.pdb)
        out_filename = os.path.join(
            output_dir, f"{base_name}_seg{seg_num}_{curr_start}_{curr_end}.pdb"
        )

        # Apply ResidueSelect filter and write segment to file
        writer.set_structure(structure)
        writer.save(out_filename, ResidueSelect(curr_start, curr_end))

        # Advance window pointer to start of next segment
        curr_start += segment_len
        seg_num += 1


# ==============================================================================
# SECTION 3: CLI ARGUMENT PARSING & WORKFLOW EXECUTION
# ==============================================================================
def main():
  # Configure ArgumentParser documentation and --help flags
  parser = argparse.ArgumentParser(
      description=(
          "Slices AlphaFold/PDB structure files into fixed-length amino acid"
          " segments (default: 200aa) for structural homology search"
          " comparisons."
      ),
      formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog="""
Usage Examples:
  # Basic usage with default 200aa slicing and 87aa tail discard limit:
  python3 fixed-length_slicer.py -i ./pdb_dir/

  # Custom segment length and tail threshold:
  python3 fixed-length_slicer.py -i ./pdb_dir/ -o ./custom_output/ -s 150 -m 87
      """,
  )

  # Define flags
  parser.add_argument(
      "-i",
      "--input",
      required=True,
      help="Path to directory containing input .pdb files.",
  )
  parser.add_argument(
      "-o",
      "--output",
      default=None,
      help=(
          "Output directory for sliced PDBs (Default: create a 'sliced_<len>aa'"
          " directory inside the input directory)."
      ),
  )
  parser.add_argument(
      "-s",
      "--segment-length",
      type=int,
      default=200,
      help="Amino acid length per slice (Default: 200).",
  )
  parser.add_argument(
      "-m",
      "--min-tail",
      type=int,
      default=87,
      help=(
          "Minimum residue length threshold to keep a leftover terminal tail"
          " (Default: 87)."
      ),
  )

  args = parser.parse_args()

  # Path Normalization: Standardize input path string and remove trailing slashes
  input_dir = os.path.abspath(os.path.normpath(args.input))

  # Input Validation: Ensure specified input folder exists
  if not os.path.exists(input_dir):
    print(f"Error: Input directory '{input_dir}' does not exist.")
    return

  # Output Setup: Assign custom path or create default subfolder inside input directory
  if args.output:
    output_dir = os.path.abspath(os.path.normpath(args.output))
  else:
    output_dir = os.path.join(input_dir, f"sliced_{args.segment_length}aa")

  os.makedirs(output_dir, exist_ok=True)

  # Target Discovery: Collect all files ending in .pdb
  pdb_files = [f for f in os.listdir(input_dir) if f.endswith(".pdb")]

  if not pdb_files:
    print(f"No .pdb files found in '{input_dir}'.")
    return

  # Print execution configuration parameters to terminal
  print("=" * 60)
  print("PDB SEGMENTATION RUN CONFIGURATION")
  print("=" * 60)
  print(f"Input Directory  : {input_dir}")
  print(f"Output Directory : {output_dir}")
  print(f"Segment Length   : {args.segment_length} aa")
  print(f"Min Tail Length  : {args.min_tail} aa")
  print(f"Found PDB Files  : {len(pdb_files)}")
  print("=" * 60)

  # Execution Loop: Process every discovered PDB file sequentially
  for pdb_file in pdb_files:
    pdb_path = os.path.join(input_dir, pdb_file)
    slice_pdb(
        pdb_path,
        output_dir,
        segment_len=args.segment_length,
        min_tail_len=args.min_tail,
    )

  print(f"\nProcessing complete! Sliced PDBs saved to: {output_dir}")


if __name__ == "__main__":
  main()
