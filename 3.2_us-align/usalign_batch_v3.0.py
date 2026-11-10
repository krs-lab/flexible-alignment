import argparse
import csv
from datetime import datetime
import os
from pathlib import Path
import re
import subprocess

# Define the absolute path to your US-align binary (avoids write privilege issues)
USALIGN_BIN = "/fast_disk/visnu/USalign/USalign"


def parse_usalign_output(stdout):
  """Extracts TM-scores, RMSD, Aligned Length, and Sequence Identity from US-align stdout using regex."""
  tm_scores = re.findall(r"TM-score=\s*([0-9.]+)", stdout)
  rmsd_match = re.search(r"RMSD=\s*([0-9.]+)", stdout)
  aln_match = re.search(r"Aligned length=\s*([0-9]+)", stdout)
  seq_match = re.search(r"Seq_ID[^=\s]*=\s*([0-9.]+)", stdout)

  tm1 = tm_scores[0] if len(tm_scores) > 0 else ""
  tm2 = tm_scores[1] if len(tm_scores) > 1 else ""
  rmsd = rmsd_match.group(1) if rmsd_match else ""
  aln_len = aln_match.group(1) if aln_match else ""
  seq_id = seq_match.group(1) if seq_match else ""

  return tm1, tm2, rmsd, aln_len, seq_id


def run_batch_alignments(query_dir, target_dir):
  # 1. Create timestamped results directory structure
  timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
  results_dir = Path(f"results_{timestamp}")
  results_dir.mkdir(parents=True, exist_ok=True)

  sup_dir = results_dir / "superimposed_pdbs"
  sup_dir.mkdir(exist_ok=True)

  csv_path = results_dir / "alignment_metrics_summary.csv"

  # 2. Find all PDB files in both directories
  q_path_dir = Path(query_dir)
  t_path_dir = Path(target_dir)

  queries = sorted(
      list(q_path_dir.glob("*.pdb")) + list(q_path_dir.glob("*.ent"))
  )
  targets = sorted(
      list(t_path_dir.glob("*.pdb")) + list(t_path_dir.glob("*.ent"))
  )

  if not queries:
    print(f"[Error] No PDB files found in query directory: {query_dir}")
    return
  if not targets:
    print(f"[Error] No PDB files found in target directory: {target_dir}")
    return

  print(f"[Info] Found {len(queries)} query(s) and {len(targets)} target(s).")
  print(f"[Info] Output directory created: {results_dir.resolve()}\n")

  fieldnames = [
      "Query",
      "Target",
      "Aligned_Length",
      "RMSD",
      "Seq_ID",
      "TM_Score_Query",
      "TM_Score_Target",
  ]

  results_rows = []
  total_pairs = len(queries) * len(targets)
  current = 0

  # 3. Nested loop: Each query vs each target
  for q_file in queries:
    for t_file in targets:
      current += 1
      q_name = q_file.stem
      t_name = t_file.stem

      # Define prefix path for US-align output files (-o option creates sup.pdb)
      sup_prefix = sup_dir / f"{q_name}_vs_{t_name}"
      cmd = [USALIGN_BIN, str(q_file), str(t_file), "-o", str(sup_prefix)]

      print(f"[{current}/{total_pairs}] Aligning {q_name} against {t_name}...")

      try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
        tm1, tm2, rmsd, aln_len, seq_id = parse_usalign_output(result.stdout)

        results_rows.append({
            "Query": q_name,
            "Target": t_name,
            "Aligned_Length": aln_len,
            "RMSD": rmsd,
            "Seq_ID": seq_id,
            "TM_Score_Query": tm1,
            "TM_Score_Target": tm2,
        })

      except subprocess.CalledProcessError as e:
        print(f"  [Error] US-align failed for {q_name} vs {t_name}: {e}")
      except Exception as e:
        print(f"  [Error] Unexpected processing error: {e}")

  # 4. Sort results rows by TM_Score_Query in descending order (highest score first)
  results_rows.sort(
      key=lambda x: float(x["TM_Score_Query"])
      if x["TM_Score_Query"] != ""
      else 0.0,
      reverse=True,
  )

  # 5. Write sorted rows out to main summary CSV
  with open(csv_path, mode="w", newline="") as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results_rows)

  # 6. Output individual CSV files for each query structure
  query_rows_dict = {}
  for row in results_rows:
    q = row["Query"]
    if q not in query_rows_dict:
      query_rows_dict[q] = []
    query_rows_dict[q].append(row)

  indiv_dir = results_dir / "individual_query_csvs"
  indiv_dir.mkdir(exist_ok=True)

  for q_name, q_rows in query_rows_dict.items():
    indiv_csv_path = indiv_dir / f"{q_name}_alignment_metrics.csv"
    with open(indiv_csv_path, mode="w", newline="") as indiv_file:
      writer = csv.DictWriter(indiv_file, fieldnames=fieldnames)
      writer.writeheader()
      writer.writerows(q_rows)

  print(f"\n[Success] Batch processing complete!")
  print(f"-> Summary CSV (Sorted by TM_Score_Query desc): {csv_path.resolve()}")
  print(f"-> Individual Query CSVs: {indiv_dir.resolve()}")
  print(f"-> Superimposed PDBs: {sup_dir.resolve()}")


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
      description=(
          "Batch structural alignment tool using US-align to compare a"
          " directory of query PDBs against a directory of target PDBs."
      ),
      epilog=(
          "Examples:\n"
          "  python3 script.py -q query_dir -t target_dir\n"
          "  python3 script.py --quary query_dir --target target_dir\n"
      ),
      formatter_class=argparse.RawDescriptionHelpFormatter,
  )
  parser.add_argument(
      "-q",
      "--query",
      "--quary",
      default=None,
      help="Path to the directory containing query PDB files.",
  )
  parser.add_argument(
      "-t",
      "--target",
      default=None,
      help="Path to the directory containing target PDB files.",
  )

  args = parser.parse_args()

  # Handles standard '--query' as well as user's '--quary' spelling
  q_dir = args.query
  t_dir = args.target

  # Fallback to interactive prompts if flags weren't provided
  if not q_dir:
    q_dir = input("Enter path to query PDB directory: ").strip()
  if not t_dir:
    t_dir = input("Enter path to target PDB directory: ").strip()

  run_batch_alignments(q_dir, t_dir)