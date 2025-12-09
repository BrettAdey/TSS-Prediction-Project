#!/usr/bin/env bash
#SBATCH --job-name=build_T2T_from100kb --output=logs/build_%A_%a.out --error=logs/build_%A_%a.err --time=6:00:00 --cpus-per-task=4 --mem=128G --array=1-24 
set -euo pipefail

### Set up ###
input_dir="./T2T_fastas_by_chr"

biasaway="$HOME/bin/biasaway"
python="$HOME/bin/python"
seqkit="$HOME/bin/seqkit"

seed=42
window=100
step=50
wrap=100000

holdout_chrs=("NC_060932.1" "NC_060933.1" "NC_060934.1")
manifest="T2T_fastas_by_chr_manifest.txt"

# Make output dirs
mkdir -p \
  T2T_forward \
  T2T_reversed \
  dinucleotide_shuffled/global \
  dinucleotide_shuffled/local \
  mononucleotide_shuffled/global \              
  mononucleotide_shuffled/local \
  repeat_mononucleotide_global_shuffled \
  nonrepeat_mononucleotide_global_shuffled \
  repeat_reversed \
  nonrepeat_reversed \
  holdout/forward \
  holdout/reversed \
  logs tmp

# Read each line into CHRS ( -t strips trailing newlines)
mapfile -t CHRS < "$manifest"

# Get number of items in CHRS
N=${#CHRS[@]}

# Read current SLURM array index
task_id=${SLURM_ARRAY_TASK_ID:-1}

# Select the chromosome fasta from the manifest corrosponding to the current array index
chr_fasta=${CHRS[$((task_id-1))]}

# Get filename stem
base=$(basename "$chr_fasta" .fasta)

# Report
echo "Task $task_id -> $base"
echo "Input: $chr_fasta"

### Functions ###
strip_ba_headers() {
  # Converts "background_seq_1 Background sequence for ID|a|b ..." to "ID|a|b"
  local in="$1"; local out="$2"
  "$seqkit" replace \
    -p '^[^ ]+\s+Background\s+sequence\s+for\s+(\S+).*' \
    -r '${1}' \
    "$in" > "$out"
}

# Append a label to the header and rewrap to single-line 100kb fasta
reheader_and_wrap() {
  local label="$1"; local in="$2"; local out="$3"
  "$seqkit" replace -p '^(.+)$' -r '${1}_'"$label" "$in" \
  | "$seqkit" seq -w "$wrap" > "$out"
}

# Reverse only lowercase (repeat) runs, or only uppercase (non-repeat) runs
py_rev_lower() {
$python - "$1" <<'python'
import sys,re
from Bio import SeqIO
from Bio.Seq import Seq

# Parse fastas into Biopython SeqRecord
# Regex pattern selects for any contiguous string of upper or lowercase characters
fasta=sys.argv[1]; pat=re.compile(r'[A-Z]+|[a-z]+')
for rec in SeqIO.parse(fasta,"fasta"):
  # Convert the sequence in to a string and save it as variable 's'
  # Initialise empty list names 'out'
  # Note rec.seq is the sequence info only (separated thanks to Biopython), so no header issues
    s=str(rec.seq); out=[]
  # Finditer returns an iterator of match objects - one for each substring of s which matches the pattern.
    for m in pat.finditer(s):
  # Group(0) selects the sequence itself from the match object (and this saves it to 'seg')
        seg=m.group(0)
  # Reverse and append to 'out' list if it's lower, otherwise append the original
        out.append(seg[::-1] if seg.islower() else seg)
  # Join the out list into one long string separated by ""
  # Turn that into a Biopython seq object
  # Save this as the sequence of the current record (i.e. overwrite the input rec.seq)
  # Write the output as fasta file to (standard out) - or in this case, the output_buffer to view it
    rec.seq=Seq("".join(out)); SeqIO.write(rec, sys.stdout, "fasta")
python
}
py_rev_upper() {
$python - "$1" <<'python'
import sys,re
from Bio import SeqIO
from Bio.Seq import Seq
fasta=sys.argv[1]; pat=re.compile(r'[A-Z]+|[a-z]+')

for rec in SeqIO.parse(fasta,"fasta"):
    s=str(rec.seq); out=[]
    for m in pat.finditer(s):
        seg=m.group(0)
        out.append(seg[::-1] if seg.isupper() else seg)
    rec.seq=Seq("".join(out)); SeqIO.write(rec, sys.stdout, "fasta")
python
}

# Extract lowercase runs (repeats) as separate fasta records with start/stop coordinates
py_extract_lower() {
$python - "$1" <<'python'
import sys, re
from Bio import SeqIO
pat = re.compile(r'[a-z]+')
for rec in SeqIO.parse(sys.argv[1], "fasta"):
    s = str(rec.seq)
    for m in pat.finditer(s):
        # Get start and end coordinates for this match
        a, b = m.start(), m.end()
        # Print header + match sequence ('|' as the delimiter, id|start|stop)
        print(f">{rec.id}|{a}|{b}\n{m.group(0)}")
python
}

# Extract uppercase runs (non-repeats) as separate fasta records with start/stop coordinates
py_extract_upper() {
$python - "$1" <<'python'
import sys, re
from Bio import SeqIO
pat = re.compile(r'[A-Z]+')
for rec in SeqIO.parse(sys.argv[1], "fasta"):
    s = str(rec.seq)
    for m in pat.finditer(s):
        a, b = m.start(), m.end()
        print(f">{rec.id}|{a}|{b}\n{m.group(0)}")
python
}

# Reinsert shuffled runs back into the original mixed-case chromosome
py_reinsert_runs_ordered() {
$python - "$1" "$2" "$3" <<'python'
import sys, collections, re
from Bio import SeqIO
from Bio.Seq import Seq

fasta, target, shuf_fasta = sys.argv[1], sys.argv[2].lower(), sys.argv[3]

### Load shuffled fragments grouped by window ID, preserving order ###
# Using python's deque (double-ended queue, 'deck') object. Allows you to pop in and out sequences in order
# By_win: each row is a 100kb sequence and a deque of shuffled sequences, deque(['ACTTTAGGAAC','ATTTGGACCTTAAG']) 
by_win = collections.defaultdict(collections.deque)
loaded = 0

# Get ID Start and End from each header
# Returns a Tuple with three elements
# Description is where the header is stored in the rec object
def parse_key(rec):
  # Split by white space and select just the first token to trim off any extra header description
    tok = rec.description.split()[0]
  # Then split by '|', limiting to three parts (ID, Start, End)
    parts = tok.split("|", 2)
    if len(parts) < 3:
        return None
    return parts[0], parts[1], parts[2]

for rec in SeqIO.parse(shuf_fasta, "fasta"):
  # Save 3 element tuple as 'key'
    key = parse_key(rec)
    if key is None:
        continue
    # Convert this record's sequence to a string
    # Append it to by_win, using the ID part of the key - note it's being appended in the same order
    by_win[key[0]].append(str(rec.seq))
    # Increment Counter
    loaded += 1

sys.stderr.write(f"[reinsert-ordered] loaded={loaded} windows={len(by_win)}\n")

### Using the new by_win object to reinsert shuffled sequences back in place

seg_re = re.compile(r'[A-Z]+|[a-z]+')
replaced = 0
skipped  = 0
for rec in SeqIO.parse(fasta, "fasta"):
    winid = rec.id
  # Lookup this record ID in the by_win dictionary
  # If it exists, set dq to it
  # If not, set it to a new deque object
    dq = by_win.get(winid, collections.deque())
  # Convert this record's sequence to a literal string
    s  = str(rec.seq)
    out = []
    for m in seg_re.finditer(s):
        seg = m.group(0)
    # If we're looking for lowercase, and the segment is lower, and dq = true (not empty)/the record is a match, i.e. a sequence we want to replace:
        if target == "lower" and seg.islower() and dq:
    # Replace by popping out the first element in dq
            repl = dq.popleft()
    # Then appending that, lowercased, to the output (and incrementing the counter)
            out.append(repl.lower()); replaced += 1
    # Same way for uppercase non-repeats
        elif target == "upper" and seg.isupper() and dq:
            repl = dq.popleft()
            out.append(repl.upper()); replaced += 1
    # If the matched segment is the opposite case to what we want to shuffle, or there isn't a match, we just append the segment as-is
        else:
            out.append(seg); skipped += 1
    # As before, replace rec.seq with concatenated 'out' string, converted to Seq object
    rec.seq = Seq(''.join(out))
    SeqIO.write(rec, sys.stdout, "fasta")

left = sum(len(v) for v in by_win.values())
sys.stderr.write(f"[reinsert-ordered] replaced={replaced} skipped={skipped} leftover_fragments={left}\n")
python
}

### Create T2T Variants ###

# Forward
reheader_and_wrap "forward" "$chr_fasta" "T2T_forward/${base}.fasta"

# Reversed
"$seqkit" seq -r "$chr_fasta" > "tmp/${base}.rev.fa"
reheader_and_wrap "reversed" "tmp/${base}.rev.fa" "T2T_reversed/${base}.reversed.fasta"
rm -f "tmp/${base}.rev.fa"

# Dinucleotide shuffles
# global k=2
"$biasaway" k -f "$chr_fasta" -k 2 --seed "$seed" 1>"tmp/${base}.gdi.raw.fa" 2>"logs/${base}.biasaway_gdi.log"
strip_ba_headers "tmp/${base}.gdi.raw.fa" "tmp/${base}.gdi.fa"
reheader_and_wrap "globally_dinucleotide_shuffled" "tmp/${base}.gdi.fa" \
  "dinucleotide_shuffled/global/${base}.global_dinuc_shuffled.fasta"
rm -f "tmp/${base}.gdi.raw.fa" "tmp/${base}.gdi.fa"

# local w k=2
"$biasaway" w -f "$chr_fasta" -k 2 -w "$window" -s "$step" --seed "$seed" 1>"tmp/${base}.ldi.raw.fa" 2>"logs/${base}.biasaway_ldi.log"
strip_ba_headers "tmp/${base}.ldi.raw.fa" "tmp/${base}.ldi.fa"
reheader_and_wrap "locally_dinucleotide_shuffled" "tmp/${base}.ldi.fa" \
  "dinucleotide_shuffled/local/${base}.local_dinuc_shuffled.fasta"
rm -f "tmp/${base}.ldi.raw.fa" "tmp/${base}.ldi.fa"

# Mononucleotide shuffles
# global k=1
"$biasaway" k -f "$chr_fasta" -k 1 --seed "$seed" 1>"tmp/${base}.gmono.raw.fa" 2>"logs/${base}.biasaway_gmono.log"
strip_ba_headers "tmp/${base}.gmono.raw.fa" "tmp/${base}.gmono.fa"
reheader_and_wrap "globally_mononucleotide_shuffled" "tmp/${base}.gmono.fa" \
  "mononucleotide_shuffled/global/${base}.global_mononuc_shuffled.fasta"
rm -f "tmp/${base}.gmono.raw.fa" "tmp/${base}.gmono.fa"

# local w k=1
"$biasaway" w -f "$chr_fasta" -k 1 -w "$window" -s "$step" --seed "$seed" 1>"tmp/${base}.lmono.raw.fa" 2>"logs/${base}.biasaway_lmono.log"
strip_ba_headers "tmp/${base}.lmono.raw.fa" "tmp/${base}.lmono.fa"
reheader_and_wrap "locally_mononucleotide_shuffled" "tmp/${base}.lmono.fa" \
  "mononucleotide_shuffled/local/${base}.local_mononuc_shuffled.fasta"
rm -f "tmp/${base}.lmono.raw.fa" "tmp/${base}.lmono.fa"

# Repeat-only / Nonrepeat-only mono shuffle (global k=1)
# Extract lowercase (repeats), shuffle, strip headers, reinsert, window + label
py_extract_lower "$chr_fasta" > "tmp/${base}_rep.fa"
if [[ -s "tmp/${base}_rep.fa" ]]; then
  "$biasaway" k -f "tmp/${base}_rep.fa" -k 1 --seed "$seed" 1>"tmp/${base}.rep_mono.raw.fa" 2>"logs/${base}.biasaway_repmono.log"
  strip_ba_headers "tmp/${base}.rep_mono.raw.fa" "tmp/${base}.rep_mono.shuf.fa"
  py_reinsert_runs_ordered "$chr_fasta" lower "tmp/${base}.rep_mono.shuf.fa" \
  > "tmp/${base}.rep_mono.reins.fa" 2> "logs/${base}.reinsertion_repeats.log"
else
  cp "$chr_fasta" "tmp/${base}.rep_mono.reins.fa"
fi
reheader_and_wrap "globally_mononucleotide_shuffled_repeats" \
  "tmp/${base}.rep_mono.reins.fa" \
  "repeat_mononucleotide_global_shuffled/${base}.global_mononuc_shuffled_repeats.fasta"
rm -f "tmp/${base}_rep.fa" "tmp/${base}.rep_mono.raw.fa" "tmp/${base}.rep_mono.shuf.fa" "tmp/${base}.rep_mono.reins.fa"

# Extract uppercase (non-repeats), shuffle, strip headers, reinsert, window + label
py_extract_upper "$chr_fasta" > "tmp/${base}_nonrep.fa"
if [[ -s "tmp/${base}_nonrep.fa" ]]; then
  "$biasaway" k -f "tmp/${base}_nonrep.fa" -k 1 --seed "$seed" 1>"tmp/${base}.nonrep_mono.raw.fa" 2>"logs/${base}.biasaway_nonrepmono.log"
  strip_ba_headers "tmp/${base}.nonrep_mono.raw.fa" "tmp/${base}.nonrep_mono.shuf.fa"
  py_reinsert_runs_ordered "$chr_fasta" upper "tmp/${base}.nonrep_mono.shuf.fa" \
  > "tmp/${base}.nonrep_mono.reins.fa" 2> "logs/${base}.reinsertion_nonrepeats.log"
else
  cp "$chr_fasta" "tmp/${base}.nonrep_mono.reins.fa"
fi
reheader_and_wrap "globally_mononucleotide_shuffled_nonrepeats" \
  "tmp/${base}.nonrep_mono.reins.fa" \
  "nonrepeat_mononucleotide_global_shuffled/${base}.global_mononuc_shuffled_nonrepeats.fasta"
rm -f "tmp/${base}_nonrep.fa" "tmp/${base}.nonrep_mono.raw.fa" "tmp/${base}.nonrep_mono.shuf.fa" "tmp/${base}.nonrep_mono.reins.fa"

# Check: Should be >0 and leftover_fragments=0 if counts match
grep reinsert-ordered "logs/${base}.reinsertion_repeats.log"
grep reinsert-ordered "logs/${base}.reinsertion_nonrepeats.log"


# Repeat-only / Nonrepeat-only reversal
py_rev_lower "$chr_fasta" > "tmp/${base}.rep_rev.fa"
reheader_and_wrap "reversed_repeats" "tmp/${base}.rep_rev.fa" \
  "repeat_reversed/${base}.reversed_repeats.fasta"
rm -f "tmp/${base}.rep_rev.fa"

py_rev_upper "$chr_fasta" > "tmp/${base}.nonrep_rev.fa"
reheader_and_wrap "reversed_nonrepeats" "tmp/${base}.nonrep_rev.fa" \
  "nonrepeat_reversed/${base}.reversed_nonrepeats.fasta"
rm -f "tmp/${base}.nonrep_rev.fa"

# Holdout copies
for h in "${holdout_chrs[@]}"; do
  if [[ "$base" == "$h" ]]; then
    cp -f "T2T_forward/${base}.fasta"           "holdout/forward/${base}.fasta"
    cp -f "T2T_reversed/${base}.reversed.fasta" "holdout/reversed/${base}.reversed.fasta"
  fi
done

echo "Done $base"

