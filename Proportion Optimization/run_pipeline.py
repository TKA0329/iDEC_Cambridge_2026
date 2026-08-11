#!/usr/bin/env python3
"""
run_pipeline.py

Orchestrates the full degenerate-motif optimization workflow for the iDEC
pipeline, sweeping the --proportion parameter of filter_sequences_pt2.py and,
for each value, running the whole downstream chain to see how many resulting
protein sequences pass the final property filters:

    filter_sequences_pt2.py  (build degenerate protein motif at proportion p)
            |
            v
    reverse_translation_V5_StopSafe.py  (--mode direct, --bins K; produces one
                                          degenerate DNA/IUPAC sequence)
            |
            v
    sample_degenerate.py  (sample N concrete DNA sequences, translate to
                            protein)
            |
            v
    batch_analyzer.py  (compute ProtParam / charge / CIDER properties)
            |
            v
    filter_sequences.py  (apply the SAME thresholds at every proportion --
                          if you didn't pass --*-max/--*-min flags, you're
                          prompted interactively once, for the first
                          proportion, and that choice is reused for the rest)

For every proportion this writes a per-proportion folder with every
intermediate file, and at the end writes a single summary.csv with one row
per proportion (proportion, degenerate_sequence, encoded_space, samples
requested/generated, rows analyzed, rows passed, pct pass, and
est_successful_sequences), sorted by est_successful_sequences descending, so
you can see directly which proportion produced the best degenerate design.

est_successful_sequences = (pct_pass / 100) * encoded_space -- i.e. the
estimated number of distinct sequences in that proportion's FULL degenerate
motif space that would pass the filter, not just the pass rate among the
samples drawn. pct_pass alone isn't comparable across proportions on its
own: a tiny, tightly-constrained motif can post a high pct_pass just because
it barely explores any sequence space, while a much larger, more permissive
motif with a lower pct_pass may still contain far more usable designs.
est_successful_sequences (an estimated count, extrapolated from the sampled
pass rate) is the metric this script optimizes for when picking "best".

All defaults below can be overridden from the command line.

USAGE (minimal, using all defaults):
    python run_pipeline.py input.fasta

USAGE (with the filter thresholds you actually care about, and a custom
sweep range):
    python run_pipeline.py input.fasta my_run \\
        --prop-min 0.1 --prop-max 0.9 --prop-step 0.1 \\
        --num-samples 100000 --seed 42 \\
        --gravy-max 0.5 --instability-max 40 --ncpr-min -0.3 --ncpr-max 0.3

USAGE (explicit list of proportions instead of a min/max/step sweep):
    python run_pipeline.py input.fasta my_run --proportions 0.2,0.35,0.5,0.8

All the filter_sequences.py threshold flags (--gravy-max, --instability-max,
--helix-sheet-max, --turn-fraction-min, --charge-sum-min,
--charge-ph7-exclude-min/--charge-ph7-exclude-max, --fcr-max,
--ncpr-min/--ncpr-max, --cider-kappa-max, --pace-scholtz-sum-max,
--pace-scholtz-mean-max, --hydrophobic-moment-min/--hydrophobic-moment-max,
--helical-face-occupancy-min) are entered ONCE and reused
identically for every proportion in the sweep, so the resulting %-pass
numbers are directly comparable across proportions. You can supply them on
the command line, or -- if you omit them -- run_pipeline.py will prompt you
interactively (the same prompts filter_sequences.py shows when run by hand)
exactly once, for the first proportion, then lock in and reuse whatever you
chose for every remaining proportion. Pass --skip-filter-prompt to instead
keep the old headless behavior of passing every row through unfiltered when
no thresholds are given (useful for cron/CI where there's no terminal to
prompt on).
"""

import argparse
import csv
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent


def load_hotspot_module():
    """Import hotspot_detection.py from the same directory as this script,
    if present. Kept as a separate optional file (not vendored in) so
    --hotspot-detection is an opt-in drop-in dependency rather than a hard
    requirement of run_pipeline.py. Returns None (with a warning printed by
    the caller) if the file isn't found."""
    path = SCRIPT_DIR / "hotspot_detection.py"
    if not path.exists():
        return None
    spec = importlib.util.spec_from_file_location("hotspot_detection", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fmt_prop(p):
    """Filesystem-friendly, stable label for a proportion value."""
    return f"{p:.4f}".rstrip("0").rstrip(".")


def parse_proportions(args):
    if args.proportions:
        try:
            props = [float(x.strip()) for x in args.proportions.split(",") if x.strip()]
        except ValueError:
            raise SystemExit(f"--proportions must be a comma-separated list of numbers, got: {args.proportions}")
    else:
        props = []
        p = args.prop_min
        # Use round() to avoid float drift (0.1 + 0.1 + ... != 0.3 exactly)
        n_steps = round((args.prop_max - args.prop_min) / args.prop_step)
        for i in range(n_steps + 1):
            props.append(round(args.prop_min + i * args.prop_step, 10))
    for p in props:
        if not (0.0 <= p <= 1.0):
            raise SystemExit(f"Proportion values must be in [0, 1]; got {p}")
    return props


def run(cmd, log_path=None, description=""):
    """Run a subprocess, streaming/saving combined stdout+stderr. Raises on
    non-zero exit (caller decides whether to catch)."""
    print(f"    $ {' '.join(str(c) for c in cmd)}")
    with open(log_path, "w", encoding="utf-8") if log_path else open(
        subprocess.os.devnull, "w"
    ) as logf:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        logf.write(proc.stdout or "")
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout or "").splitlines()[-25:])
        raise RuntimeError(
            f"{description} failed (exit {proc.returncode}). Last output:\n{tail}\n"
            f"(full log: {log_path})"
        )
    return proc.stdout or ""


def run_interactive(cmd, description=""):
    """Run a subprocess with stdin/stdout/stderr inherited from the parent
    (i.e. NOT captured), so the user can see and answer prompts exactly as
    if they'd typed the command themselves. Used for filter_sequences.py's
    interactive mode, which the plain `run()` helper can't support since it
    pipes stdout to a log file (any input() call inside the subprocess would
    have nothing visible to prompt against)."""
    print(f"    $ {' '.join(str(c) for c in cmd)}")
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError(f"{description} failed (exit {proc.returncode}).")


def prompt_for_hotspot_k(ks, default_k):
    """Ask the user which k value (from the just-printed sweep) to lock in
    for the rest of the proportion sweep. Blank input keeps --hotspot-k's
    current default rather than forcing a choice."""
    while True:
        raw = input(
            f"Choose a --hotspot-k value to use for this and all remaining "
            f"proportions [default {default_k}]: "
        ).strip()
        if not raw:
            return default_k
        try:
            return float(raw)
        except ValueError:
            print(f"'{raw}' isn't a number -- enter one of {ks}, any other numeric "
                  f"value, or leave blank for the default ({default_k}).")


def stdin_is_interactive():
    """True if there's a real terminal to prompt on. False in headless
    contexts (cron, CI, piped input) where input() would just raise
    EOFError."""
    try:
        return sys.stdin.isatty()
    except (AttributeError, ValueError):
        return False


def metadata_to_filter_flags(meta):
    """Convert the threshold values recorded in filter_sequences.py's
    metadata comments (written by write_output_with_metadata()) back into
    --*-max/--*-min CLI flags, so the thresholds a user chose interactively
    can be reapplied identically, non-interactively, to every other
    proportion in the sweep. Any threshold that was skipped (blank at the
    prompt) is recorded as the literal string "None" and is omitted here."""
    def _get(key):
        val = meta.get(key)
        if val is None or val == "" or val == "None":
            return None
        return val

    flags = []
    single_value_map = [
        ("--gravy-max", "gravy_max"),
        ("--instability-max", "instability_index_max"),
        ("--helix-sheet-max", "helix_sheet_sum_max"),
        ("--turn-fraction-min", "turn_fraction_min"),
        ("--charge-sum-min", "charge_sum_min"),
        ("--fcr-max", "fcr_max"),
        ("--cider-kappa-max", "cider_kappa_max"),
        ("--pace-scholtz-sum-max", "pace_scholtz_sum_max"),
        ("--pace-scholtz-mean-max", "pace_scholtz_mean_max"),
        ("--helical-face-occupancy-min", "helical_face_occupancy_min"),
    ]
    for flag, key in single_value_map:
        val = _get(key)
        if val is not None:
            flags += [flag, str(val)]

    ph7_min, ph7_max = _get("charge_at_pH7_min"), _get("charge_at_pH7_max")
    if ph7_min is not None and ph7_max is not None:
        flags += ["--charge-ph7-exclude-min", str(ph7_min), "--charge-ph7-exclude-max", str(ph7_max)]

    ncpr_min, ncpr_max = _get("ncpr_min"), _get("ncpr_max")
    if ncpr_min is not None and ncpr_max is not None:
        flags += ["--ncpr-min", str(ncpr_min), "--ncpr-max", str(ncpr_max)]

    hm_min, hm_max = _get("hydrophobic_moment_min"), _get("hydrophobic_moment_max")
    if hm_min is not None and hm_max is not None:
        flags += ["--hydrophobic-moment-min", str(hm_min), "--hydrophobic-moment-max", str(hm_max)]

    return flags


def parse_metadata_comments(csv_path):
    """Read the leading '# key: value' comment lines written by
    filter_sequences.py's write_output_with_metadata()."""
    meta = {}
    with open(csv_path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                break
            line = line[1:].strip()
            if ":" not in line:
                continue
            key, _, value = line.partition(":")
            meta[key.strip()] = value.strip()
    return meta


def count_csv_rows(csv_path):
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.reader(line for line in f if not line.startswith("#"))
        rows = list(reader)
    return max(0, len(rows) - 1)  # minus header


def read_degenerate_sequence(fasta_path):
    seq_lines = []
    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
            break  # only the first record's sequence (bins=1 in this pipeline)
    return seq_lines[0] if seq_lines else ""


def parse_encoded_space(rt_log_text):
    m = re.search(r"k=\d+ total encoded space\s*:\s*(\d+)", rt_log_text)
    return int(m.group(1)) if m else None


def build_parser():
    p = argparse.ArgumentParser(
        description="Sweep filter_sequences_pt2 proportions through the full "
                     "RT -> sample -> analyze -> filter chain and report %% pass per proportion.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("input_fasta", type=Path, help="Aligned protein FASTA (input to filter_sequences_pt2.py)")
    p.add_argument("output_dir", type=Path, nargs="?", default=Path("pipeline_output"),
                    help="Root output directory (default: pipeline_output)")

    # ── Proportion sweep ─────────────────────────────────────────
    sweep = p.add_argument_group("proportion sweep (filter_sequences_pt2.py)")
    sweep.add_argument("--prop-min", type=float, default=0.1)
    sweep.add_argument("--prop-max", type=float, default=0.9)
    sweep.add_argument("--prop-step", type=float, default=0.1)
    sweep.add_argument("--proportions", type=str, default=None,
                        help="Comma-separated explicit list of proportions, overrides --prop-min/max/step")
    sweep.add_argument("--random-base", action="store_true",
                        help="Passthrough to filter_sequences_pt2.py --random-base")

    # ── Reverse translation ──────────────────────────────────────
    rt = p.add_argument_group("reverse translation (reverse_translation_V5_StopSafe.py)")
    rt.add_argument("--rt-bins", type=int, default=1,
                     help="Passthrough to RT --bins (default: 1). Mode is fixed to 'direct' for this pipeline.")
    rt.add_argument("--rt-sample-size", type=int, default=1000,
                     help="Passthrough to RT --sample-size (default: 1000; 0 = no sampling/use all sequences)")
    rt.add_argument("--allow-stop-codons", action="store_true",
                     help="Passthrough to RT --allow-stop-codons")

    # ── Sampling ─────────────────────────────────────────────────
    samp = p.add_argument_group("degenerate sampling (sample_degenerate.py)")
    samp.add_argument("--num-samples", type=int, default=100_000,
                       help="Number of concrete sequences to sample per proportion (default: 100000)")
    samp.add_argument("--no-unique", action="store_true",
                       help="Passthrough to sample_degenerate.py --no-unique (allow duplicate DNA samples)")

    # ── Batch analysis ───────────────────────────────────────────
    ba = p.add_argument_group("batch analysis (batch_analyzer.py)")
    ba.add_argument("--workers", type=int, default=None)
    ba.add_argument("--processes", action="store_true")
    ba.add_argument("--cider-kappa", action="store_true")
    ba.add_argument("--no-fcr-ncpr", action="store_true")
    ba.add_argument("--cahs-motif1", action="store_true",
                     help="Passthrough to batch_analyzer.py --cahs-motif1 (score the CAHS "
                          "motif-1 helical segment: Pace & Scholtz helix propensity, "
                          "Eisenberg hydrophobic moment, helical face occupancy, salt-bridge "
                          "count, Pro/Gly helix-breaker flags; default off)")

    # ── Final filter thresholds (entered once, applied to every proportion) ──
    filt = p.add_argument_group(
        "final filter thresholds (filter_sequences.py, applied identically at every proportion)"
    )
    filt.add_argument("--gravy-max", type=float, default=None)
    filt.add_argument("--instability-max", type=float, default=None)
    filt.add_argument("--helix-sheet-max", type=float, default=None)
    filt.add_argument("--turn-fraction-min", type=float, default=None)
    filt.add_argument("--charge-sum-min", type=float, default=None)
    filt.add_argument("--charge-ph7-exclude-min", type=float, default=None)
    filt.add_argument("--charge-ph7-exclude-max", type=float, default=None)
    filt.add_argument("--fcr-max", type=float, default=None)
    filt.add_argument("--ncpr-min", type=float, default=None)
    filt.add_argument("--ncpr-max", type=float, default=None)
    filt.add_argument("--cider-kappa-max", type=float, default=None)
    filt.add_argument("--pace-scholtz-sum-max", type=float, default=None)
    filt.add_argument("--pace-scholtz-mean-max", type=float, default=None)
    filt.add_argument("--hydrophobic-moment-min", type=float, default=None)
    filt.add_argument("--hydrophobic-moment-max", type=float, default=None)
    filt.add_argument("--helical-face-occupancy-min", type=float, default=None)
    filt.add_argument(
        "--skip-filter-prompt", action="store_true",
        help="If no --*-max/--*-min threshold flags are supplied, skip the "
             "interactive filter_sequences.py prompts and pass every row "
             "through unfiltered at every proportion instead (old behavior). "
             "Use this for unattended/headless runs (cron, CI) where no "
             "terminal is available to prompt on.",
    )

    # ── Hotspot detection (optional; restricts variant positions) ──
    hs = p.add_argument_group(
        "hotspot detection (optional; runs ONCE on the input FASTA before the "
        "proportion sweep, and restricts filter_sequences_pt2.py's variant output "
        "to the flagged hotspot region(s) at EVERY proportion)"
    )
    hs.add_argument(
        "--hotspot-detection", dest="hotspot_detection", action="store_true", default=True,
        help="Run positional-variability hotspot detection once on the input FASTA "
             "(requires hotspot_detection.py in the same directory as this script), "
             "before the proportion sweep starts. The flagged hotspot region(s) are "
             "then passed as --hotspot-positions to filter_sequences_pt2.py for "
             "every proportion, so its variant rows (>2, >3, ...) are only emitted "
             "for positions inside those regions. Also adds n_hotspot_regions/"
             "hotspot_regions columns to summary.csv (same value at every "
             "proportion, since detection runs once). On by default; use "
             "--no-hotspot-detection to disable.",
    )
    hs.add_argument(
        "--no-hotspot-detection", dest="hotspot_detection", action="store_false",
        help="Disable hotspot detection: every position is treated as a candidate "
             "variant position (the old pre-hotspot behavior), and no "
             "hotspot_table.tsv/hotspot_regions.tsv are written.",
    )
    hs.add_argument("--hotspot-metric", choices=["noise", "entropy"], default="noise",
                     help="Per-position score: noise = 100 - top_aa%%, or Shannon entropy "
                          "over the full residue distribution (default: noise)")
    hs.add_argument("--hotspot-detector", choices=["mad", "plain"], default="mad",
                     help="mad = robust median/MAD z-score thresholding (default); "
                          "plain = naive mean/stdev z-score, for comparison")
    hs.add_argument("--hotspot-k", type=float, default=0.5,
                     help="Z-score threshold for flagging a position (default: 0.5)")
    hs.add_argument("--hotspot-max-gap", type=int, default=1,
                     help="Max position gap allowed when merging adjacent flagged "
                          "positions into one region (default: 1)")
    hs.add_argument(
        "--hotspot-k-sweep", type=str, default=None,
        help="Comma-separated k values (e.g. 0.5,1.0,1.5,2.0,2.5,3.0) to sweep on "
             "the input FASTA. Prints region/flagged-position counts at each k, "
             "then prompts you to choose one value. Whatever you choose is locked "
             "in and used (as --hotspot-k) for the single up-front detection run "
             "whose regions then apply to every proportion. Requires "
             "--hotspot-detection and a terminal to prompt on; if no terminal is "
             "available, falls back to --hotspot-k with a warning.")

    # ── Misc ─────────────────────────────────────────────────────
    p.add_argument("--seed", type=int, default=42,
                    help="Random seed, reused for filter_sequences_pt2 (--random-base only), "
                         "RT, and sample_degenerate.py (default: 42)")
    p.add_argument("--continue-on-error", action="store_true", default=True,
                    help="Keep sweeping remaining proportions if one fails (default: on)")
    p.add_argument("--stop-on-error", dest="continue_on_error", action="store_false",
                    help="Abort the whole sweep on the first failed proportion")
    p.add_argument("--python", type=str, default=sys.executable,
                    help="Python interpreter to invoke the sub-scripts with")
    return p


def main():
    args = build_parser().parse_args()

    if (args.charge_ph7_exclude_min is None) != (args.charge_ph7_exclude_max is None):
        raise SystemExit("--charge-ph7-exclude-min and --charge-ph7-exclude-max must be supplied together.")
    if (args.ncpr_min is None) != (args.ncpr_max is None):
        raise SystemExit("--ncpr-min and --ncpr-max must be supplied together.")
    if (args.hydrophobic_moment_min is None) != (args.hydrophobic_moment_max is None):
        raise SystemExit("--hydrophobic-moment-min and --hydrophobic-moment-max must be supplied together.")

    if not args.input_fasta.exists():
        raise SystemExit(f"Input FASTA not found: {args.input_fasta}")

    proportions = parse_proportions(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    filter_cli_flags = []
    for flag, val in [
        ("--gravy-max", args.gravy_max),
        ("--instability-max", args.instability_max),
        ("--helix-sheet-max", args.helix_sheet_max),
        ("--turn-fraction-min", args.turn_fraction_min),
        ("--charge-sum-min", args.charge_sum_min),
        ("--charge-ph7-exclude-min", args.charge_ph7_exclude_min),
        ("--charge-ph7-exclude-max", args.charge_ph7_exclude_max),
        ("--fcr-max", args.fcr_max),
        ("--ncpr-min", args.ncpr_min),
        ("--ncpr-max", args.ncpr_max),
        ("--cider-kappa-max", args.cider_kappa_max),
        ("--pace-scholtz-sum-max", args.pace_scholtz_sum_max),
        ("--pace-scholtz-mean-max", args.pace_scholtz_mean_max),
        ("--hydrophobic-moment-min", args.hydrophobic_moment_min),
        ("--hydrophobic-moment-max", args.hydrophobic_moment_max),
        ("--helical-face-occupancy-min", args.helical_face_occupancy_min),
    ]:
        if val is not None:
            filter_cli_flags += [flag, str(val)]

    # If no thresholds were supplied on the command line, retain
    # filter_sequences.py's normal interactive-prompt behavior instead of
    # silently passing every row through unfiltered: prompt once (via a real,
    # stdio-inherited subprocess call) for the very first proportion, then
    # reuse whatever thresholds the user chose -- as fixed --*-max/--*-min
    # flags -- for every remaining proportion. This keeps the sweep
    # comparable across proportions while still giving the user real control.
    filter_prompt_pending = False
    if filter_cli_flags:
        pass  # thresholds already supplied explicitly; nothing to prompt for
    elif args.skip_filter_prompt:
        print("Warning: no final-filter thresholds were supplied and "
              "--skip-filter-prompt was set, so the final filter step will "
              "pass every row through unfiltered (100% pass at every "
              "proportion).\n")
    elif not stdin_is_interactive():
        print("Warning: no final-filter thresholds were supplied and no "
              "terminal is available to prompt on (stdin is not a TTY), so "
              "the final filter step will pass every row through unfiltered "
              "(100% pass at every proportion). Pass --gravy-max/"
              "--instability-max/etc. explicitly, or run interactively, to "
              "make the sweep meaningful.\n")
    else:
        filter_prompt_pending = True
        print("No final-filter thresholds were supplied on the command "
              "line -- you'll be prompted interactively (same prompts as "
              "running filter_sequences.py by hand) once, for the first "
              "proportion. Whatever thresholds you choose will then be "
              "reused automatically, non-interactively, for every "
              "remaining proportion so the sweep stays comparable.\n")

    print(f"Proportions to sweep ({len(proportions)}): {proportions}\n")

    # Hotspot detection runs ONCE, up front, on the raw input FASTA -- not
    # per-proportion -- since the set of aligned positions/residues doesn't
    # change with --proportion (only the allowed-set thresholding does).
    # Whatever region(s) it flags are then reused, as a fixed
    # --hotspot-positions value, for every proportion in the sweep.
    hotspot_regions_label = ""
    n_hotspot_regions = ""
    hotspot_positions_arg = None  # e.g. "8-11,16-21", passed straight through to pt2

    if args.hotspot_detection:
        hotspot_mod = load_hotspot_module()
        if hotspot_mod is None:
            print("Warning: --hotspot-detection was set but hotspot_detection.py "
                  "was not found next to run_pipeline.py -- skipping hotspot "
                  "analysis for this run.\n")
        else:
            print(f"Hotspot detection enabled (running once on input FASTA): "
                  f"metric={args.hotspot_metric} detector={args.hotspot_detector} "
                  f"k={args.hotspot_k} max_gap={args.hotspot_max_gap}")

            records = hotspot_mod.parse_fasta(args.input_fasta)
            try:
                hs_table = hotspot_mod.build_table(records)
            except ValueError as exc:
                raise SystemExit(f"Hotspot detection failed: {exc}")

            table_path = args.output_dir / "hotspot_table.tsv"
            sweep_path = args.output_dir / "hotspot_k_sweep.tsv"
            regions_path = args.output_dir / "hotspot_regions.tsv"

            chosen_k = args.hotspot_k
            if args.hotspot_k_sweep:
                ks = [float(x) for x in args.hotspot_k_sweep.split(",")]
                if stdin_is_interactive():
                    hotspot_mod.sweep_k(
                        hs_table, args.hotspot_metric, args.hotspot_detector, ks,
                        max_gap=args.hotspot_max_gap, save_csv=sweep_path,
                    )
                    chosen_k = prompt_for_hotspot_k(ks, args.hotspot_k)
                    print(f"  -> using k={chosen_k} for the hotspot run applied to "
                          f"every proportion\n")
                else:
                    print("Warning: --hotspot-k-sweep requires a terminal to prompt "
                          f"on; falling back to --hotspot-k={args.hotspot_k}.\n")

            hs_result = hotspot_mod.run(
                hs_table, args.hotspot_metric, args.hotspot_detector,
                chosen_k, max_gap=args.hotspot_max_gap,
            )

            # Save the annotated per-position table (score/z/hotspot columns).
            hotspot_mod.write_table(
                hotspot_mod.annotate_table(hs_table, hs_result), table_path
            )

            region_strs = [f"{s}-{e}" if s != e else f"{s}" for s, e in hs_result["regions"]]
            n_hotspot_regions = len(hs_result["regions"])
            hotspot_regions_label = ";".join(region_strs)
            hotspot_positions_arg = ",".join(region_strs)  # matches parse_position_spec's format

            # Save the final chosen-k regions table too (independent of
            # whether --hotspot-k-sweep was used -- sweep_k's --save-regions
            # only records the sweep itself, not the k that got locked in).
            with open(regions_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(["region_start", "region_end", "flagged_positions"])
                for start, end in hs_result["regions"]:
                    idxs = [i for i, p in enumerate(hs_result["positions"]) if start <= p <= end]
                    flagged_positions = [hs_result["positions"][i] for i in idxs if hs_result["flagged"][i]]
                    writer.writerow([start, end, ",".join(map(str, flagged_positions))])

            if region_strs:
                print(f"  {n_hotspot_regions} hotspot region(s) at k={chosen_k}: "
                      f"{hotspot_regions_label}")
            else:
                print(f"  No hotspot regions flagged at k={chosen_k} -- every "
                      f"proportion's pt2 step will emit NO variant rows at all "
                      f"(only the base/consensus sequence), since --hotspot-positions "
                      f"will be an empty set.")
            print(f"  Table saved to:   {table_path}")
            if args.hotspot_k_sweep and stdin_is_interactive():
                print(f"  K-sweep saved to: {sweep_path}")
            print(f"  Regions saved to: {regions_path}")
            print("  These flagged region(s) will be applied (as --hotspot-positions) "
                  "to filter_sequences_pt2.py for every proportion in the sweep.\n")

    summary_rows = []

    for p in proportions:
        label = fmt_prop(p)
        print(f"=== proportion {p} ===")
        pdir = args.output_dir / f"prop_{label}"
        pdir.mkdir(parents=True, exist_ok=True)
        row = {"proportion": p}

        try:
            # 1. filter_sequences_pt2.py -----------------------------------
            print("  [1/5] filter_sequences_pt2.py (build degenerate motif)")
            pt2_cmd = [
                args.python, str(SCRIPT_DIR / "filter_sequences_pt2.py"),
                str(args.input_fasta), str(pdir / "variants.fasta"),
                "--report", str(pdir / "pt2_report.tsv"),
                "--proportion", str(p),
            ]
            if args.random_base:
                pt2_cmd += ["--random-base", "--seed", str(args.seed)]
            if hotspot_positions_arg is not None:
                pt2_cmd += ["--hotspot-positions", hotspot_positions_arg]
            run(pt2_cmd, pdir / "01_pt2.log", "filter_sequences_pt2.py")

            # Hotspot regions were computed ONCE, up front, on the raw input
            # FASTA (see above) -- the same regions apply to every proportion,
            # so just record that fixed result here rather than recomputing.
            row["n_hotspot_regions"] = n_hotspot_regions
            row["hotspot_regions"] = hotspot_regions_label

            # 2. reverse_translation_V5_StopSafe.py -------------------------
            print("  [2/5] reverse_translation_V5_StopSafe.py (--mode direct)")
            rt_cmd = [
                args.python, str(SCRIPT_DIR / "reverse_translation_V5_StopSafe.py"),
                str(pdir / "variants.fasta"),
                "--mode", "direct",
                "--bins", str(args.rt_bins),
                "--sample-size", str(args.rt_sample_size),
                "--seed", str(args.seed),
                "--out-seq", str(pdir / "degenerate.fasta"),
                "--log", str(pdir / "02_rt.log"),
            ]
            if args.allow_stop_codons:
                rt_cmd.append("--allow-stop-codons")
            rt_stdout = run(rt_cmd, pdir / "02_rt_subprocess.log", "reverse_translation_V5_StopSafe.py")

            degenerate_seq = read_degenerate_sequence(pdir / "degenerate.fasta")
            row["degenerate_sequence"] = degenerate_seq
            row["encoded_space"] = parse_encoded_space(rt_stdout)

            # 3. sample_degenerate.py ---------------------------------------
            print(f"  [3/5] sample_degenerate.py (--num {args.num_samples})")
            samp_cmd = [
                args.python, str(SCRIPT_DIR / "sample_degenerate.py"),
                "--seq-file", str(pdir / "degenerate.fasta"),
                "--num", str(args.num_samples),
                "--output", str(pdir / "samples.csv"),
                "--seed", str(args.seed),
            ]
            if args.no_unique:
                samp_cmd.append("--no-unique")
            run(samp_cmd, pdir / "03_sample.log", "sample_degenerate.py")
            row["samples_generated"] = count_csv_rows(pdir / "samples.csv")

            # 4. batch_analyzer.py -------------------------------------------
            print("  [4/5] batch_analyzer.py")
            ba_cmd = [
                args.python, str(SCRIPT_DIR / "batch_analyzer.py"),
                str(pdir / "samples.csv"), str(pdir / "analyzed.csv"),
            ]
            if args.workers is not None:
                ba_cmd += ["-j", str(args.workers)]
            if args.processes:
                ba_cmd.append("--processes")
            if args.cider_kappa:
                ba_cmd.append("--cider-kappa")
            if args.no_fcr_ncpr:
                ba_cmd.append("--no-fcr-ncpr")
            if args.cahs_motif1:
                ba_cmd.append("--cahs-motif1")
            run(ba_cmd, pdir / "04_batch_analyzer.log", "batch_analyzer.py")
            row["rows_analyzed"] = count_csv_rows(pdir / "analyzed.csv")

            # 5. filter_sequences.py ------------------------------------------
            if filter_prompt_pending:
                # First proportion, no thresholds supplied yet: run for real,
                # with stdin/stdout inherited so the user sees and answers
                # the actual filter_sequences.py prompts. No --non-interactive
                # and no output capture here (a captured/piped subprocess
                # can't be prompted against).
                print("  [5/5] filter_sequences.py (interactive -- choose thresholds now)")
                filt_cmd = [
                    args.python, str(SCRIPT_DIR / "filter_sequences.py"),
                    str(pdir / "analyzed.csv"), str(pdir / "filtered"),
                ]
                run_interactive(filt_cmd, "filter_sequences.py")
                print("  (05_filter.log skipped for this proportion: ran interactively, output went to your terminal)")
            else:
                print("  [5/5] filter_sequences.py --non-interactive")
                filt_cmd = [
                    args.python, str(SCRIPT_DIR / "filter_sequences.py"),
                    str(pdir / "analyzed.csv"), str(pdir / "filtered"),
                    "--non-interactive",
                ] + filter_cli_flags
                run(filt_cmd, pdir / "05_filter.log", "filter_sequences.py")

            filtered_csv = pdir / "filtered" / f"filtered_analyzed.csv"
            meta = parse_metadata_comments(filtered_csv)

            if filter_prompt_pending:
                # Lock in whatever thresholds the user just chose, and reuse
                # them (non-interactively) for every remaining proportion so
                # pct_pass stays comparable across the whole sweep.
                filter_cli_flags = metadata_to_filter_flags(meta)
                filter_prompt_pending = False
                if filter_cli_flags:
                    print(f"  (thresholds chosen interactively will be reused for remaining proportions: "
                          f"{' '.join(filter_cli_flags)})\n")
                else:
                    print("  (no thresholds were chosen at the prompts -- every row will pass through "
                          "unfiltered for the remaining proportions too)\n")
            row["rows_before_filter"] = meta.get("rows_before")
            row["rows_after_filter"] = meta.get("rows_after")
            row["pct_pass"] = meta.get("pass_rate_pct")
            row["status"] = "ok"

            # est_successful_sequences: extrapolate the sampled pass rate out
            # across the WHOLE degenerate motif's combinatorial space, not
            # just the samples drawn. pct_pass on its own rewards tiny,
            # over-constrained motifs that happen to score well on the few
            # sequences they can even produce; this is the number that
            # actually answers "how many usable designs does this proportion
            # give me".
            enc = row.get("encoded_space")
            try:
                pct = float(row["pct_pass"]) if row.get("pct_pass") not in (None, "") else None
            except (TypeError, ValueError):
                pct = None
            if pct is not None and enc:
                row["est_successful_sequences"] = (pct / 100.0) * enc
                print(f"  -> pct_pass = {row['pct_pass']}%  |  encoded_space = {enc:,}  |  "
                      f"est_successful_sequences ~= {row['est_successful_sequences']:,.0f}\n")
            else:
                row["est_successful_sequences"] = ""
                print(f"  -> pct_pass = {row['pct_pass']}%  "
                      f"(encoded_space unavailable -- est_successful_sequences can't be computed)\n")

        except Exception as exc:
            row["status"] = f"FAILED: {exc}"
            print(f"  !! {exc}\n")
            if not args.continue_on_error:
                summary_rows.append(row)
                break

        summary_rows.append(row)

    # ── Write summary CSV ────────────────────────────────────────────
    fieldnames = [
        "proportion", "status", "degenerate_sequence", "encoded_space",
        "samples_generated", "rows_analyzed", "rows_before_filter",
        "rows_after_filter", "pct_pass", "est_successful_sequences",
        "n_hotspot_regions", "hotspot_regions",
    ]
    summary_path = args.output_dir / "summary.csv"
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    print("=" * 60)
    print(f"Summary written to: {summary_path}")

    ok_rows = [
        r for r in summary_rows
        if r.get("status") == "ok" and r.get("est_successful_sequences") not in (None, "")
    ]
    if ok_rows:
        # Optimize for est_successful_sequences, per the module docstring:
        # pct_pass alone rewards tiny, over-constrained motifs that barely
        # explore any sequence space, while est_successful_sequences
        # extrapolates the sampled pass rate across the full degenerate
        # motif's combinatorial space -- the number that actually answers
        # "how many usable designs does this proportion give me".
        best = max(ok_rows, key=lambda r: float(r["est_successful_sequences"]))
        print(f"Best proportion so far: {best['proportion']}  "
              f"(~{float(best['est_successful_sequences']):,.0f} est. successful sequences, "
              f"{best['pct_pass']}% pass, {best['rows_after_filter']}/{best['rows_before_filter']} rows)")
        print("(This is the best point directly sampled by the sweep -- it is not "
              "a curve fit / interpolated optimum. Narrow --prop-min/--prop-max/"
              "--prop-step around this value and re-run for a finer search, or "
              "ask about adding a regression-based optimizer once you've seen "
              "the shape of this sweep.)")
    else:
        print("No proportion completed successfully, or encoded_space was "
              "unavailable for every proportion (so est_successful_sequences "
              "could not be computed) -- check the per-proportion logs.")
    print("=" * 60)


if __name__ == "__main__":
    main()