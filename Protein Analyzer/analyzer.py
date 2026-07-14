"""
Batch Protein Property Analyzer - Streamlit App
Install: pip install biopython pandas streamlit
Run:     streamlit run analyzer.py
"""

import pandas as pd
import streamlit as st

from protein_pipeline import (
    CONSERVATIVE_GROUPS,
    HAS_LOCALCIDER,
    read_csv_safe,
    analyze_sequences_parallel,
    count_conservative_combos,
    count_random_combos,
    ensure_analysis_columns,
    expand_rows_with_mutations,
    max_unique_conservative_variants,
    merge_analysis_results,
    parse_mutation_regions,
    screen_analysis_columns,
)
from ui_constants import APP_CSS, EXPECTED_FORMAT_EXAMPLE


def render_header():
    st.title("Protein property analyzer")
    st.caption("Upload a CSV with sequence, mutation region, and copy count columns.")
    st.divider()


def render_upload_panel():
    col_up, col_hint = st.columns([2, 1])
    with col_up:
        uploaded = st.file_uploader("Upload CSV", type=["csv"], label_visibility="collapsed")
    with col_hint:
        st.markdown("**Expected format**")
        st.code(EXPECTED_FORMAT_EXAMPLE, language="text")
        st.caption(
            "First column = sequence, second = mutation region (e.g. `5-7` or `3,6,10`), "
            "third = number of copies (ignored in scan mode)."
        )
    return uploaded


def read_input_df(uploaded):
    try:
        return read_csv_safe(uploaded, comment="#")
    except Exception as e:
        st.error(f"Could not read file: {e}")
        return None


def resolve_required_columns(df):
    if df.shape[1] < 3:
        st.error("CSV must contain at least 3 columns: sequence, mutation_region, num_copies.")
        return None
    seq_col, region_col, copies_col = df.columns[:3]
    st.success(
        f"Using columns: sequence=`{seq_col}`, "
        f"mutation region=`{region_col}`, copies=`{copies_col}`"
    )
    return seq_col, region_col, copies_col


def render_combo_counter(df, seq_col, region_col, mutation_mode):
    st.markdown("#### Possible combinations in your dataset")
    rows_info = []
    for i, (_, row) in enumerate(df.iterrows()):
        seq = str(row[seq_col]).strip().upper().replace(" ", "")
        region_text = "" if pd.isna(row[region_col]) else str(row[region_col]).strip()
        try:
            regions = parse_mutation_regions(region_text)
            if mutation_mode in ("conservative", "scan"):
                total, breakdown = count_conservative_combos(seq, regions)
                per_pos = ", ".join(
                    f"pos {p}: {aa}->[{','.join(subs)}] ({n} choices)"
                    for p, aa, subs, n in breakdown
                )
            else:
                total, n_pos = count_random_combos(seq, regions)
                per_pos = f"{n_pos} positions x 20 AAs each"
            rows_info.append(
                {
                    "Row": i + 1,
                    "Sequence (truncated)": seq[:30] + ("..." if len(seq) > 30 else ""),
                    "Region": region_text,
                    "Total combinations": f"{total:,}",
                    "Breakdown": per_pos,
                }
            )
        except Exception as e:
            rows_info.append(
                {
                    "Row": i + 1,
                    "Sequence (truncated)": seq[:30],
                    "Region": region_text,
                    "Total combinations": "error",
                    "Breakdown": str(e),
                }
            )

    st.dataframe(pd.DataFrame(rows_info), width="stretch", hide_index=True)
    st.caption(
        "Conservative/scan combos count each position's substitute pool size + 1 "
        "(keeping original). Random combos = 20^n_positions."
    )


def render_mutation_controls():
    st.divider()
    st.subheader("Mutation settings")

    mode = st.radio(
        "Substitution mode",
        options=["Random", "Conservative", "Scan (single-position exhaustive)"],
        horizontal=True,
        help=(
            "**Random**: each position mutates to any other AA at random.\n\n"
            "**Conservative**: each position swaps to a physicochemically similar AA "
            "(polar->polar, nonpolar->nonpolar, etc.). Number of copies set per row.\n\n"
            "**Scan**: for each position in the region, generate every possible "
            "conservative substitute, one mutation at a time. "
            "Ignores the copies column. Best for identifying which residue matters."
        ),
    )

    mode_key = "scan" if "Scan" in mode else mode.lower()

    if mode_key in ("conservative", "scan"):
        with st.expander("View conservative substitution groups"):
            rows = [
                {"Residue": aa, "Conservative substitutes": ", ".join(subs)}
                for aa, subs in sorted(CONSERVATIVE_GROUPS.items())
            ]
            st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

    if mode_key == "scan":
        st.info(
            "Scan mode generates all single-position conservative variants. "
            "Each output row differs from the original by exactly one residue. "
            "The `copies` column in your CSV is ignored."
        )

    random_seed = st.number_input(
        "Random seed (used in random/conservative modes)",
        min_value=0,
        max_value=999999,
        value=42,
        step=1,
        disabled=(mode_key == "scan"),
    )

    st.divider()
    st.subheader("Analysis settings")
    include_extended_charge_metrics = st.checkbox(
        "Include FCR and NCPR",
        value=True,
        help="These are fast and usually fine to leave on.",
    )
    include_cider_kappa = st.checkbox(
        "Include CIDER kappa (slow)",
        value=False,
        help="This is the most expensive metric by far on large batches.",
    )

    if include_cider_kappa and not HAS_LOCALCIDER:
        st.warning(
            "`localcider` is not installed in this environment, so the `cider_kappa` "
            "column will be blank."
        )

    run_analysis = st.button("Generate variants and analyze")
    return mode_key, random_seed, include_extended_charge_metrics, include_cider_kappa, run_analysis


def render_summary(results_df, total_count, mutation_mode):
    n_ok = results_df["error"].fillna("").eq("").sum()
    n_err = total_count - n_ok
    ok = results_df[results_df["error"].fillna("") == ""]

    st.subheader(f"Results - {n_ok} of {total_count} sequences processed")
    if n_err:
        st.warning(f"{n_err} sequence(s) had errors. Check the `error` column in the download.")

    if ok.empty:
        return

    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Variants", n_ok)
    m2.metric("Avg MW", f"{ok['molecular_weight'].mean() / 1000:.1f} kDa")
    m3.metric("Avg GRAVY", f"{ok['gravy'].mean():.3f}")
    m4.metric("Avg pI", f"{ok['isoelectric_point'].mean():.2f}")
    m5.metric("Stable", f"{ok['stable'].sum()} / {n_ok}")

    st.divider()
    st.subheader("Distributions")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**GRAVY score**")
        st.bar_chart(ok["gravy"], height=160, width="stretch")
        st.markdown("**Isoelectric point (pI)**")
        st.bar_chart(ok["isoelectric_point"], height=160, width="stretch")
    with c2:
        st.markdown("**Molecular weight (Da)**")
        st.bar_chart(ok["molecular_weight"], height=160, width="stretch")
        st.markdown("**Instability index**")
        st.bar_chart(ok["instability_index"], height=160, width="stretch")
    st.divider()


def render_outputs(out_df):
    st.subheader("Full results table")
    display_df = out_df.copy()
    if "stable" in display_df.columns:
        display_df["stable"] = display_df["stable"].map({True: "stable", False: "unstable"}).fillna(
            display_df["stable"]
        )

    max_display_rows = 5000
    if len(display_df) > max_display_rows:
        st.info(
            f"Showing first {max_display_rows:,} of {len(display_df):,} results "
            "(to avoid browser timeout). Download CSV for the full dataset."
        )
        display_df = display_df.head(max_display_rows)

    st.dataframe(display_df, width="stretch", height=400)

    csv_bytes = out_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download results CSV",
        data=csv_bytes,
        file_name="protein_results.csv",
        mime="text/csv",
    )


def main():
    st.set_page_config(page_title="Protein Analyzer", page_icon="P", layout="wide")
    st.markdown(APP_CSS, unsafe_allow_html=True)
    render_header()

    uploaded = render_upload_panel()
    if not uploaded:
        st.info("Upload a CSV file above to get started.")
        return

    df = read_input_df(uploaded)
    if df is None:
        return

    cols = resolve_required_columns(df)
    if cols is None:
        return
    seq_col, region_col, copies_col = cols

    (
        mutation_mode,
        random_seed,
        include_extended_charge_metrics,
        include_cider_kappa,
        run_analysis,
    ) = render_mutation_controls()

    with st.expander("How many combinations are possible for your sequences?", expanded=False):
        render_combo_counter(df, seq_col, region_col, mutation_mode)

    if not run_analysis:
        st.info("Configure settings above and click **Generate variants and analyze**.")
        return

    working_df, skipped_zero_copy_rows, duplicate_warnings = expand_rows_with_mutations(
        df=df,
        seq_col=seq_col,
        region_col=region_col,
        copies_col=copies_col,
        random_seed=random_seed,
        mutation_mode=mutation_mode,
    )

    if working_df.empty:
        st.error("No sequences to analyze after expansion. Check your region and copies values.")
        return
    if skipped_zero_copy_rows:
        st.warning(f"{skipped_zero_copy_rows} input row(s) skipped (copy count was 0).")

    if duplicate_warnings:
        for w in duplicate_warnings:
            if w["mode"] == "random":
                n_pos = len(w["per_position"])
                st.warning(
                    f"Row {w['row']}: you requested **{w['requested']} copies** but only "
                    f"**{w['max_unique']:,} unique** random variants exist "
                    f"({n_pos} position{'s' if n_pos != 1 else ''} x 19 possible substitutes each). "
                    f"The extra copies will be duplicates."
                )
            else:
                pos_detail = ", ".join(
                    f"pos {p} ({aa}: {n} substitute{'s' if n != 1 else ''})"
                    for p, aa, n in w["per_position"]
                )
                st.warning(
                    f"Row {w['row']}: you requested **{w['requested']} copies** but only "
                    f"**{w['max_unique']} unique** conservative variants exist for this region "
                    f"({pos_detail}). The extra copies will be duplicates."
                )

    force_recompute_mask = pd.Series(False, index=working_df.index)
    if "original_sequence" in working_df.columns:
        force_recompute_mask = (
            working_df["original_sequence"].astype(str).str.strip()
            != working_df[seq_col].astype(str).str.strip()
        )

    analysis_plan = screen_analysis_columns(
        working_df,
        include_extended_charge_metrics=include_extended_charge_metrics,
        include_cider_kappa=include_cider_kappa,
        force_recompute_mask=force_recompute_mask,
    )
    rows_to_update_mask = analysis_plan["rows_to_update_mask"]
    rows_to_update_count = analysis_plan["rows_to_update_count"]

    if analysis_plan["missing_columns"]:
        st.caption(f"Missing analyzer columns to add: {', '.join(analysis_plan['missing_columns'])}")
    if analysis_plan["partially_empty_columns"]:
        st.caption(
            f"Analyzer columns with blanks to fill: {', '.join(analysis_plan['partially_empty_columns'])}"
        )
    if analysis_plan["forced_recompute_count"]:
        st.caption(
            f"{analysis_plan['forced_recompute_count']:,} row(s) will be recomputed because the "
            "sequence changed during mutation expansion."
        )

    if rows_to_update_count == 0:
        st.info("All requested analyzer columns are already populated. Nothing needed.")
        out_df = ensure_analysis_columns(
            working_df,
            include_extended_charge_metrics=include_extended_charge_metrics,
            include_cider_kappa=include_cider_kappa,
        )
    else:
        with st.spinner(f"Analyzing {rows_to_update_count:,} sequences..."):
            progress_bar = st.progress(0)
            status_text = st.empty()

            def update_progress(completed, total):
                progress_bar.progress(completed / total)
                status_text.text(f"Analyzed {completed:,} / {total:,} sequences")

            results = analyze_sequences_parallel(
                working_df.loc[rows_to_update_mask, seq_col].astype(str),
                progress_callback=update_progress,
                include_extended_charge_metrics=include_extended_charge_metrics,
                include_cider_kappa=include_cider_kappa,
            )

            progress_bar.progress(1.0)
            status_text.text(f"Completed {rows_to_update_count:,} sequences")

        fresh_results_df = pd.DataFrame(results, index=working_df.index[rows_to_update_mask])
        out_df = merge_analysis_results(
            working_df,
            fresh_results_df,
            include_extended_charge_metrics=include_extended_charge_metrics,
            include_cider_kappa=include_cider_kappa,
            force_recompute_mask=force_recompute_mask,
        )
        out_df = ensure_analysis_columns(
            out_df,
            include_extended_charge_metrics=include_extended_charge_metrics,
            include_cider_kappa=include_cider_kappa,
        )

    results_df = out_df[analysis_plan["requested_columns"]].copy()

    render_summary(results_df, len(working_df), mutation_mode)
    render_outputs(out_df)


if __name__ == "__main__":
    main()
