#!/usr/bin/env python3
"""Collect per-sample QC JSONs from pipeline stages into MultiQC custom content.

Decision 5C: Full pipeline summary — read QC, segmentation, annotation, spatial
metrics from each step.

Decision 6C: Table overview + expandable per-metric plots for cross-sample view.

Decision 7B: Each module emits its own QC JSON; this collector process renders
the summary for MultiQC.

Input:  directory of per-sample QC JSON files (one per stage per sample).
Output: MultiQC custom content YAML files — one table + one bargraph section.
"""

import argparse
import json
import sys
from pathlib import Path

import yaml


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_qc_jsons(paths: list[Path]) -> list[dict]:
    """Load QC JSON files from an explicit list of paths, returning a flat list."""
    results = []
    for p in sorted(paths):
        try:
            data = json.loads(p.read_text())
            if isinstance(data, list):
                results.extend(data)
            else:
                results.append(data)
        except (json.JSONDecodeError, OSError) as e:
            print(f"WARNING: skipping {p}: {e}", file=sys.stderr)
    return results


def group_by_stage(records: list[dict]) -> dict[str, list[dict]]:
    """Group QC records by their 'stage' field."""
    groups: dict[str, list[dict]] = {}
    for rec in records:
        stage = rec.get("stage", "unknown")
        groups.setdefault(stage, []).append(rec)
    return groups


def build_general_stats_table(records: list[dict]) -> dict:
    """Build a MultiQC general_stats table from all QC records.

    Returns a dict suitable for writing as a *_mqc.yaml custom content file.
    """
    # Collect all metrics across all stages per sample
    sample_metrics: dict[str, dict] = {}
    for rec in records:
        sample = rec.get("sample_id", "unknown")
        metrics = rec.get("metrics", {})
        stage = rec.get("stage", "unknown")
        if sample not in sample_metrics:
            sample_metrics[sample] = {}
        for key, val in metrics.items():
            # Prefix metric with stage to avoid collisions
            sample_metrics[sample][f"{stage}_{key}"] = val

    # Build headers and data
    # Discover all metric keys
    all_keys = set()
    for m in sample_metrics.values():
        all_keys.update(m.keys())

    headers = {}
    for key in sorted(all_keys):
        # Determine header attributes
        stage_prefix = key.split("_")[0] if "_" in key else ""
        metric_name = key.split("_", 1)[1] if "_" in key else key
        headers[key] = {
            "title": metric_name.replace("_", " ").title(),
            "description": f"{stage_prefix}: {metric_name}",
            "namespace": stage_prefix.replace("_", " ").title() if stage_prefix else "Pipeline",
            "scale": "RdYlGn" if "fraction" in key or "rate" in key else "Blues",
        }
        # Format percentages
        if "fraction" in key or "rate" in key:
            headers[key]["max"] = 1.0
            headers[key]["min"] = 0.0
            headers[key]["format"] = "{:,.4f}"
        elif "count" in key or "total" in key or "num" in key:
            headers[key]["format"] = "{:,.0f}"

    return {
        "id": "spatialaxe_pipeline_summary",
        "section_name": "Pipeline Summary",
        "description": "Cross-sample summary of all spatialaxe pipeline stages.",
        "plot_type": "generalstats",
        "pconfig": [headers],
        "data": sample_metrics,
    }


def build_bargraph_section(grouped: dict[str, list[dict]]) -> dict:
    """Build a MultiQC bargraph section showing key metrics per stage."""
    categories = {}
    data = {}

    for stage, records in grouped.items():
        for rec in records:
            sample = rec.get("sample_id", "unknown")
            metrics = rec.get("metrics", {})
            if sample not in data:
                data[sample] = {}
            for key, val in metrics.items():
                cat_key = f"{stage}_{key}"
                if cat_key not in categories:
                    categories[cat_key] = {"name": f"{stage}: {key.replace('_', ' ')}"}
                if isinstance(val, (int, float)):
                    data[sample][cat_key] = val

    return {
        "id": "spatialaxe_stage_metrics",
        "section_name": "Stage Metrics",
        "description": "Per-stage metrics across all samples in the pipeline run.",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "spatialaxe_stage_metrics_plot",
            "title": "spatialaxe: Stage Metrics",
            "ylab": "Value",
        },
        "categories": categories,
        "data": data,
    }


def build_per_stage_tables(grouped: dict[str, list[dict]]) -> list[dict]:
    """Build one MultiQC table section per pipeline stage."""
    sections = []
    for stage, records in sorted(grouped.items()):
        # Collect all metric keys for this stage
        all_keys = set()
        for rec in records:
            all_keys.update(rec.get("metrics", {}).keys())

        headers = {}
        for key in sorted(all_keys):
            headers[key] = {
                "title": key.replace("_", " ").title(),
                "description": key,
            }
            if "fraction" in key or "rate" in key:
                headers[key]["max"] = 1.0
                headers[key]["min"] = 0.0
                headers[key]["format"] = "{:,.4f}"
                headers[key]["scale"] = "RdYlGn"
            elif "count" in key or "total" in key or "num" in key:
                headers[key]["format"] = "{:,.0f}"

        data = {}
        for rec in records:
            sample = rec.get("sample_id", "unknown")
            data[sample] = rec.get("metrics", {})

        sections.append({
            "id": f"spatialaxe_{stage}",
            "section_name": f"{stage.replace('_', ' ').title()}",
            "description": f"Detailed metrics from the {stage} stage.",
            "plot_type": "table",
            "pconfig": {
                "id": f"spatialaxe_{stage}_table",
                "title": f"spatialaxe: {stage.replace('_', ' ').title()}",
            },
            "headers": headers,
            "data": data,
        })
    return sections


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Collect QC JSONs and produce MultiQC custom content."
    )
    parser.add_argument(
        "--qc-files",
        required=True,
        nargs="+",
        help="One or more QC JSON file paths from pipeline stages.",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory to write MultiQC custom content YAML files.",
    )
    parser.add_argument(
        "--merged-output",
        default="all_qc_records.json",
        help="Path to write the merged, flattened QC records as JSON.",
    )
    args = parser.parse_args()

    qc_files = [Path(f) for f in args.qc_files]
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records = load_qc_jsons(qc_files)
    Path(args.merged_output).write_text(json.dumps(records, indent=2))

    if not records:
        print("WARNING: no QC JSON files found; writing empty summary.", file=sys.stderr)
        # Write an empty placeholder
        (out_dir / "spatialaxe_pipeline_summary_mqc.yaml").write_text(
            yaml.dump({"id": "spatialaxe_pipeline_summary", "section_name": "Pipeline Summary",
                       "description": "No QC data collected.", "data": {}})
        )
        return

    grouped = group_by_stage(records)

    # 1. General stats table (Decision 6C: table overview)
    general_stats = build_general_stats_table(records)
    (out_dir / "spatialaxe_pipeline_summary_mqc.yaml").write_text(
        yaml.dump(general_stats, default_flow_style=False, sort_keys=False)
    )

    # 2. Bargraph overview (Decision 6C: expandable per-metric plots)
    bargraph = build_bargraph_section(grouped)
    (out_dir / "spatialaxe_stage_metrics_mqc.yaml").write_text(
        yaml.dump(bargraph, default_flow_style=False, sort_keys=False)
    )

    # 3. Per-stage detail tables
    for section in build_per_stage_tables(grouped):
        stage_id = section["id"]
        (out_dir / f"{stage_id}_mqc.yaml").write_text(
            yaml.dump(section, default_flow_style=False, sort_keys=False)
        )

    print(f"Wrote {2 + len(grouped)} MultiQC custom content files to {out_dir}")


if __name__ == "__main__":
    main()
