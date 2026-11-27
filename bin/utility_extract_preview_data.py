#!/usr/bin/env python3
"""
Extract preview data from Baysor preview HTML reports.

Parses embedded Vega-Lite spec variables and base64 PNG images from the
Baysor preview.html file, writing MultiQC-compatible TSV and PNG files.
"""

import argparse
import base64
import html
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from bs4 import BeautifulSoup


def get_png_files(soup: BeautifulSoup, outdir: Path) -> None:
    """Get png base64 images following specific h1 tags in preview.html"""
    target_ids = ["Transcript_Plots", "Noise_Level"]
    outdir.mkdir(parents=True, exist_ok=True)

    for h1_id in target_ids:
        h1_tag = soup.find("h1", id=h1_id)
        if not h1_tag:
            print(f"[WARN] No <h1> with id {h1_id} found")
            continue

        # Look for the first <img> after the h1 in the DOM
        img_tag = h1_tag.find_next("img")
        if not img_tag or not img_tag.get("src"):
            print(f"[WARN] No <img> found after h1#{h1_id}")
            continue

        img_src = img_tag["src"]
        if img_src.startswith("data:image/png;base64,"):
            base64_data = img_src.split(",", 1)[1]
            data = base64.b64decode(base64_data)
        else:
            print(f"[WARN] img src is not base64 PNG for h1#{h1_id}")
            continue

        # save png files with _mqc suffix for MultiQC integration
        img_name = f"{h1_id}_mqc.png".lower()
        out_path = outdir / img_name
        with open(out_path, "wb") as f:
            f.write(data)

        print(f"[INFO] Saved {img_name}")

    return None


def extract_js_object(text: str, start_idx: int) -> Tuple[Optional[str], int]:
    """Extract json-like object starting at start_idx."""
    if start_idx >= len(text) or text[start_idx] != "{":
        return None, start_idx

    stack, in_str, escape, quote = [], False, False, None
    for i in range(start_idx, len(text)):
        ch = text[i]
        if in_str:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == quote:
                in_str = False
        else:
            if ch in ('"', "'"):
                in_str, quote = True, ch
            elif ch == "{":
                stack.append("{")
            elif ch == "}":
                stack.pop()
                if not stack:
                    return text[start_idx : i + 1], i + 1
            elif ch == "/" and i + 1 < len(text):
                # skip js comments
                nxt = text[i + 1]
                if nxt == "/":
                    end = text.find("\n", i + 2)
                    i = len(text) - 1 if end == -1 else end
                elif nxt == "*":
                    end = text.find("*/", i + 2)
                    if end == -1:
                        break
                    i = end + 1

    return None, start_idx


def js_to_json(js: str) -> str:
    """Convert a JS object string to valid JSON."""
    # Remove comments
    js = re.sub(r"/\*.*?\*/", "", js, flags=re.S)
    js = re.sub(r"//[^\n]*", "", js)

    # Convert single-quoted strings to double-quoted strings
    js = re.sub(
        r"'((?:\\.|[^'\\])*)'",
        lambda m: '"' + m.group(1).replace('"', '\\"') + '"',
        js,
    )

    # Remove trailing commas
    js = re.sub(r",\s*(?=[}\]])", "", js)
    js = re.sub(r",\s*,+", ",", js)

    return js.strip()


def find_variables(script_text: str) -> Dict[str, str]:
    """Find all 'var|let|const specN =' declarations and extract their objects."""
    specs: Dict[str, str] = {}
    script_text = html.unescape(script_text)
    pattern = re.compile(r"(?:var|let|const)\s+(spec\d+)\s*=\s*{", re.I)

    for match in pattern.finditer(script_text):
        var = match.group(1)
        obj, _ = extract_js_object(script_text, match.end() - 1)
        if obj:
            specs[var] = obj
        else:
            print(f"[WARN] Could not extract object for {var}")
    return specs


def write_tsvs(specs: Dict[str, str], outdir: Path) -> List[Path]:
    """Convert extracted json to tsv."""
    outdir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    for var, js_obj in specs.items():
        try:
            data = json.loads(js_to_json(js_obj))
            values = data.get("data", {}).get("values", [])
            if not values:
                print(f"[WARN] No data.values found in {var}")
                continue

            df = pd.DataFrame(values)
            outpath = outdir / f"{var}_mqc.tsv"

            with open(outpath, "w") as f:
                f.write("# plot_type: linegraph\n")
                f.write(f"# section_name: {var}\n")
                f.write("# description: Extracted preview data\n")
                df.to_csv(f, sep="\t", index=False)

            written.append(outpath)
            print(f"[INFO] Wrote {outpath} ({len(df)} rows × {len(df.columns)} cols)")
        except Exception as e:
            print(f"[ERROR] Failed to process {var}: {e}")

    return written


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract preview data from Baysor preview HTML reports."
    )
    parser.add_argument(
        "--preview-html",
        required=True,
        help="Path to Baysor preview HTML file",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Output directory prefix (sample ID)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    input_path: Path = Path(args.preview_html)
    outdir: Path = Path(args.prefix)

    text = input_path.read_text(encoding="utf-8", errors="ignore")
    soup = BeautifulSoup(text, "html.parser")

    # get the script section
    if "<script" in text.lower():
        script_text = "\n".join(s.get_text() for s in soup.find_all("script"))
    else:
        script_text = text

    spec_variables = find_variables(script_text)
    if not spec_variables:
        print("[ERROR] No variables (spec1, spec2, spec3) found.")
        sys.exit(1)

    # write tsv files for multiqc
    written = write_tsvs(spec_variables, outdir)
    if not written:
        print("[ERROR] No TSVs written.")
        sys.exit(1)

    # get png files
    get_png_files(soup=soup, outdir=outdir)
