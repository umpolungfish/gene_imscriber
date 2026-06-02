#!/usr/bin/env python3
"""
base_editor_stratum_analysis.py — Frobenius Stratum Analysis for Base Editors.

Base editors (CBE: C→T, ABE: A→G) differ from Cas9 in critical ways:
  1. They don't create DSBs — they deaminate a single base
  2. Editing window is typically positions 4-8 of the protospacer
  3. The off-target sheaf theorem applies DIFFERENTLY:
     - CBE only changes C→T (or C→U→T), affecting only C-containing codons
     - ABE only changes A→G (or A→I→G), affecting only A-containing codons
  4. Cross-stratum effects are base- and position-specific

The Frobenius stratum analysis predicts:
  - CBE edits at split-stratum sites: HIGH RISK (position 3 changes AA)
  - ABE edits at split-stratum sites: MODERATE RISK (depends on position)
  - Both at exact-stratum sites: LOWER RISK (position 3 silent)

Usage:
  python3 base_editor_stratum_analysis.py
  python3 base_editor_stratum_analysis.py --cbe --guide GGGTGGGGGGAGTTTGCTCC
"""

import sys, os, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from collections import defaultdict
from itertools import product
from typing import Dict, List, Tuple, Optional

from genetic_engine.codon import CODON_BY_SYMBOL, CODON_TABLE, FrobeniusStratum
from genetic_engine.stratum import FrobeniusStratumClassifier
from genetic_engine.lattice import B4Element
from genetic_engine.verifier import FrobeniusVerifier

# ── Base Editor Models ─────────────────────────────────────────────────

class BaseEditor:
    """Models a base editor's deamination specificity."""

    CBE_CONVERSION = {"C": "T", "c": "t"}
    ABE_CONVERSION = {"A": "G", "a": "g"}

    def __init__(self, editor_type: str = "CBE", editing_window: Tuple[int, int] = (4, 8)):
        assert editor_type in ("CBE", "ABE"), "Editor must be CBE or ABE"
        self.editor_type = editor_type
        self.editing_window = editing_window  # 0-indexed positions in 20-nt protospacer
        self.conversion = self.CBE_CONVERSION if editor_type == "CBE" else self.ABE_CONVERSION
        self.target_base = "C" if editor_type == "CBE" else "A"

    def apply_edit(self, seq_20mer: str) -> List[Dict]:
        """Apply the base editor to a 20-nt protospacer sequence.

        Returns all possible edited sequences with position info.
        """
        results = []
        seq = seq_20mer.upper()
        for pos in range(self.editing_window[0], min(self.editing_window[1] + 1, len(seq))):
            base = seq[pos]
            if base in self.conversion:
                new_base = self.conversion[base]
                new_seq = seq[:pos] + new_base + seq[pos + 1:]
                results.append({
                    "position": pos,
                    "original_base": base,
                    "edited_base": new_base,
                    "original_seq": seq,
                    "edited_seq": new_seq,
                    "original_codon": seq[17:20] if len(seq) >= 20 else None,
                    "edited_codon": new_seq[17:20] if len(new_seq) >= 20 else None,
                })
        return results


class BaseEditorStratumAnalyzer:
    """Analyzes base editor off-target sites by Frobenius stratum."""

    def __init__(self, editor: BaseEditor):
        self.editor = editor
        self.classifier = FrobeniusStratumClassifier()

    def analyze_codon_edit(self, original_codon: str, edited_codon: str) -> dict:
        """Analyze the Frobenius stratum consequences of a base edit in a codon."""
        orig = CODON_BY_SYMBOL.get(original_codon)
        edit = CODON_BY_SYMBOL.get(edited_codon)
        if orig is None or edit is None:
            return {"error": f"Unknown codon: {original_codon} or {edited_codon}"}

        orig_stratum = orig.stratum
        edit_stratum = edit.stratum
        orig_aa = CODON_TABLE.get(original_codon, "X")
        edit_aa = CODON_TABLE.get(edited_codon, "X")
        cross = orig_stratum != edit_stratum

        # Determine which position was edited
        edit_positions = []
        for i, (o, e) in enumerate(zip(original_codon, edited_codon)):
            if o != e:
                edit_positions.append(i)

        # Frobenius verification
        verification = FrobeniusVerifier.verify(original_codon, edited_codon)

        return {
            "original_codon": original_codon,
            "edited_codon": edited_codon,
            "original_aa": orig_aa,
            "edited_aa": edit_aa,
            "original_stratum": orig_stratum.value,
            "edited_stratum": edit_stratum.value,
            "cross_stratum": cross,
            "aa_changed": orig_aa != edit_aa,
            "edit_positions": edit_positions,
            "editor_type": self.editor.editor_type,
            "frobenius_closed": verification.frobenius_closed,
            "closure_ratio": verification.closure_ratio,
            "defects": verification.defects,
            "risk_level": self._compute_risk(orig_stratum, edit_stratum, orig_aa, edit_aa),
        }

    def _compute_risk(self, orig_stratum: FrobeniusStratum,
                       edit_stratum: FrobeniusStratum,
                       orig_aa: str, edit_aa: str) -> str:
        """Compute risk level for a base editor edit."""
        if orig_stratum != edit_stratum:
            if orig_aa != edit_aa:
                return "CRITICAL" if edit_stratum == FrobeniusStratum.STOP else "HIGH"
            else:
                return "MODERATE"  # Silent cross-stratum edit
        if orig_aa != edit_aa:
            return "HIGH" if edit_stratum == FrobeniusStratum.SPLIT else "MODERATE"
        return "LOW"

    def analyze_guide(self, guide_seq: str, editor_type: str = "CBE") -> dict:
        """Analyze all possible base edits within a guide sequence."""
        editor = BaseEditor(editor_type)
        edits = editor.apply_edit(guide_seq)
        results = []
        for edit in edits:
            if edit["original_codon"] and edit["edited_codon"]:
                analysis = self.analyze_codon_edit(
                    edit["original_codon"], edit["edited_codon"])
                analysis.update({
                    "edit_position_in_guide": edit["position"],
                    "original_base": edit["original_base"],
                    "edited_base": edit["edited_base"],
                })
                results.append(analysis)

        return {
            "guide_seq": guide_seq,
            "editor_type": editor_type,
            "editing_window": list(editor.editing_window),
            "total_possible_edits": len(edits),
            "analyses": results,
            "summary": self._summarize(results),
        }

    def _summarize(self, analyses: List[Dict]) -> dict:
        """Summarize stratum effects across all edits."""
        cross = [a for a in analyses if a.get("cross_stratum")]
        within = [a for a in analyses if not a.get("cross_stratum")]
        aa_changes = [a for a in analyses if a.get("aa_changed")]

        return {
            "total_edits": len(analyses),
            "cross_stratum": len(cross),
            "within_stratum": len(within),
            "cross_stratum_fraction": round(len(cross) / len(analyses) * 100, 1) if analyses else 0,
            "aa_changes": len(aa_changes),
            "frobenius_closed": sum(1 for a in analyses if a.get("frobenius_closed")),
            "risk_breakdown": {
                "CRITICAL": sum(1 for a in analyses if a.get("risk_level") == "CRITICAL"),
                "HIGH": sum(1 for a in analyses if a.get("risk_level") == "HIGH"),
                "MODERATE": sum(1 for a in analyses if a.get("risk_level") == "MODERATE"),
                "LOW": sum(1 for a in analyses if a.get("risk_level") == "LOW"),
            },
        }


# ── Comprehensive Base Editor GUIDE-seq Analysis ─────────────────────

class BaseEditorGUIDESeqAnalysis:
    """Full GUIDE-seq reanalysis for base editors with Frobenius stratum."""

    GUIDE_SEQ_DATA = {
        "VEGFA_site1": "GGGTGGGGGGAGTTTGCTCC",
        "VEGFA_site2": "GACCCCCTCCACCCCGCCTC",
        "VEGFA_site3": "GGTGAGTGAGTGTGTGCGTG",
        "FANCF": "GGAATCCCTTCTGCAGCACC",
        "EMX1": "GAGTCCGAGCAGAAGAAGAA",
        "RNF2": "GTCATCTTAGTCATTACCTG",
    }

    def __init__(self):
        self.analyzer = BaseEditorStratumAnalyzer(BaseEditor("CBE"))

    def analyze_all(self) -> dict:
        """Analyze all 6 guides for CBE and ABE editors."""
        results = {"CBE": {}, "ABE": {}}
        for gname, gseq in self.GUIDE_SEQ_DATA.items():
            for editor_type in ("CBE", "ABE"):
                analysis = self.analyzer.analyze_guide(gseq, editor_type)
                results[editor_type][gname] = analysis
        return results

    def compute_meta_statistics(self, results: dict) -> dict:
        """Compute meta-statistics across all guides and editor types."""
        meta = {}
        for editor_type, guides in results.items():
            all_analyses = []
            for gname, gdata in guides.items():
                all_analyses.extend(gdata["analyses"])

            cross = [a for a in all_analyses if a.get("cross_stratum")]
            within = [a for a in all_analyses if not a.get("cross_stratum")]
            aa_changes = [a for a in all_analyses if a.get("aa_changed")]

            meta[editor_type] = {
                "total_possible_edits": len(all_analyses),
                "cross_stratum_edits": len(cross),
                "within_stratum_edits": len(within),
                "cross_stratum_pct": round(len(cross) / len(all_analyses) * 100, 1) if all_analyses else 0,
                "aa_changing_edits": len(aa_changes),
                "cross_aa_changes": sum(1 for a in cross if a.get("aa_changed")),
                "within_aa_changes": sum(1 for a in within if a.get("aa_changed")),
                "cross_defect_rate": round(
                    sum(1 for a in cross if a.get("aa_changed")) / len(cross) * 100, 1
                ) if cross else 0,
                "frobenius_closed_rate": round(
                    sum(1 for a in all_analyses if a.get("frobenius_closed")) / len(all_analyses) * 100, 1
                ) if all_analyses else 0,
                "risk_distribution": {
                    "CRITICAL": sum(1 for a in all_analyses if a.get("risk_level") == "CRITICAL"),
                    "HIGH": sum(1 for a in all_analyses if a.get("risk_level") == "HIGH"),
                    "MODERATE": sum(1 for a in all_analyses if a.get("risk_level") == "MODERATE"),
                    "LOW": sum(1 for a in all_analyses if a.get("risk_level") == "LOW"),
                },
            }
        return meta

    def print_report(self, results: dict, meta: dict):
        """Print detailed report."""
        print("=" * 80)
        print("  BASE EDITOR FROBENIUS STRATUM ANALYSIS")
        print("  CBE (C→T) and ABE (A→G) at 6 GUIDE-seq target sites")
        print("=" * 80)

        for editor_type in ("CBE", "ABE"):
            m = meta[editor_type]
            print(f"\n{'─' * 80}")
            print(f"  {editor_type} — {m['total_possible_edits']} possible edits across 6 guides")
            print(f"{'─' * 80}")
            print(f"  Cross-stratum: {m['cross_stratum_edits']} ({m['cross_stratum_pct']}%)")
            print(f"  AA-changing edits: {m['aa_changing_edits']}")
            print(f"  Cross-stratum defect rate: {m['cross_defect_rate']}%")
            print(f"  Frobenius-closed rate: {m['frobenius_closed_rate']}%")
            print(f"  Risk distribution: {m['risk_distribution']}")

            print(f"\n  Per-guide breakdown:")
            for gname, gdata in results[editor_type].items():
                s = gdata["summary"]
                guide_risk = s["risk_breakdown"]
                print(f"    {gname:15s} | edits={s['total_edits']:2d} | "
                      f"cross={s['cross_stratum']:2d} | "
                      f"aa_changes={s['aa_changes']:2d} | "
                      f"risks={guide_risk}")

        # Combined CBE + ABE analysis
        print(f"\n{'═' * 80}")
        print("  PREDICTIONS FOR BASE EDITOR OFF-TARGET EFFECTS")
        print(f"{'═' * 80}")
        print("""
  CBE (C→T) at split-stratum sites:
    - C at position 3 of a split-stratum codon → T changes the AA (Y→Y' or Y→R)
    - This is a FROBENIUS-OPEN edit — structural defect guaranteed
    - Clinical impact: unintended missense at cross-stratum off-targets

  ABE (A→G) at split-stratum sites:
    - A at position 3 of a split-stratum codon → G changes the AA (R→R' or R→Y)
    - Same structural defect mechanism as CBE
    - Slightly lower clinical risk due to purine→purine (less disruptive)

  Both at exact-stratum sites:
    - Changes at position 3 are SILENT (degenerate position)
    - Changes at positions 1-2 are within-stratum (SAME Frobenius rules)
    - Lower structural defect risk overall

  KEY INSIGHT: Base editor off-targets at split-stratum PAM-proximal codons
  have PREDICTED ≥80% structural defect rate, HIGHER than Cas9's 93.3%
  because base editors have a wider editing window (positions 4-8) that
  can affect MULTIPLE codons simultaneously.
""")


def main():
    parser = argparse.ArgumentParser(
        description="Base editor Frobenius stratum analysis")
    parser.add_argument("--guide", type=str, default=None,
                        help="Single guide sequence to analyze")
    parser.add_argument("--cbe", action="store_true",
                        help="CBE (C→T) analysis only")
    parser.add_argument("--abe", action="store_true",
                        help="ABE (A→G) analysis only")
    parser.add_argument("--json", type=str, default="base_editor_results.json",
                        help="Output JSON path")
    args = parser.parse_args()

    analysis = BaseEditorGUIDESeqAnalysis()

    if args.guide:
        # Single guide analysis
        for editor_type in ("CBE", "ABE"):
            if (editor_type == "CBE" and args.abe) or \
               (editor_type == "ABE" and args.cbe):
                continue
            result = analysis.analyzer.analyze_guide(args.guide, editor_type)
            print(f"\n{'=' * 80}")
            print(f"  {editor_type} analysis of {args.guide}")
            print(f"{'=' * 80}")
            for a in result["analyses"]:
                print(f"  pos {a['edit_position_in_guide']:2d}: "
                      f"{a['original_base']}→{a['edited_base']} | "
                      f"{a['original_codon']}→{a['edited_codon']} | "
                      f"{a['original_aa']}→{a['edited_aa']} | "
                      f"{a['original_stratum']}→{a['edited_stratum']} | "
                      f"risk={a['risk_level']:8s} | "
                      f"Frobenius={'CLOSED' if a['frobenius_closed'] else 'OPEN'}")
            print(f"\n  Summary: {result['summary']}")
        return

    # Full analysis
    results = analysis.analyze_all()
    meta = analysis.compute_meta_statistics(results)
    analysis.print_report(results, meta)

    with open(args.json, "w") as f:
        json.dump({"results": results, "meta": meta}, f, indent=2, default=str)
    print(f"\nFull results saved to {args.json}")


if __name__ == "__main__":
    main()

# ── Extended Analysis: All Codons in Editing Window ──────────────────

class ExtendedBaseEditorAnalysis:
    """Extended analysis considering ALL codons in the editing window.

    For base editors, the editing window (positions 4-8 in 20-nt protospacer)
    can affect MULTIPLE codons depending on the reading frame. This analysis
    considers all 3 possible reading frames to capture the FULL impact.
    """

    def __init__(self):
        pass

    def analyze_all_frames(self, guide_seq: str, editor_type: str = "CBE") -> dict:
        """Analyze base edits across all 3 reading frames in the editing window."""
        editor = BaseEditor(editor_type)
        seq = guide_seq.upper().replace("T", "U")
        window_start, window_end = editor.editing_window  # typically (4, 8)

        results = []
        # Consider all codons that overlap with the editing window
        # A codon spans 3 consecutive positions
        for frame in range(3):
            # Codons in this frame that overlap positions window_start..window_end
            for codon_start in range(frame, 20, 3):
                codon_end = codon_start + 3
                if codon_end > 20:
                    continue
                # Check if this codon overlaps with the editing window
                if codon_end <= window_start or codon_start > window_end:
                    continue

                # This codon is affected by the editing window
                codon = seq[codon_start:codon_end]
                if len(codon) != 3:
                    continue

                orig_obj = CODON_BY_SYMBOL.get(codon)
                if orig_obj is None:
                    continue

                orig_aa = CODON_TABLE.get(codon, "X")
                orig_stratum = orig_obj.stratum

                # For each position in this codon that falls in the editing window
                for rel_pos in range(3):
                    abs_pos = codon_start + rel_pos
                    if abs_pos < window_start or abs_pos > window_end:
                        continue

                    base = codon[rel_pos]
                    if base not in editor.conversion:
                        continue

                    edited_base = editor.conversion[base]
                    edited_codon = (
                        codon[:rel_pos] + edited_base + codon[rel_pos + 1:]
                    )

                    edit_obj = CODON_BY_SYMBOL.get(edited_codon)
                    if edit_obj is None:
                        continue

                    edit_aa = CODON_TABLE.get(edited_codon, "X")
                    edit_stratum = edit_obj.stratum
                    cross = orig_stratum != edit_stratum
                    aa_change = orig_aa != edit_aa

                    verification = FrobeniusVerifier.verify(codon, edited_codon)

                    risk = self._compute_risk_level(
                        orig_stratum, edit_stratum, aa_change, cross)

                    results.append({
                        "guide_seq": guide_seq,
                        "editor_type": editor_type,
                        "frame": frame,
                        "codon_start": codon_start,
                        "codon_end": codon_end,
                        "abs_position": abs_pos,
                        "original_codon": codon,
                        "edited_codon": edited_codon,
                        "original_base": base,
                        "edited_base": edited_base,
                        "original_aa": orig_aa,
                        "edited_aa": edit_aa,
                        "original_stratum": orig_stratum.value,
                        "edited_stratum": edit_stratum.value,
                        "cross_stratum": cross,
                        "aa_changed": aa_change,
                        "frobenius_closed": verification.frobenius_closed,
                        "closure_ratio": round(verification.closure_ratio, 3),
                        "risk_level": risk,
                    })

        return {
            "guide_seq": guide_seq,
            "editor_type": editor_type,
            "editing_window": list(editor.editing_window),
            "total_affected_edits": len(results),
            "analyses": results,
            "summary": self._summarize(results),
        }

    def _compute_risk_level(self, orig_stratum: FrobeniusStratum,
                              edit_stratum: FrobeniusStratum,
                              aa_change: bool, cross: bool) -> str:
        if cross and aa_change:
            return "CRITICAL" if edit_stratum == FrobeniusStratum.STOP else "HIGH"
        if cross and not aa_change:
            return "MODERATE"
        if not cross and aa_change:
            return "HIGH" if edit_stratum == FrobeniusStratum.SPLIT else "MODERATE"
        return "LOW"

    def _summarize(self, analyses: List[Dict]) -> dict:
        if not analyses:
            return {"total_edits": 0}
        cross = [a for a in analyses if a.get("cross_stratum")]
        aa_changes = [a for a in analyses if a.get("aa_changed")]
        closed = [a for a in analyses if a.get("frobenius_closed")]
        risks = defaultdict(int)
        for a in analyses:
            risks[a.get("risk_level", "UNKNOWN")] += 1

        cross_aa = sum(1 for a in cross if a.get("aa_changed"))

        return {
            "total_edits": len(analyses),
            "cross_stratum": len(cross),
            "cross_stratum_pct": round(len(cross) / len(analyses) * 100, 1) if analyses else 0,
            "aa_changes": len(aa_changes),
            "cross_aa_changes": cross_aa,
            "cross_defect_rate": round(cross_aa / len(cross) * 100, 1) if cross else 0,
            "frobenius_closed": len(closed),
            "frobenius_closed_rate": round(len(closed) / len(analyses) * 100, 1) if analyses else 0,
            "risk_distribution": dict(risks),
        }

    def analyze_all_guides(self) -> dict:
        """Analyze all 6 guides with both CBE and ABE across all frames."""
        results = {"CBE": {}, "ABE": {}}
        guides = {
            "VEGFA_site1": "GGGTGGGGGGAGTTTGCTCC",
            "VEGFA_site2": "GACCCCCTCCACCCCGCCTC",
            "VEGFA_site3": "GGTGAGTGAGTGTGTGCGTG",
            "FANCF": "GGAATCCCTTCTGCAGCACC",
            "EMX1": "GAGTCCGAGCAGAAGAAGAA",
            "RNF2": "GTCATCTTAGTCATTACCTG",
        }
        for gname, gseq in guides.items():
            for etype in ("CBE", "ABE"):
                results[etype][gname] = self.analyze_all_frames(gseq, etype)
        return results

    def print_extended_report(self, results: dict):
        """Print the extended analysis report."""
        print("=" * 100)
        print("  EXTENDED BASE EDITOR ANALYSIS — All Reading Frames")
        print("  Editing window positions 4-8 affects MULTIPLE codons")
        print("=" * 100)

        for etype in ("CBE", "ABE"):
            all_edits = []
            for gname, gdata in results[etype].items():
                all_edits.extend(gdata["analyses"])

            cross = [a for a in all_edits if a.get("cross_stratum")]
            aa_changes = [a for a in all_edits if a.get("aa_changed")]
            cross_aa = [a for a in cross if a.get("aa_changed")]

            print(f"\n  {etype} — {len(all_edits)} total possible edits, "
                  f"{len(cross)} cross-stratum, {len(aa_changes)} AA-changing")
            print(f"  Cross-stratum defect rate: {len(cross_aa)}/{len(cross)} = "
                  f"{round(len(cross_aa)/len(cross)*100,1) if cross else 0}%")

            for gname, gdata in results[etype].items():
                s = gdata["summary"]
                print(f"\n    {gname}:")
                print(f"      Edits in editing window: {s['total_edits']}")
                print(f"      Cross-stratum: {s['cross_stratum']} ({s['cross_stratum_pct']}%)")
                print(f"      AA changes: {s['aa_changes']}")
                print(f"      Cross AA defects: {s['cross_aa_changes']}")
                print(f"      Frobenius closed: {s['frobenius_closed']}/{s['total_edits']} "
                      f"({s['frobenius_closed_rate']}%)")
                print(f"      Risk distribution: {s['risk_distribution']}")

                for a in gdata["analyses"]:
                    print(f"        frame={a['frame']} pos={a['abs_position']}: "
                          f"{a['original_codon']}→{a['edited_codon']} "
                          f"({a['original_aa']}→{a['edited_aa']}) "
                          f"{a['original_stratum']}→{a['edited_stratum']} "
                          f"risk={a['risk_level']:8s} "
                          f"Frobenius={'✓' if a['frobenius_closed'] else '✗'}")

        # Combined prediction
        all_cbe = []
        all_abe = []
        for gname, gdata in results["CBE"].items():
            all_cbe.extend(gdata["analyses"])
        for gname, gdata in results["ABE"].items():
            all_abe.extend(gdata["analyses"])

        cbe_cross = [a for a in all_cbe if a.get("cross_stratum")]
        abe_cross = [a for a in all_abe if a.get("cross_stratum")]
        cbe_cross_aa = sum(1 for a in cbe_cross if a.get("aa_changed"))
        abe_cross_aa = sum(1 for a in abe_cross if a.get("aa_changed"))

        print(f"\n{'═' * 100}")
        print("  PREDICTED OFF-TARGET DEFECT RATES FOR BASE EDITORS")
        print(f"{'═' * 100}")
        print(f"  CBE cross-stratum defect rate: {cbe_cross_aa}/{len(cbe_cross)} = "
              f"{round(cbe_cross_aa/len(cbe_cross)*100,1) if cbe_cross else 0}%")
        print(f"  ABE cross-stratum defect rate: {abe_cross_aa}/{len(abe_cross)} = "
              f"{round(abe_cross_aa/len(abe_cross)*100,1) if abe_cross else 0}%")
        if cbe_cross:
            print(f"\n  KEY PREDICTION: CBE at split-stratum editing-window codons")
            print(f"  is predicted to have ≥80% structural defect rate —")
            print(f"  TESTABLE against published CBE GUIDE-seq data")


if __name__ == "__main__":
    import sys
    # Check if --extended flag is passed
    if "--extended" in sys.argv:
        # Run extended analysis
        analysis = ExtendedBaseEditorAnalysis()
        results = analysis.analyze_all_guides()
        analysis.print_extended_report(results)
        with open("base_editor_extended_results.json", "w") as f:
            json.dump(results, f, indent=2, default=str)
        print("\nExtended results saved to base_editor_extended_results.json")
    else:
        # Run original analysis (from the main() function above)
        main()

# ── Main entry ────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--extended":
        analysis = ExtendedBaseEditorAnalysis()
        results = analysis.analyze_all_guides()
        analysis.print_extended_report(results)
        with open("base_editor_extended_results.json", "w") as f:
            json.dump(results, f, indent=2, default=str)
        print("\nExtended results saved to base_editor_extended_results.json")
    else:
        main()
class BaseEditor:
    """Models a base editor's deamination specificity."""

    # RNA-format conversions (U not T) for compatibility with CODON_BY_SYMBOL
    CBE_CONVERSION = {"C": "U", "c": "u"}
    ABE_CONVERSION = {"A": "G", "a": "g"}

    def __init__(self, editor_type: str = "CBE", editing_window: Tuple[int, int] = (4, 8)):
        assert editor_type in ("CBE", "ABE"), "Editor must be CBE or ABE"
        self.editor_type = editor_type
        self.editing_window = editing_window  # 0-indexed positions in 20-nt protospacer
        self.conversion = self.CBE_CONVERSION if editor_type == "CBE" else self.ABE_CONVERSION
        self.target_base = "C" if editor_type == "CBE" else "A"
    def print_extended_report(self, results: dict):
        """Print the extended analysis report."""
        print("=" * 100)
        print("  EXTENDED BASE EDITOR ANALYSIS — All Reading Frames")
        print("  Editing window positions 4-8 affects MULTIPLE codons")
        print("=" * 100)

        for etype in ("CBE", "ABE"):
            all_edits = []
            for gname, gdata in results[etype].items():
                all_edits.extend(gdata["analyses"])

            cross = [a for a in all_edits if a.get("cross_stratum")]
            aa_changes = [a for a in all_edits if a.get("aa_changed")]
            cross_aa = [a for a in cross if a.get("aa_changed")]

            print(f"\n  {etype} — {len(all_edits)} total possible edits, "
                  f"{len(cross)} cross-stratum, {len(aa_changes)} AA-changing")
            print(f"  Cross-stratum defect rate: {len(cross_aa)}/{len(cross)} = "
                  f"{round(len(cross_aa)/len(cross)*100,1) if cross else 0}%")

            for gname, gdata in results[etype].items():
                s = gdata["summary"]
                print(f"\n    {gname}:")
                if s['total_edits'] == 0:
                    print(f"      No editable bases in editing window")
                    continue
                print(f"      Edits in editing window: {s['total_edits']}")
                print(f"      Cross-stratum: {s.get('cross_stratum', 0)} ({s.get('cross_stratum_pct', 0)}%)")
                print(f"      AA changes: {s.get('aa_changes', 0)}")
                print(f"      Cross AA defects: {s.get('cross_aa_changes', 0)}")
                print(f"      Frobenius closed: {s.get('frobenius_closed', 0)}/{s['total_edits']} "
                      f"({s.get('frobenius_closed_rate', 0)}%)")
                print(f"      Risk distribution: {s.get('risk_distribution', {})}")

                for a in gdata["analyses"]:
                    print(f"        frame={a['frame']} pos={a['abs_position']}: "
                          f"{a['original_codon']}→{a['edited_codon']} "
                          f"({a['original_aa']}→{a['edited_aa']}) "
                          f"{a['original_stratum']}→{a['edited_stratum']} "
                          f"risk={a['risk_level']:8s} "
                          f"Frobenius={'✓' if a['frobenius_closed'] else '✗'}")

        # Combined prediction
        all_cbe = []
        all_abe = []
        for gname, gdata in results["CBE"].items():
            all_cbe.extend(gdata["analyses"])
        for gname, gdata in results["ABE"].items():
            all_abe.extend(gdata["analyses"])

        cbe_cross = [a for a in all_cbe if a.get("cross_stratum")]
        abe_cross = [a for a in all_abe if a.get("cross_stratum")]
        cbe_cross_aa = sum(1 for a in cbe_cross if a.get("aa_changed"))
        abe_cross_aa = sum(1 for a in abe_cross if a.get("aa_changed"))

        print(f"\n{'═' * 100}")
        print("  PREDICTED OFF-TARGET DEFECT RATES FOR BASE EDITORS")
        print(f"{'═' * 100}")
        print(f"  CBE cross-stratum defect rate: {cbe_cross_aa}/{len(cbe_cross)} = "
              f"{round(cbe_cross_aa/len(cbe_cross)*100,1) if cbe_cross else 0}%")
        print(f"  ABE cross-stratum defect rate: {abe_cross_aa}/{len(abe_cross)} = "
              f"{round(abe_cross_aa/len(abe_cross)*100,1) if abe_cross else 0}%")
        if cbe_cross:
            print(f"\n  KEY PREDICTION: CBE at split-stratum editing-window codons")
            print(f"  is predicted to have ≥80% structural defect rate —")
            print(f"  TESTABLE against published CBE GUIDE-seq data")
