#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TrioConcordance: A Comprehensive Tri-way Genotype Concordance Analyzer

This software performs a rigorous, locus-by-locus and sample-by-sample concordance 
analysis across three independent genotyping datasets. It implements two main analytical 
modules:
1. 10-Way State Classification: Classifies allelic states into 10 mutually exclusive 
   biological categories based on identical-by-state (IBS) logic.
2. Error Transition Profiling: Extracts specific mismatch patterns (e.g., Target Capture 
   Errors where WGS == Reference != Capture) and profiles the transition frequencies.
3. Standard command to extract genotype matrix from VCF
   `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' input.vcf > extracted_genotypes.txt`

Author: Xuzhen Li
Version: 1.0.0
"""

import sys
import argparse
import logging
import pandas as pd
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

def setup_logger():
    """Configure standard terminal logging."""
    logger = logging.getLogger("TrioConcordance")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

LOGGER = setup_logger()

class GenotypeAnalyzer:
    def __init__(self, wgs_path, cap_path, crr_path, out_dir, prefix):
        self.wgs_path = Path(wgs_path)
        self.cap_path = Path(cap_path)
        self.crr_path = Path(crr_path)
        self.out_dir = Path(out_dir)
        self.prefix = prefix
        self.out_dir.mkdir(parents=True, exist_ok=True)
        
    @staticmethod
    def standardize_genotype(gt_str):
        """Standardize genotype format and handle missing alleles."""
        if pd.isna(gt_str) or gt_str == '.' or str(gt_str).startswith('./.'):
            return "./."
        gt_str = str(gt_str).split('=')[-1].replace('|', '/')
        try:
            parts = gt_str.split('/')
            if len(parts) < 2: 
                return "./."
            alleles = [1 if int(p) > 0 else 0 for p in parts]
            alleles.sort()
            return f"{alleles[0]}/{alleles[1]}"
        except Exception:
            return "./."

    def load_data(self):
        """Safely load, validate, and standardize input matrices."""
        LOGGER.info("Loading genotype matrices into memory...")
        try:
            self.df_wgs = pd.read_csv(self.wgs_path, sep='\t', header=None)
            self.df_cap = pd.read_csv(self.cap_path, sep='\t', header=None)
            self.df_crr = pd.read_csv(self.crr_path, sep='\t', header=None)
        except FileNotFoundError as e:
            LOGGER.error(f"Input file missing: {e}")
            sys.exit(1)
            
        if not (self.df_wgs.shape == self.df_cap.shape == self.df_crr.shape):
            LOGGER.error("Matrix dimension mismatch. Inputs must have identical structures.")
            sys.exit(1)

        self.position_info = self.df_wgs.iloc[:, 0:2].copy()
        self.position_info.columns = ['CHROM', 'POS']

        LOGGER.info("Standardizing genotypic states across all samples...")
        def safe_map(df):
            return df.map(self.standardize_genotype) if hasattr(df, 'map') else df.applymap(self.standardize_genotype)

        self.wgs = safe_map(self.df_wgs.iloc[:, 4:])
        self.cap = safe_map(self.df_cap.iloc[:, 4:])
        self.crr = safe_map(self.df_crr.iloc[:, 4:])

        self.sample_names = [f"Sample_{i+1}" for i in range(self.wgs.shape[1])]
        self.wgs.columns = self.cap.columns = self.crr.columns = self.sample_names

    def run_10way_classification(self):
        """Execute the 10-state mutually exclusive Boolean evaluation."""
        LOGGER.info("Module 1: Computing Boolean intersections for 10-way exclusivity...")
        
        m_WC = (self.wgs == self.cap)
        m_WR = (self.wgs == self.crr)
        m_CR = (self.cap == self.crr)
        m_miss = (self.wgs == "./.")
        m_pres = (self.wgs != "./.")

        st_all_eq = m_WC & m_WR
        st_wc_only = m_WC & (~m_WR)
        st_wr_only = m_WR & (~m_WC)
        st_cr_only = m_CR & (~m_WC)
        st_all_diff = (~m_WC) & (~m_WR) & (~m_CR)

        self.res_10way = {
            'A1_Miss_All_Eq': st_all_eq & m_miss,
            'A2_Miss_WC_Only': st_wc_only & m_miss,
            'A3_Miss_WR_Only': st_wr_only & m_miss,
            'A4_Miss_CR_Only': st_cr_only & m_miss,
            'A5_Miss_All_Diff': st_all_diff & m_miss,
            'B1_Pres_All_Eq': st_all_eq & m_pres,
            'B2_Pres_WC_Only': st_wc_only & m_pres,
            'B3_Pres_WR_Only': st_wr_only & m_pres,
            'B4_Pres_CR_Only': st_cr_only & m_pres,
            'B5_Pres_All_Diff': st_all_diff & m_pres
        }

    def run_error_transition_profiling(self):
        """Extract transition frequencies for Capture Errors (WGS == CRR != Cap)."""
        LOGGER.info("Module 2: Extracting genotypic transitions for Capture Errors...")
        
        mask_capture_error = (self.wgs == self.crr) & (self.wgs != self.cap)
        results = []

        for i, sample in enumerate(self.sample_names):
            wgs_err = self.wgs.iloc[:, i][mask_capture_error.iloc[:, i]]
            cap_err = self.cap.iloc[:, i][mask_capture_error.iloc[:, i]]
            
            df_pair = pd.DataFrame({'WGS': wgs_err, 'Cap': cap_err})
            counts = df_pair.value_counts().reset_index(name='Count')
            
            for _, row in counts.iterrows():
                results.append({
                    'Sample': sample,
                    'Transition': f"{row['WGS']}->{row['Cap']}",
                    'Count': row['Count']
                })

        self.df_transitions = pd.DataFrame(results)
        if not self.df_transitions.empty:
            self.df_transitions = self.df_transitions.sort_values(by=['Transition', 'Sample'])
        
    def export_tables(self):
        """Compile and write standardized output files."""
        LOGGER.info("Generating comprehensive statistical reports...")
        
        # 1. Export 10-Way Global Stats
        out_global = self.out_dir / f"{self.prefix}_Global_Stats.tsv"
        global_stats = [{'Category': k, 'Count': v.sum().sum()} for k, v in self.res_10way.items()]
        pd.DataFrame(global_stats).to_csv(out_global, sep='\t', index=False)

        # 2. Export 10-Way Sample Stats
        out_sample = self.out_dir / f"{self.prefix}_Sample_Stats.tsv"
        sample_stats_list = []
        for s in self.sample_names:
            row = {'Sample': s}
            row.update({k: v[s].sum() for k, v in self.res_10way.items()})
            sample_stats_list.append(row)
        self.df_sample = pd.DataFrame(sample_stats_list)
        self.df_sample.to_csv(out_sample, sep='\t', index=False)
        
        # 3. Export Capture Error Transitions
        out_transition = self.out_dir / f"{self.prefix}_Capture_Error_Transitions.tsv"
        self.df_transitions.to_csv(out_transition, sep='\t', index=False)

    def export_plot(self):
        """Generate publication-ready vector graphics."""
        LOGGER.info("Rendering publication-ready figures...")
        
        mpl.rcParams.update({
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "axes.linewidth": 1.2,
            "axes.spines.top": False,
            "axes.spines.right": False
        })

        categories = list(self.res_10way.keys())
        colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', 
                  '#FFFF33', '#A65628', '#F781BF', '#999999', '#66C2A5']
        
        fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
        
        data_subset = self.df_sample[categories].values
        row_sums = data_subset.sum(axis=1)
        data_pct = (data_subset.T / row_sums).T * 100
        
        bottom = [0] * len(self.sample_names)
        for i, cat in enumerate(categories):
            ax.bar(self.sample_names, data_pct[:, i], bottom=bottom, 
                   color=colors[i], edgecolor='black', linewidth=0.5, label=cat)
            bottom = bottom + data_pct[:, i]

        ax.set_ylabel("Proportion of Loci (%)", fontweight='bold')
        ax.set_xlabel("Samples", fontweight='bold')
        ax.set_xticks(range(len(self.sample_names)))
        ax.set_xticklabels(self.sample_names, rotation=45, ha='right', fontsize=9)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=9)

        plt.tight_layout()
        out_pdf = self.out_dir / f"{self.prefix}_Concordance_Distribution.pdf"
        plt.savefig(out_pdf, transparent=True, bbox_inches='tight')
        LOGGER.info(f"Analysis successfully completed. All results saved to {self.out_dir}/")

def parse_args():
    parser = argparse.ArgumentParser(
        description="TrioConcordance: Evaluate genotype concordance and extract error transitions.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-w", "--wgs", required=True, help="Path to WGS genotypes matrix (tab-separated)")
    parser.add_argument("-c", "--cap", required=True, help="Path to Capture genotypes matrix (tab-separated)")
    parser.add_argument("-r", "--crr", required=True, help="Path to Reference genotypes matrix (tab-separated)")
    parser.add_argument("-o", "--outdir", default="Results", help="Output directory path")
    parser.add_argument("-p", "--prefix", default="Trio", help="Prefix for output files")
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    analyzer = GenotypeAnalyzer(args.wgs, args.cap, args.crr, args.outdir, args.prefix)
    
    # Core Pipeline
    analyzer.load_data()
    analyzer.run_10way_classification()
    analyzer.run_error_transition_profiling()
    analyzer.export_tables()
    analyzer.export_plot()