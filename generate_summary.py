#!/usr/bin/env python3
"""
Viral Screening Summary Generator

Generates a comprehensive summary of viral detection results from the pipeline.
"""

import argparse
import json
import pandas as pd
import pysam
from pathlib import Path
import sys


def parse_kraken_report(report_file):
    """Parse KrakenUniq report and extract viral hits."""
    viral_hits = {}
    
    try:
        with open(report_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    percentage = float(parts[0])
                    reads = int(parts[1])
                    taxid = parts[4]
                    name = parts[5]
                    
                    # Focus on viral hits (percentage > 0.01%)
                    if percentage > 0.01 and reads > 0:
                        viral_hits[taxid] = {
                            'name': name,
                            'reads': reads,
                            'percentage': percentage
                        }
    except FileNotFoundError:
        print(f"Warning: Kraken report not found: {report_file}")
    
    return viral_hits


def parse_bracken_output(bracken_file):
    """Parse Bracken output for abundance estimates."""
    abundance_data = {}
    
    try:
        with open(bracken_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    name = parts[0]
                    taxid = parts[1]
                    reads = int(parts[5])
                    abundance = float(parts[6])
                    
                    if reads > 0:
                        abundance_data[taxid] = {
                            'name': name,
                            'reads': reads,
                            'abundance': abundance
                        }
    except FileNotFoundError:
        print(f"Warning: Bracken output not found: {bracken_file}")
    
    return abundance_data


def analyze_viral_bam(bam_file):
    """Analyze viral BAM for coverage and depth statistics."""
    stats = {
        'total_viral_reads': 0,
        'viral_coverage': {},
        'mean_depth': 0
    }
    
    try:
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            # Count total viral reads
            stats['total_viral_reads'] = bam.count()
            
            # Analyze coverage by reference
            coverage_by_ref = {}
            for ref in bam.references:
                coverage_by_ref[ref] = {
                    'reads': 0,
                    'bases_covered': 0
                }
            
            for read in bam.fetch(until_eof=True):
                if not read.is_unmapped:
                    ref_name = bam.get_reference_name(read.reference_id)
                    if ref_name in coverage_by_ref:
                        coverage_by_ref[ref_name]['reads'] += 1
                        # Simple coverage calculation
                        coverage_by_ref[ref_name]['bases_covered'] += read.query_length
            
            stats['viral_coverage'] = coverage_by_ref
            
    except FileNotFoundError:
        print(f"Warning: Viral BAM not found: {bam_file}")
    
    return stats


def parse_fastp_qc(qc_file):
    """Parse fastp QC metrics."""
    qc_data = {}
    
    try:
        with open(qc_file, 'r') as f:
            data = json.load(f)
            qc_data = {
                'total_reads': data.get('summary', {}).get('before_filtering', {}).get('total_reads', 0),
                'total_bases': data.get('summary', {}).get('before_filtering', {}).get('total_bases', 0),
                'q20_rate': data.get('summary', {}).get('after_filtering', {}).get('q20_rate', 0),
                'q30_rate': data.get('summary', {}).get('after_filtering', {}).get('q30_rate', 0),
                'gc_content': data.get('summary', {}).get('after_filtering', {}).get('gc_content', 0)
            }
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Warning: Could not parse fastp QC: {e}")
    
    return qc_data


def generate_summary(sample_id, kraken_report, bracken_output, viral_bam, qc_json, output_file):
    """Generate comprehensive viral detection summary."""
    
    # Parse all inputs
    viral_hits = parse_kraken_report(kraken_report)
    abundance_data = parse_bracken_output(bracken_output)
    viral_stats = analyze_viral_bam(viral_bam)
    qc_data = parse_fastp_qc(qc_json)
    
    # Create summary data
    summary_rows = []
    
    # Add viral hits
    for taxid, hit_data in viral_hits.items():
        abundance = abundance_data.get(taxid, {})
        row = {
            'sample_id': sample_id,
            'taxid': taxid,
            'virus_name': hit_data['name'],
            'kraken_reads': hit_data['reads'],
            'kraken_percentage': hit_data['percentage'],
            'bracken_reads': abundance.get('reads', 0),
            'bracken_abundance': abundance.get('abundance', 0),
            'detection_method': 'kraken+bracken'
        }
        summary_rows.append(row)
    
    # Add QC metrics
    qc_row = {
        'sample_id': sample_id,
        'taxid': 'QC',
        'virus_name': 'Quality_Control',
        'total_reads': qc_data.get('total_reads', 0),
        'q20_rate': qc_data.get('q20_rate', 0),
        'q30_rate': qc_data.get('q30_rate', 0),
        'gc_content': qc_data.get('gc_content', 0),
        'viral_reads_total': viral_stats['total_viral_reads'],
        'detection_method': 'qc_metrics'
    }
    summary_rows.append(qc_row)
    
    # Create DataFrame and save
    df = pd.DataFrame(summary_rows)
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Summary generated: {output_file}")
    print(f"Found {len(viral_hits)} viral hits")
    print(f"Total viral reads: {viral_stats['total_viral_reads']}")
    
    return df


def main():
    parser = argparse.ArgumentParser(description='Generate viral screening summary')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--kraken-report', required=True, help='KrakenUniq report file')
    parser.add_argument('--bracken-output', required=True, help='Bracken output file')
    parser.add_argument('--viral-bam', required=True, help='Viral alignment BAM file')
    parser.add_argument('--qc-json', required=True, help='fastp QC JSON file')
    parser.add_argument('--output', required=True, help='Output TSV file')
    
    args = parser.parse_args()
    
    # Generate summary
    summary_df = generate_summary(
        args.sample,
        args.kraken_report,
        args.bracken_output,
        args.viral_bam,
        args.qc_json,
        args.output
    )
    
    # Print summary to stdout
    print("\n=== VIRAL DETECTION SUMMARY ===")
    viral_hits = summary_df[summary_df['detection_method'] == 'kraken+bracken']
    if not viral_hits.empty:
        print("Viral hits detected:")
        for _, row in viral_hits.iterrows():
            print(f"  {row['virus_name']}: {row['kraken_reads']} reads ({row['kraken_percentage']:.3f}%)")
    else:
        print("No viral hits detected")
    
    print("=== END SUMMARY ===")


if __name__ == '__main__':
    main() 