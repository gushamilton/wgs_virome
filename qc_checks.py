#!/usr/bin/env python3
"""
Viral Screening QC Checks

Performs quality control checks on viral screening pipeline results.
"""

import argparse
import json
import sys
from pathlib import Path


def check_read_counts(candidate_reads, qc_reads, nonhuman_reads):
    """Check read count consistency and flag potential issues."""
    qc_flags = []
    qc_metrics = {}
    
    # Calculate ratios
    qc_ratio = qc_reads / candidate_reads if candidate_reads > 0 else 0
    nonhuman_ratio = nonhuman_reads / qc_reads if qc_reads > 0 else 0
    
    qc_metrics['candidate_reads'] = candidate_reads
    qc_metrics['qc_passed_reads'] = qc_reads
    qc_metrics['nonhuman_reads'] = nonhuman_reads
    qc_metrics['qc_pass_rate'] = qc_ratio
    qc_metrics['nonhuman_rate'] = nonhuman_ratio
    
    # QC checks
    if candidate_reads < 5000:
        qc_flags.append('LOW_CANDIDATE_READS')
        qc_metrics['low_coverage_flag'] = True
    else:
        qc_metrics['low_coverage_flag'] = False
    
    if qc_ratio < 0.8:
        qc_flags.append('LOW_QC_PASS_RATE')
        qc_metrics['low_qc_flag'] = True
    else:
        qc_metrics['low_qc_flag'] = False
    
    if nonhuman_ratio > 0.95:
        qc_flags.append('HIGH_NONHUMAN_RATE')
        qc_metrics['high_nonhuman_flag'] = True
    else:
        qc_metrics['high_nonhuman_flag'] = False
    
    return qc_flags, qc_metrics


def check_fastp_qc(qc_json_file):
    """Check fastp QC metrics."""
    qc_flags = []
    qc_metrics = {}
    
    try:
        with open(qc_json_file, 'r') as f:
            data = json.load(f)
            
        # Extract metrics
        before_filtering = data.get('summary', {}).get('before_filtering', {})
        after_filtering = data.get('summary', {}).get('after_filtering', {})
        
        qc_metrics['total_reads_before'] = before_filtering.get('total_reads', 0)
        qc_metrics['total_bases_before'] = before_filtering.get('total_bases', 0)
        qc_metrics['total_reads_after'] = after_filtering.get('total_reads', 0)
        qc_metrics['total_bases_after'] = after_filtering.get('total_bases', 0)
        qc_metrics['q20_rate'] = after_filtering.get('q20_rate', 0)
        qc_metrics['q30_rate'] = after_filtering.get('q30_rate', 0)
        qc_metrics['gc_content'] = after_filtering.get('gc_content', 0)
        
        # QC checks
        if qc_metrics['q20_rate'] < 0.8:
            qc_flags.append('LOW_Q20_RATE')
            qc_metrics['low_q20_flag'] = True
        else:
            qc_metrics['low_q20_flag'] = False
        
        if qc_metrics['q30_rate'] < 0.6:
            qc_flags.append('LOW_Q30_RATE')
            qc_metrics['low_q30_flag'] = True
        else:
            qc_metrics['low_q30_flag'] = False
        
        if qc_metrics['gc_content'] < 0.3 or qc_metrics['gc_content'] > 0.7:
            qc_flags.append('ATYPICAL_GC_CONTENT')
            qc_metrics['atypical_gc_flag'] = True
        else:
            qc_metrics['atypical_gc_flag'] = False
            
    except (FileNotFoundError, json.JSONDecodeError) as e:
        qc_flags.append('FASTP_QC_PARSE_ERROR')
        qc_metrics['fastp_error'] = str(e)
    
    return qc_flags, qc_metrics


def check_contamination(kraken_report_file):
    """Check for common contaminants."""
    qc_flags = []
    qc_metrics = {}
    
    try:
        with open(kraken_report_file, 'r') as f:
            phiX_found = False
            phiX_percentage = 0.0
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    percentage = float(parts[0])
                    name = parts[5].lower()
                    
                    # Check for PhiX (common lab contaminant)
                    if 'phix' in name or 'phi-x174' in name:
                        phiX_found = True
                        phiX_percentage = percentage
                        break
            
            qc_metrics['phix_detected'] = phiX_found
            qc_metrics['phix_percentage'] = phiX_percentage
            
            if phiX_percentage > 0.1:
                qc_flags.append('HIGH_PHIX_CONTAMINATION')
                qc_metrics['high_phix_flag'] = True
            else:
                qc_metrics['high_phix_flag'] = False
                
    except FileNotFoundError:
        qc_flags.append('KRAKEN_REPORT_NOT_FOUND')
        qc_metrics['kraken_error'] = 'File not found'
    
    return qc_flags, qc_metrics


def generate_qc_summary(sample_id, candidate_reads, qc_reads, nonhuman_reads, qc_json_file, kraken_report_file):
    """Generate comprehensive QC summary."""
    
    # Perform all QC checks
    read_flags, read_metrics = check_read_counts(candidate_reads, qc_reads, nonhuman_reads)
    fastp_flags, fastp_metrics = check_fastp_qc(qc_json_file)
    contam_flags, contam_metrics = check_contamination(kraken_report_file)
    
    # Combine all flags and metrics
    all_flags = read_flags + fastp_flags + contam_flags
    all_metrics = {**read_metrics, **fastp_metrics, **contam_metrics}
    
    # Determine overall QC status
    if not all_flags:
        qc_status = 'PASS'
    elif len(all_flags) <= 2:
        qc_status = 'WARNING'
    else:
        qc_status = 'FAIL'
    
    # Create QC summary
    qc_summary = {
        'sample_id': sample_id,
        'qc_status': qc_status,
        'qc_flags': all_flags,
        'qc_metrics': all_metrics,
        'timestamp': str(Path().cwd() / 'timestamp.txt')  # Placeholder for actual timestamp
    }
    
    return qc_summary


def main():
    parser = argparse.ArgumentParser(description='Perform QC checks on viral screening results')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--candidate-reads', type=int, required=True, help='Number of candidate reads extracted')
    parser.add_argument('--qc-reads', type=int, required=True, help='Number of reads passing QC')
    parser.add_argument('--nonhuman-reads', type=int, required=True, help='Number of non-human reads')
    parser.add_argument('--qc-json', required=True, help='fastp QC JSON file')
    parser.add_argument('--kraken-report', required=True, help='KrakenUniq report file')
    parser.add_argument('--output', required=True, help='Output JSON file')
    
    args = parser.parse_args()
    
    # Generate QC summary
    qc_summary = generate_qc_summary(
        args.sample,
        args.candidate_reads,
        args.qc_reads,
        args.nonhuman_reads,
        args.qc_json,
        args.kraken_report
    )
    
    # Save QC summary
    with open(args.output, 'w') as f:
        json.dump(qc_summary, f, indent=2)
    
    # Print QC summary to stdout
    print(f"\n=== QC SUMMARY FOR {args.sample} ===")
    print(f"QC Status: {qc_summary['qc_status']}")
    
    if qc_summary['qc_flags']:
        print("QC Flags:")
        for flag in qc_summary['qc_flags']:
            print(f"  - {flag}")
    else:
        print("No QC flags raised")
    
    print(f"Key Metrics:")
    print(f"  - Candidate reads: {qc_summary['qc_metrics']['candidate_reads']:,}")
    print(f"  - QC pass rate: {qc_summary['qc_metrics']['qc_pass_rate']:.3f}")
    print(f"  - Non-human rate: {qc_summary['qc_metrics']['nonhuman_rate']:.3f}")
    print(f"  - Q20 rate: {qc_summary['qc_metrics'].get('q20_rate', 'N/A'):.3f}")
    print(f"  - Q30 rate: {qc_summary['qc_metrics'].get('q30_rate', 'N/A'):.3f}")
    
    if qc_summary['qc_metrics'].get('phix_detected', False):
        print(f"  - PhiX contamination: {qc_summary['qc_metrics']['phix_percentage']:.3f}%")
    
    print("=== END QC SUMMARY ===")
    
    # Exit with appropriate code
    if qc_summary['qc_status'] == 'FAIL':
        sys.exit(1)
    elif qc_summary['qc_status'] == 'WARNING':
        sys.exit(2)
    else:
        sys.exit(0)


if __name__ == '__main__':
    main() 