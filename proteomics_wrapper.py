#!/usr/bin/env python

"""
Wrapper script that orchestrates multiple protein localization prediction extract scripts.
Runs each extract tool with the same input files and manages output naming consistently.
"""

import argparse
import sys
import os
import re
import json
from datetime import datetime
from pathlib import Path

from cello_extract import run as run_ce
from signalp_extract import run as run_sp
#from bepipred_extract import run as run_bp
from lipop_extract import run as run_lp
from psortb_extract import run as run_pb
from tmhmm_extract import run as run_tm

#from vaxijen_extract import run as run_vj
#from virulentpred_extract import run as run_vp


# Tool definitions with their signatures
TOOLS = {
    'cello': {
        'run': run_ce,
        'params': ['input_path', 'idlist', 'output'],
        'output_prefix': 'cello_results'
    },
    'signalp': {
        'run': run_sp,
        'params': ['input_path', 'idlist', 'output'],
        'output_prefix': 'signalp5_results'
    },
    #'bepipred': {
    #    'run': run_bp,
    #    'params': ['input_path', 'sep', 'output'],
    #    'output_prefix': 'bepipred_results',
    #    'sep': ','  # Default separator for bepipred
    #},
    'lipop': {
        'run': run_lp,
        'params': ['input_path', 'idlist', 'output'],
        'output_prefix': 'lipop_results'
    },
    'psortb': {
        'run': run_pb,
        'params': ['input_path', 'idlist', 'output'],
        'output_prefix': 'psortb_results'
    },
    'tmhmm': {
        'run': run_tm,
        'params': ['input_path', 'idlist', 'output'],
        'output_prefix': 'tmhmm_results'
    },
}


def guess_input_files(input_dir: str) -> dict:
    """
    Intelligently guess input files based on naming conventions.
    Looks for common patterns: *bepipred*, *cello*, *signalp*, *lipop*, *psortb*, *tmhmm*

    Args:
        input_dir (str): Directory to search for input files

    Returns:
        dict: Mapping of tool names to input file paths
    """
    input_map = {}

    if not os.path.isdir(input_dir):
        raise ValueError(f"Input directory does not exist: {input_dir}")

    files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

    # Define patterns for each tool (case-insensitive)
    patterns = {
        'bepipred': r'(bepipred|epitope)',
        'signalp': r'(signalp|signal_p)',
        'cello': r'cello',
        'lipop': r'(lipop|lipo_p)',
        'psortb': r'(psortb|psort)',
        'tmhmm': r'(tmhmm|tm_hmm)',
    }

    for filename in files:
        for tool, pattern in patterns.items():
            if re.search(pattern, filename, re.IGNORECASE):
                if tool not in input_map:
                    input_map[tool] = os.path.join(input_dir, filename)
                    print(f"[INFO] Matched {tool}: {filename}")

    return input_map


def generate_output_name(tool_prefix: str) -> str:
    """
    Generate a standardized output filename with timestamp.
    
    Args:
        tool_prefix (str): The tool name prefix
    
    Returns:
        str: Output filename with timestamp
    """
    date = datetime.now()
    timestamp = date.strftime('%Y-%m-%d_%H-%M-%S')
    return f"{tool_prefix}_{timestamp}.tsv"


def run_tools(inputs: dict, idlist: str, output_dir: str = '.', tools_to_run: list = None):
    """
    Run specified extraction tools with provided inputs.
    Handles tools with alternative input formats by running them multiple times with unique output names.
    
    Args:
        inputs (dict): Mapping of tool names (or tool_index variants) to input file paths
        idlist (str): Path to the SeqID list file
        output_dir (str): Directory to write output files (default: current dir)
        tools_to_run (list): Specific tools to run (default: all available)
    
    Returns:
        dict: Results mapping tool names to success status
    """
    if tools_to_run is None:
        tools_to_run = list(TOOLS.keys())

    results = {}

    for tool_name in tools_to_run:
        if tool_name not in TOOLS:
            print(f"[WARN] Unknown tool: {tool_name}")
            continue

        if tool_name not in inputs:
            print(f"[SKIP] No input file found for {tool_name}")
            results[tool_name] = None
            continue

        tool_config = TOOLS[tool_name]
        input_file = inputs[tool_name]
        output_file = os.path.join(output_dir, generate_output_name(tool_config['output_prefix']))

        print(f"\n[RUN] Processing {tool_name}...")
        print(f"  Input: {input_file}")
        print(f"  Output: {output_file}")

        try:
            # Build arguments based on tool signature
#			if tool_name == 'bepipred':
#			else:
            tool_config['run'](input_file, idlist, output_file)

            if os.path.exists(output_file):
                print(f"[OK] {tool_name} completed successfully")
                results[tool_name] = output_file
            else:
                print(f"[ERROR] Output file not created for {tool_name}")
                results[tool_name] = False

        except Exception as e:
            print(f"[ERROR] {tool_name} failed: {e}")
            results[tool_name] = False

    return results


def main():
    """Main entry point for the wrapper."""
    parser = argparse.ArgumentParser(
        prog='proteomics_wrapper.py',
        description='Wrapper to run multiple extraction scripts for protein localization prediction tools.'
    )
    
    parser.add_argument('-i', '--input-dir', 
                       help='Directory containing input files. Script will intelligently match files to tools.',
                       required=True)
    parser.add_argument('-l', '--idlist',
                       help='Path to SeqID list file (single column, plain text).',
                       required=True)
    parser.add_argument('-o', '--output-dir',
                       help='Output directory for results (default: current directory).',
                       default='.')
    parser.add_argument('-t', '--tools',
                       help='Comma-separated list of tools to run (default: all). Options: ' + ', '.join(TOOLS.keys()),
                       default=None)
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    try:
        # Guess input files
        print("[INFO] Scanning for input files...")
        input_map = guess_input_files(args.input_dir)
        
        if not input_map:
            print("[ERROR] No input files matched to tools. Please check your input directory.")
            return 1
        
        print(f"[INFO] Found {len(input_map)} input file(s)")
        
        # Parse tools to run
        tools_to_run = None
        if args.tools:
            tools_to_run = [t.strip() for t in args.tools.split(',')]
        
        # Run tools
        results = run_tools(input_map, args.idlist, args.output_dir, tools_to_run)
        
        # Summary
        print("\n" + "="*60)
        print("SUMMARY")
        print("="*60)
        success_count = sum(1 for r in results.values() if r and r is not None)
        print(f"Completed: {success_count}/{len(results)}")
        for tool, result in results.items():
            if result is True or (isinstance(result, str) and os.path.exists(result)):
                print(f"  ✓ {tool}")
            elif result is False:
                print(f"  ✗ {tool}")
            else:
                print(f"  - {tool} (skipped)")
        
        return 0 if success_count > 0 else 1
    
    except Exception as e:
        print(f'[ERROR] {e}', file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
