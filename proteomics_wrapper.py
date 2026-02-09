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

from io_helpers import guess_input_files, generate_output_name, run_tools, TOOLS


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
