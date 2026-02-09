#!/usr/bin/env python

"""IO helper utilities copied from proteomics_wrapper for guessing input
files and running extraction tools. Kept minimal and focused for use by
other scripts (e.g., table_builder2.py).
"""

import os
import re
from datetime import datetime
from pathlib import Path

from cello_extract import run as run_ce
from signalp_extract import run as run_sp
from lipop_extract import run as run_lp
from psortb_extract import run as run_pb
from tmhmm_extract import run as run_tm
from vaxijen_extract import run as run_vj
from virulentpred_extract import run as run_vp


# Tool definitions with their signatures (keeps compatibility with run_tools)
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
    'vaxijen': {
        'run': run_vj,
        'params': ['input_path', 'output'],
        'output_prefix': 'vaxijen_results'
    },
    'virulentpred': {
        'run': run_vp,
        'params': ['input_path', 'output'],
        'output_prefix': 'virulentpred_results'
    },
}


def guess_input_files(input_dir: str) -> dict:
    """Guess input files in a directory by common naming patterns.

    Returns a mapping: tool name -> file path
    """
    input_map = {}

    if not os.path.isdir(input_dir):
        raise ValueError(f"Input directory does not exist: {input_dir}")

    files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

    patterns = {
        'bepipred': r'(bepipred|epitope)',
        'signalp': r'(signalp|signal_p)',
        'cello': r'cello',
        'lipop': r'(lipop|lipo_p)',
        'psortb': r'(psortb|psort)',
        'tmhmm': r'(tmhmm|tm_hmm)',
        'vaxijen': r'vaxijen',
        'virulentpred': r'(virulentpred|virulen|virulen[t]?)',
    }

    for filename in files:
        for tool, pattern in patterns.items():
            if re.search(pattern, filename, re.IGNORECASE):
                if tool not in input_map:
                    input_map[tool] = os.path.join(input_dir, filename)

    return input_map


def generate_output_name(tool_prefix: str) -> str:
    """Generate a standardized output filename with timestamp."""
    date = datetime.now()
    timestamp = date.strftime('%Y-%m-%d_%H-%M-%S')
    return f"{tool_prefix}_{timestamp}.tsv"


def run_tools(inputs: dict, idlist: str, output_dir: str = '.', tools_to_run: list = None):
    """Run specified extraction tools with provided inputs.

    Returns a dict mapping tool name -> output path or status
    """
    if tools_to_run is None:
        tools_to_run = list(TOOLS.keys())

    results = {}

    for tool_name in tools_to_run:
        if tool_name not in TOOLS:
            results[tool_name] = None
            continue

        if tool_name not in inputs:
            results[tool_name] = None
            continue

        tool_config = TOOLS[tool_name]
        input_file = inputs[tool_name]
        output_file = os.path.join(output_dir, generate_output_name(tool_config['output_prefix']))

        # If no runner is implemented for this tool, skip execution
        runner = tool_config.get('run')
        if not callable(runner):
            results[tool_name] = None
            continue

        try:
            runner(input_file, idlist, output_file)

            if os.path.exists(output_file):
                results[tool_name] = output_file
            else:
                results[tool_name] = False

        except Exception:
            results[tool_name] = False

    return results
