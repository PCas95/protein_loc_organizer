# protein_loc_organizer

*Data extraction of protein localization from online prediction tools*

A suite of lightweight Python tools for extracting and standardizing output from various protein localization and property prediction servers. These tools convert non-tabular, inconsistent web-tool outputs into consistent, analysis-ready tables.

This repository aims at making available `protein_loc_organizer`, a small suite of tools used to produce data tables for the research article "I. Krasteva et al., Cold adaptation drives metabolic specialization and virulence potential in Listeria monocytogenes clonal complexes CC1 and CC9 (**PLACEHOLDER: DOI**)" (**submitted to Food Microbiology, article currently awaiting publication**).

# ⚠️ WARNING

The scripts in this repo are published for transparency of publication procedures, but are currently "beta" and not final: be defensive and do not consider them to be fully automatic yet.

## Quick Index

- [Overview](#overview)
- [Features](#features)
- [Repo Structure](#repo-structure)
- [Getting Started](#getting-started)
- [Notice / Disclaimer](#notice--disclaimer)
- [References](#references)
- [Acknowledgements](#acknowledgements)
- [Authors](#authors)
- [License Notice](#license-notice)
- [Additional Development Notes](#additional-development-notes)

## Overview

Many online protein prediction tools often present results in formats that are difficult to parse or batch-analyze.

`protein_loc_organizer` provides a unified set of extractors that:

- Parse raw text, HTML, or downloaded reports from prediction servers;
- Extract key biological features (*e.g.*, signal peptides, localization scores, transmembrane segments) reported by the prediction software;
- Convert them into standardized plain-text tables (TSV), both human- and machine-readable;
- Facilitate downstream analyses and reproducibility of data organization.

## Features

- Standardized output format;
- Standalone extractors for each prediction server/software;
- Minimal dependencies;
- Unified and consistent approach;
- Data extraction based on smart matching using a list of IDs;
- Wrapper to run multiple extractors **[BETA]**.

## Repo Structure

```
.
├── bepipred_extract.py       # Parser for BepiPred outputs
├── cello_extract.py          # Parser for CELLO subcellular localization results
├── lipop_extract.py          # Parser for LipoP predictions
├── psortb_extract.py         # Parser for PSORTb localization results
├── signalp_extract.py        # Parser for SignalP outputs
├── tmhmm_extract.py          # Parser for TMHMM transmembrane predictions
├── vaxijen_extract.py        # Parser for VaxiJen antigenicity results
├── virulentpred_extract.py   # Parser for VirulentPred outputs
├── common_utils.py           # Common utilities imported in all scripts
├── table_builder.py          # Combines parsed outputs into unified table (not automatic yet, it's undergoing revision)
├── proteomics_wrapper.py     # [NEW] Beta: Runs all parsers at once
├── LICENSE
└── README.md
```

## Getting Started

### Requirements

- Python 3.x
- Python modules:
	- `collections`
	- `datetime`
	- `itertools`
	- `operator`
	- `argparse`
	- `json`
	- `sys`
	- `os`
	- `re`
	- `pathlib`
	- `numpy`
	- `pandas`

### Installation

Just clone this repo and you are good to go:

```bash
git clone https://github.com/PCas95/protein_loc_organizer.git

cd protein_loc_organizer
```

Except for `numpy` and `pandas`, all the dependencies are standard Python modules: most likely you will only need to install those two:

```bash
pip install numpy pandas

python cello_extract.py -h
python proteomics_wrapper.py -h
```

Or, if you prefer to manage `numpy`/`pandas` in a **Conda** environment:

```bash
conda create -n protein_loc_organizer python=3.11 numpy pandas -y
conda activate protein_loc_organizer
# run any extractor inside environment:
python signalp_extract.py -h
```

> **Note:** since these tools are based on output that may be subject to change, building a fixed package distribution (conda or Docker) was deemed unnecessary and potentially misleading for this project. We may consider adding such feature in the future, if there is community interest.

### Usage

Each extractor can be run independently or it can be combined with others by using the wrapper `proteomics_wrapper.py`.

`table_builder.py` builds a final, merged table with all prediction and additional inferred classification (based on the analysis needs for the publication "I. Krasteva et al., Cold adaptation drives metabolic specialization and virulence potential in Listeria monocytogenes clonal complexes CC1 and CC9").

- **For independent runs:**

```bash
python3 protein_loc_organizer/bepipred_extract.py --input bepipred_2.csv --idlist ids.txt --output test-bepi2_results.tsv

python3 protein_loc_organizer/tmhmm_extract.py -i prediction_servers/TMHMM_ext.txt -l accession_numbers_single_col.csv
```

- **For wrapper execution:**

```bash
python3 protein_loc_organizer/proteomics_wrapper.py -i prediction_servers/ -l accession_numbers_single_col.csv -o wrapper-test/
```

> [**In-depth info in the User Guide**](./wiki/User_Guide.md)

### Supported Tools

- **BepiPred v2.0** – epitope prediction
- **CELLO** – subcellular localization
- **LipoP v1.0** – signal peptides/lipoprotein detection
- **PSORTb v3.0** – subcellular localization
- **SignalP v5.0** – signal peptide prediction
- **TMHMM v2.0** – transmembrane helices prediction
- **VaxiJen** – antigenicity prediction
- **VirulentPred v2.0** – virulence prediction

## Notice / Disclaimer

This project is **open-source** and is **not affiliated with, endorsed by, or sponsored by** the developers or maintainers of any of the prediction tools mentioned in this repository. All trademarks and tool names belong to their respective owners.

The extractors in this suite are provided solely to help users organize and standardize results for downstream analysis. They operate on outputs generated by third-party prediction servers that **may not provide standardized or tabular output formats**, and that **may not provide public code or stable programmatic interfaces**.

This project acts only as an **independent** utility to structure information from prediction servers for research workflows.

Because these prediction tools may change their underlying logic, version, webpage layout or output format/contents **without notice**, any such modification may cause the extractors in this repository to become **incompatible** with future outputs.

Therefore, users are responsible for complying with the terms of use of each prediction tool and for verifying the accuracy and continued compatibility of the extracted data.

## References

	I. Krasteva et al. "Cold adaptation drives metabolic specialization and virulence potential in Listeria monocytogenes clonal complexes CC1 and CC9" (awaiting publication).

## Acknowledgements

Many thanks to the Proteomics and Serology Unit, in particular to Ivanka Krasteva and Federica d'Onofrio, for collaboration and feedback.

## Authors

Pierluigi Castelli

## License Notice

Every script in the this repo is open source and free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The software is distributed for transparency of research methods, and in hope that it will be useful to others, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE (see ["Disclaimer" section](#notice--disclaimer)).

See the GNU Affero General Public License for more details.

## Additional Development Notes

The extractors in this repository follow a common design: they take highly unstructured, web-derived prediction results and coerce them into a consistent, machine-readable format. Because each prediction server outputs data differently (and often without a formal specification), the internal logic of each extractor must be adapted to the quirks of that tool's output.

Example: the BepiPred extractor (`bepipred_extract.py`) performs header validation, reconstructs per-position data into nested dictionaries, and applies domain-specific filtering rules (*e.g.*, probability thresholds, solvent exposure, and minimum run lengths).

Similar patterns are used throughout the suite, but each extractor necessarily implements tool-specific parsing rules depending on what information is provided by the original source.

> You can find descriptions of expected inputs in the [tools' user guide](wiki/User_Guide.md) 

### Structure and Expectations

- **Input is inherently unstructured**

	Many supported tools provide text as webpages, *ad-hoc* CSVs, or inconsistent HTML layouts. Additionally, spacing and formatting are often irregular. Because of this lack of structure, parts of the code may rely on fragile assumptions (*e.g.*, hidden formatting choices) that cannot be made more robust without official specifications.

	Starting with version **1.26.01**, the extractor scripts use an updated extraction approach that minimises dependence on unstable strings. This includes avoiding direct interaction with variable whitespace and strings known to have changed in earlier prediction outputs (strings are stripped and normalised leveraging features of Python's classes whenever possible, rather than using a fuzzy approach).

- **Extractors are intentionally lightweight**

	Each extractor is an independent script. This keeps them simple to inspect, debug, and update when upstream tools change their format. All of the extractor share the same design (core logic in a `run()` function) for easy debug, import and testing.

- **Error handling is defensive**

	Upstream tools may silently alter their output. As protection, starting from version **1.26.01**, the extractors validate headers, field counts, or expected keywords early and dinamically, rather than with hard coded strings. This helps preventing silent production of corrupted tables or errors at a string change in the output of the upstream tool. **Always check the output!**

- **Output is standardized, input never is**

	All extractors aim to emit predictable, tidy, TSV outputs. Ensuring this consistency often means compensating for missing values, renaming fields, or flattening nested patterns.

### Known Limitations

Some parts of the parsing logic could theoretically be more elegant or generalizable, but cannot be improved meaningfully because the input formats themselves have no guarantees. This cannot be changed as long as upstream tools continue to provide:

- undocumented formats
- *ad-hoc* HTML layouts
- inconsistent column structures and spacing
- hard-coded symbols or abbreviations

Starting from version **1.26.01** of the extractors, brittle or overly manual sections have been significantly reduced. Remaining areas that could be changed for more rubust functions have not been improved, since there is no reliable way to ensure correct extraction from shifting, unstructured data sources, in this specific case (see [Future Improvements](#future-improvements) below).

### Future Improvements

Coming up in the next version(s):

- Align `vaxijen_extract.py` and `virulentpred_extract.py` to new logic
- Re-factoring of `bepipred_extract.py` to comply to new, list-based logic
- `table_builder.py`:
	- Re-factoring to accept the new structured outputs and implement reliable automatic execution
	- Allow merging and production of 3 possible tables:
		- data from localisation (all tools except virulentpred and vaxijen)
		- data from virulence (only virulentpred and vaxijen)
		- final (join of the previous 2, without final score)
	- injection into wrapper as optional (or other "single-command solution") 
- Final version of `proteomics_wrapper.py`
	- more elegant error handling
	- fine-grained control over inputs
	- more elegant reports of executed and failed processes (both to `stderr` and to log file)
- Improve type hints and docstrings consistently
- ~Add tiny example inputs~ -> documentation and examples for each tool in this repo

### Contributing

If you wish to contribute extractors for additional prediction tools, you can:

1. Open an issue describing the tool you want to support, including example outputs if possible.
2. Submit a pull request with a new `*_extract.py` script following the structure and conventions used in the existing extractors.
3. Add documentation explaining the new tool and its expected output fields.

Contributors adding new extractors are encouraged to follow the general pattern:

- Validate the expected header or a minimal set of identifying fields.
- Build an internal data structure representing each protein and its per-position/per-feature metadata.
- Apply tool-specific filtering or interpretation rules.
- Produce clean CSV output compatible with table_builder.py.

Given the instability of upstream formats, **simplicity and clarity take precedence over abstraction**, at this time.

Contributions of bug fixes, improvements, and test cases are also welcome. Please feel free to discuss ideas or ask questions via the issue tracker before starting major work.

### Final Note

The main purpose of this repository is to provide access to the code used to produce data table in the manuscript "[I. Krasteva et al., Cold adaptation drives metabolic specialization and virulence potential in Listeria monocytogenes clonal complexes CC1 and CC9]" (**PLACEHOLDER: DOI**).

The scripts are released for transparency and to help the community, but the most desirable outcome would be the release of a **tested, open-source, published, reproducible and documented bioinformatics pipeline for proteomics**.
