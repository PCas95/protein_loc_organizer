# `bepipred_extract.py`

## Input

BepiPred2 already outputs a `.csv` file. The extractor's job is to identify Epitopes Exposed To Solvent (EETS). To do this, the input table is parsed and, for each protein entry, the parser retrieves lines for aminoacid positions that meet the following requirements:

1) epitope probability >= 0.5;
2) exposed/buried == E;
3) 4 or more consecutive aminoacids that meet conditions 1 and 2.

> Requirements were defined in the paper. Please read the main README of this repo for the link.

- Argument for `-i|--input` is **the `.csv` table produced by BepiPred2**.
- Argument for `-l|--idlist` is a plain text file with the protein IDs to retrieve, organised as a single column. IDs are normalised (first two lines in the example below), so even headers from fastas are supported (example from line 3).
```
tr_Q8Y841_Q8Y841_LISMO
ENT65499.1
>tr|Q8Y8E3|Q8Y8E3_LISMO Lmo0961 protein OS=Listeria monocytogenes serovar 1/2a (strain ATCC BAA-679 / EGD-e) OX=169963 GN=lmo0961 PE=3 SV=1
```

> **Note:** this tool has been used and tested on output from BepiPred **version 2**. At the time of writing, BepiPred3's `.csv` output currently doesn't have information on aminoacid exposition. For this reason, input from that version will not be accepted. Input from more recent versions of the software may become supported in the future, if there will be interest and data availability to do it. 

## Additional arguments

- Argument for `-s|--separator` is optional: default is comma (`,`), since that's the default separator of BepiPred's output. The user can provide a separator to override the default one, if the table is from older or alternate outputs, or has already been processed. This does not change the separator used for the output (always `\t`).
- Argument for `-o|--output` is also optional: by default, the output table is named automatically and saved in the current directory. This argument allows the user to provide a different path and name for the output file.

## Output

The final output of this tool is a filtered TSV table.

The input is filtered according to the criteria expressed above, then restructured and saved as a table with `\t` separator. Columns and value formats are the same as the original table used as input.
