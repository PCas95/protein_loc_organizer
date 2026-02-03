# `cello_extract.py`

## Input

CELLO's output is an unstructured html page, with a mixture of different blank space characters for spacing. Report, prediction and other relevant information are all listed together, without a tabular structure. No file download is provided by CELLO's web server.

To create the expected input, copy the contents of the whole page and save it as a simple `.txt` file. This unformatted text will be the input for `cello_extract.py`.

The extractor parses the unstructured text file and automatically removes all non-data lines and undesired characters. Prediction data and corresponding SeqIDs for each query will be retrieved, cleaned and reorganised as a TSV file.

- Argument for `-i|--input` is any plain text file containing the output of CELLO. Just copy-paste what you see on the html page and you're good to go: the extractor will do the rest.
- Argument for `-l|--idlist` is a plain text file with the protein IDs to retrieve, organised in a single column. IDs are normalised to the formats in the first two lines of the example below, so even headers from fastas are supported.
```
tr_Q8Y841_Q8Y841_LISMO
ENT65499.1
>tr|Q8Y8E3|Q8Y8E3_LISMO Lmo0961 protein OS=Listeria monocytogenes serovar 1/2a (strain ATCC BAA-679 / EGD-e) OX=169963 GN=lmo0961 PE=3 SV=1
```

## Additional arguments

- Argument for `-o|--output` is optional: by default, the output table is named automatically and saved in the current directory. This argument allows the user to provide a different path and name for the output file.

## Output

The final output of this tool is a structured TSV table, containing the columns listed below.

```
| SeqID | Protein | Cytoplasmic | Extracellular | InnerMembrane | OuterMembrane | Periplasmic | Prediction |
```

Probabilities of association to each cellular sublocalisation are stored in the middle columns.

The "Prediction" column lists all the subcellular localisations that were marked as final predictions by CELLO (*i.e.*, those with higher association scores). If more than one sublocalisation is reported by CELLO, the column will contain the full list of predictions, separated by semicolons (`;`).
