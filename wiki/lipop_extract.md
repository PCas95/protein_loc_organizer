# `lipop_extract.py`

## Input

LipoP outputs are a variety of unstructured html pages, with no tabular structure. No file download is provided by the web server.

To create the expected input, request a LipoP analysis with **SHORT output format**, then copy the contents of the whole page and save it as a simple `.txt` file. This unformatted text will be the input for `lipop_extract.py`.

The extractor parses the unstructured text file, matches the queried SeqIDs and finally splits the text strings into single values. After this processing, the initial file is reorganised as a TSV file.

- Argument for `-i|--input` is any plain text file containing the SHORT output of LipoP. Copy-paste what you see on the html page and you're good to go: the extractor will do the rest.
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
| SeqID | prediction | score | margin | cleavage | notes |
```

LipoP can provide additional, variable information to each line. Such information is still retrieved and saved as values in the "notes" column.
