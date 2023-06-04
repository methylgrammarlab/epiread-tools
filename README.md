# Package: epiread-tools

Utils and conventions for using biscuit epiread files.


## Usage

basic install
```
pip install git+https://github.com/methylgrammarlab/epiread-tools
```
alternatively, you may need to use an [access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
```
pip install git+https://<PERSONAL ACCESS TOKEN>@github.com/methylgrammarlab/epiread-tools.git
```

you could also clone the repository and manually install
```
git clone https://<PERSONAL ACCESS TOKEN>@github.com/methylgrammarlab/epiread-tools.git
cd epiread-tools
python setup.py install
```

#### Configuration Parameters

The following parameters can be specified in the `config` dictionary or overwritten with specific arguments:

- `cpg_coordinates` (str): Path to the tab-delimited, gzipped file containing the coordinates of CpG sites.
- `genomic_intervals` (str): Intervals to process, either comma-delimited chrN:start-stop or a bed file.
- `bedfile` (bool): Set to True if the genomic intervals are in a bed file.
- `header` (bool): Set to True if the genomic intervals file has a header (only relevant if `bedfile` is True).
- `epiread_files` (list): List of gzipped epiread files.
- `epiformat` (str): Format of the epiread files. Currently, only supports `old_epiread`, `old_epiread_A`, and `pat`.
- `outfile` (str): Output file path for the generated bedgraph file.
- `minimal_cpg_per_read` (int): Minimum number of CpGs required per read. Default is 1.

All parameters can be specified with a `config.json` file, and specific arguments can overwrite the JSON file.


#### convert epiread to bedgraph
```
epireadToBedgraph -j <config.json>
```
any parameter in the config file can be overwritten via the command line, e.g.:

```
epireadToBedgraph -j <config.json> --outfile=<output_file> -i chr3:50-150,ch10:440-450
```

To parse epiread with cpg coordinates use flag -A

if for some reason that doesn't work (Command not found), just find the script with
```
pip show epiread_tools
```
and run
```
python3 epireadToBedgraph.py -j <config.json>
```


