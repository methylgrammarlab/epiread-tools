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

To parse epiread files:
The config.json file contains the following parameters:
* "cpg_coordinates": tab delimited, gzipped file with coordinates of CpG sites
* "genomic_intervals": intervals to process, either comma delimited or a bed file
* "bedfile": set True if genomic_intervals are in a bed file
* "header": set True if genomic_intervals file has a header (only relevant if "bedfile" is True)
* "epiread_files": list of epiread files, each file is gzipped 
* "epiformat": format of epiread file, currently only old_epiread, old_epiread_A and pat
* "outfile": output file path
* "minimal_cpg_per_read": default 1, only process reads with at least n CpGs


convert epiread to bedgraph
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


