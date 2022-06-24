# Package: epiread-tools

Utils and conventions for using biscuit epiread files.


## Usage

basic install
```
pip install git+https://github.com/methylgrammarlab/epiread-tools
```
alternatively, you may need you use an [access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
```
pip install git+https://<PERSONAL ACCESS TOKEN>@github.com/methylgrammarlab/epiread-tools.git
```

convert epiread to bedgraph
```
python3 epireadToBedgraph.py <cpg_coordinates_file> <epiread> <run_name> <output_file> -i < regions, e.g.: chr1:205500000-205700000,chr10:17400000-17600000>
```

To parse epiread with cpg coordinates use flag -A

