# Package: epiread-tools

Utils and conventions for using biscuit epiread files.


## Usage

### basic install
'''
pip install git+https://github.com/methylgrammarlab/epiread-tools
'''

### convert epiread to bedgraph
'''
python3 epireadToBedgraph.py <cpg_coordinates_file> <epiread> <run_name> <output_file> -i < regions, e.g.: chr1:205500000-205700000,chr10:17400000-17600000>
'''