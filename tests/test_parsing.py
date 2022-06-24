import pytest
from epiread-tools.naming_conventions import *


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)

    def test_bedgraph_from_intervals(self):
        # genomic_intervals=["chr1:205598000-205603169"]
        # cpg_coordinates = "/Users/ireneu/PycharmProjects/dmr-cleanup/in-silico/epireads/hg19.CpG.bed.gz"
        # epiread_files = ["/Users/ireneu/PycharmProjects/dmr-cleanup/epiread_format/small_Pancreas-Beta-Z0000043H.after_fix_bug_dash.epiread.gz"]
        # outfile="."
        # header=False
        # bedfile=False
        # runner = BedgraphRunner(genomic_intervals, cpg_coordinates, epiread_files, outfile, header, bedfile)
        # runner.tobedgraph()
        pass

    def test_bedgraph_from_bedfile(self):
        pass
        #with & without header

    def test_old_epiread(self):
        pass

    def test_coord_epiread(self):
        pass

    def test_mapper(self):
        pass


if __name__ == '__main__':
    unittest.main()
