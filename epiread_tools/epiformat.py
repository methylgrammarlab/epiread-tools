from  epiread_tools.epiparser import EpireadReader, CoordsEpiread, PatReader
from epiread_tools.ModbamReader import ModbamReader

epiformat_to_reader = {"old_epiread": EpireadReader, "old_epiread_A": CoordsEpiread, "pat": PatReader, "modbam": ModbamReader}
