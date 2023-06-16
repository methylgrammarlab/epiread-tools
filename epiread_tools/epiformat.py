from  epiread_tools.epiparser import EpireadReader, CoordsEpiread, PatReader

epiformat_to_reader = {"old_epiread": EpireadReader, "old_epiread_A": CoordsEpiread, "pat": PatReader}
