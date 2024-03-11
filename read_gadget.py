import bigfile
from fake_spectra import spectra
bf = bigfile.BigFile(“PART 001”) pos = bf[“1/Position”][:]
rates = GasProperties(redshift, snap, hubble, units=unit)
temp = rates.get_temp(0, -1)
dens = rates.get_code_rhoH(0, -1)
