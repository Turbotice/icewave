import glob
import load

date = '0221'

base='/media/turbots/Hublot24/Share_hublot/Data'+date+'/Telephones/'

folders = glob.glob(base+'Bic24*')
for folder in folders:
	load.extract_all(folder)
