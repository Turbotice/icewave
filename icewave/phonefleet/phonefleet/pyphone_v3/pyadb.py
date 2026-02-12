

def get_adb_liste():
    id_list = {}
    filename = 'adb_usb_liste'
    with open(filename,'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for i,row in enumerate(spamreader):
            if i>0 and len(row)==2:
                num,id=row
                id_list[int(num)]=id
    return id_list
