import glob
import load

date = '0223'
global base
base='/media/turbots/Hublot24/Share_hublot/Data'+date+'/Telephones/'

def gen_parser():    
    parser = argparse.ArgumentParser(description="Manipulate smartphone data")
    parser.add_argument('-folder', dest='folder', type=str,default=base,help='select date to process data')
#    print(parser)   
    args = parser.parse_args()
    print(args)
    return args


def main(args):
    folder = args.folder
    load.extract_all(folder)

if __name__=='__main__':
    args = gen_parser()
    main(args)
