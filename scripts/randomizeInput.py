from argparse import ArgumentParser
from utils.Protein import Protein
import os


def processProtein(inputPath, outputPath, clean=True):
    p = Protein()
    p.fromFile(inputPath)
    p.prepare()
    p.makeRandomRotation(inplace=True)
    p.makeRandomTranslation(inplace=True)
    p.toFile(outputPath)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', required=True, help='PDB file or folder containing PDB files to be randomized')
    parser.add_argument('-o', required=True, help='Folder or file where randomized structures are being saved')
    parser.add_argument('-clean', default=True, type=bool, help='Whether to include a cleaning step removing ligands and adding H-Atoms')
    args = parser.parse_args()


    if os.path.isdir(args.i):
        inputPath = args.i if args.i.endswith('/') else '{}/'.format(args.i)
        outputPath = args.o if args.i.endswith('/') else '{}/'.format(args.o)
        for f in os.listdir(inputPath):
            print('Processing: ', '{}{}'.format(inputPath, f))
            fileType = os.path.splitext(f)[-1]
            if fileType != '.pdb':
                continue

            processProtein('{}{}'.format(inputPath, f), '{}{}'.format(outputPath, f), clean=args.clean)
            print('Saved to: ', '{}{}'.format(outputPath, f))

    else:
        print('Processing: ', args.i)
        processProtein(args.i, args.o, clean=args.clean)
        print('Saved to: ', args.o)