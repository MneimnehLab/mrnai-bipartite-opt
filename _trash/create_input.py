#!/usr/bin/python

def main():
    print 'mrnai-opt requires input files in a particular format.\nCreate input files for your sequences here.\n'
    
    num_rnas = raw_input('How many RNAs do you have? ')

    try:
        num_rnas = int(num_rnas)
        sequences = []
        for i in range(num_rnas):
            rna_name = raw_input('RNA Name: ')
            sequence = raw_input('RNA Sequence: ')
            even_odd = int(raw_input('Even (0) or Odd (1): '))

            



    except:
        print 'Invalid input'


if __name__ == "__main__":
    main()
