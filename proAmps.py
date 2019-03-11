'''				-- Jose Santiago Sanchez Fragoso --

    proAmps   - Antiicrobial propetide cleavage site predictor and Antimicrobial peptides classificator

        USAGE: 
            
            proAmps.py -i <FASTA> [OPTIONS]

        OPTIONS:
            -i <INPUT FASTA>; mandatory
            -c ; predict mature sequeces with out AMP prediction
            -p ; predict mature sequences and predict AMP with mature sequences; 
                    if not specified prediction of AMP id done assuming mature peptides are given
            -o <OUTPUT FILE>; indicate the name and path of the output fasta qith the anotation given by proAmps, if not specified results are directed to STDOUT 

        DEPENDENCIES:
            Keras 2.2 (Python)
            SciKitlearn (Phython)
            Pandas (Python)
            Numpy (Python)
            importr library (Python)
            R       v > 3.0
            Peptides package (R)
      
        OTHERS:
            proAmps.py, proAmpsLib.py, N_ter_proAmp.h5, and proAmp.h5 must be in the same directory when executed.

        EXAMPLES:
            #Predict mature AMPs and direct results to file 
            python proAmps.py -i test.fasta -c -o ./mature_sequences.fasta
     
            #Predict AMPs with previous mature AMP prediction, print result to STDOUT
            python proAmps.py -i test.fasta -p

            #Predict AMPs with no mature sequence prediction, direct results to file
            python proAmps.py -i test.fasta -o ./annotated_sequences.fasta 
'''

import proAmps_module
import sys, os, getopt

libpath = os.path.dirname(os.path.realpath(sys.argv[0]))
def main(argv):
    try:
       opts, args = getopt.getopt(argv,"hcpi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
       print('proAmps_module.py -i <FASTA> [OPTIONS]')
       sys.exit(2)
    flag = 0
    flagO = 0
    flagG = 0
    for opt, arg in opts:
       if opt == "-i":
          flagG = 1
          fasta = arg
       elif opt == '-h':
          print('proAmps_module.py -i <FASTA> [OPTIONS]')
          sys.exit()
       elif opt == "-c":
           cs = proAmps_module.CSpredictor(fasta,libpath + "/N_ter_proAmp.h5")         
           flag = "c"
       elif opt == "-p":
           ap = proAmps_module.AMPpredictor(fasta,libpath + "/proAmp.h5")
           ap.pred_from_propeps(libpath + "/N_ter_proAmp.h5")
           flag = "p"
       elif opt == "-o":
           flagO = 1
           outfile = arg
    if flagG == 0:
        print('proAmps_module.py -i <FASTA> [OPTIONS]')
        sys.exit()
    if flag == "c":
        if flagO == 1:
            tem = sys.stdout
            sys.stdout = f = open(outfile, 'w')
            cs.print_fasta()
            sys.stdout = tem
            f.close()
            sys.exit()
        else:
            cs.print_fasta()
            sys.exit()
    elif flag == "p":
        if flagO == 1:
            tem = sys.stdout
            sys.stdout = f = open(outfile, 'w')
            ap.print_fasta()
            sys.stdout = tem
            f.close()
            sys.exit()
        else:
            ap.print_fasta() 
            sys.exit()
    elif flag == 0:
        ap = proAmps_module.AMPpredictor(fasta,libpath + "/proAmp.h5")
        ap.predict_amps()
        if flagO == 1:
            tem = sys.stdout
            sys.stdout = f = open(outfile, 'w')
            ap.print_fasta()
            sys.stdout = tem
            f.close()
            sys.exit()
        else:
            ap.print_fasta()
            sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])
