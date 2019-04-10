'''				-- Jose Santiago Sanchez Fragoso --

    proAmps   - Antimicrobial propetide cleavage site predictor and Antimicrobial peptides classificator

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

import sys, os, getopt,textwrap

libpath = os.path.dirname(os.path.realpath(sys.argv[0]))
def main(argv):
    try:
       opts, args = getopt.getopt(argv,"hcpi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
       print_help()
       sys.exit(2)
    flag = 0
    flagO = 0
    flagG = 0
    for opt, arg in opts:
       if opt == "-h":
           print_help()
           sys.exit()
       elif opt == "-i":
          flagG = 1
          fasta = arg
       elif opt == "-c":
           flag = "c"
       elif opt == "-p":
           flag = "p"
       elif opt == "-o":
           flagO = 1
           outfile = arg
    if flagG == 0:
        print_help()
        sys.exit()

    import proAmps_module
    if flag == "c":
        cs = proAmps_module.CSpredictor(fasta, libpath + "/CS_N_dnn.h5",libpath + "/CS_N_svm.pkl" , libpath+ "/CS_N_svm.pkl",  libpath + "/CS_C_dnn.h5", libpath + "/CS_C_rf.pkl", libpath + "/CS_C_svm.pkl")         
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
        ap = proAmps_module.AMPpredictor(fasta, libpath + "/multi_rf.pkl" ,libpath + "/proAmps.h5",libpath + "/proAmps_rf.pkl")
        ap.pred_from_propeps( libpath + "/CS_N_dnn.h5",libpath + "/CS_N_svm.pkl" , libpath+ "/CS_N_svm.pkl",  libpath + "/CS_C_dnn.h5", libpath + "/CS_C_rf.pkl", libpath + "/CS_C_svm.pkl")
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
        ap = proAmps_module.AMPpredictor(fasta,libpath + "/multi_rf.pkl",libpath + "/proAmps.h5",libpath + "/proAmps_rf.pkl")
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



def print_help():
        print(textwrap.dedent("""\    proAmps   - Antimicrobial propetide cleavage site predictor and Antimicrobial peptides classificator

        USAGE: 
            proAmps.py -i <FASTA> [OPTIONS]

        OPTIONS:

            -h ; help
            -i <INPUT FASTA>; mandatory
            -c ; predict mature sequeces with out AMP prediction
            -p ; predict mature sequences and predict AMP with mature sequences; 
                    if not specified prediction of AMP id done assuming mature peptides are given
            -o <OUTPUT FILE>; indicate the name and path of the output fasta qith the anotation given by proAmps, DEFAULT: results are directed to STDOUT 

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
     
            #Predict AMPs with no mature sequence prediction, direct results to file
            #(no AMPs will be predicted as they are not mature sequences) 
            python proAmps.py -i test.fasta -o ./annotated_sequences.fasta 
            
            #Predict AMPs predictiong also mature sequences AMP
            #print result to STDOUT
            python proAmps.py -i test.fasta -p

		"""))

if __name__ == "__main__":
   main(sys.argv[1:])
