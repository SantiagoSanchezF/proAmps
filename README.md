# proAmps  
Antiicrobial propetide cleavage site predictor and Antimicrobial peptides classificator   
     
Dependencies for instalation:   
   linux/OSX 64X   
   Anaconda/Miniconda (https://conda.io/docs/user-guide/install/index.html), for enviroment instalation   
     
     
INSTALLATION:   
   from this folder:   
        conda env create -f proAmps.yml   
        source activate proAmps   
        Rscript install.peptides.R   
            
USAGE:   
   activate the proAmps enviroment:   
        $ source activate proAmps   
   from this folder execute:   
        $ python proAmps -i test.fasta   
             
EXAMPLES:   
   #Predict mature AMPs and direct results to file    
   python proAmps.py -i test.fasta -c -o ./mature_sequences.fasta   
      
   #Predict AMPs with previous mature AMP prediction, print result to STDOUT    
   python proAmps.py -i test.fasta -p    
    
   #Predict AMPs with no mature sequence prediction, direct results to file    
   python proAmps.py -i test.fasta -o ./annotated_sequences.fasta    
    
