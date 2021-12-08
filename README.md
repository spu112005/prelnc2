Usage: 
1 . If you want to train your model

    Example: prelnc2.py -l ath_lnc.fa -p ath_pc.fa -t train -o ath_model.pkl 0 100 0.3 10
    
    -l lncRNAfile    -p pcRNAfile    -train for trainning   -o outputfile
    
                         0 100 0.3 10         
                         
    Default max_depth=None, min_samples_leaf=100, min_weight_fraction_leaf=0.3,n_estimators=10,
    
2 . If you want to predict your RNA

    Example:prelnc2.py -f ath_lnc.fa -m ath_model.pkl -t predict -o result.txt
    
    -f  file_in for predict  -m modelfile   -predict   -o outputfile
    
