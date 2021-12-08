#!/usr/bin/python
# coding=utf-8
# Required: Biopython, Statsmodel and txCdsPredict

from __future__ import print_function
import numpy
import time
import sys
import getopt
from sklearn.ensemble import RandomForestClassifier
import pickle
from lib.funclibs import *
######################### Functions #########################


def get_feature(lncFA, pcFA, outFILE):
    f_out = open(outFILE+".feature.txt", "w")
    get_feat(lncFA,f_out,islncs=1)
    get_feat(pcFA,f_out,islncs=0)
    f_out.close()


#############################################################
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:l:p:t:o:m:f:")
    except getopt.GetoptError:
        print("Your Type is wrong!")
        print(help_docs)
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print(help_docs)
            sys.exit()
        elif opt in ("-l",):
            lnc_fasta = arg
        elif opt in ("-p",):
            pc_fasta = arg
        elif opt in ("-o",):
            outputfile = arg
        elif opt in ("-t",):
            traintype = arg
        elif opt in ("-f",):
            predictfile = arg
        elif opt in ("-m",):
            modelfile = arg
    if args == []:
        args = [None,100,0.3,10]
    r_para = args
    
    # max_depth=None, min_samples_leaf=100, min_weight_fraction_leaf=0.3,n_estimators=10,
    
    print(opts,args)
    if traintype == "train":
        f_out = open(outputfile + ".feature.txt", "w")
        start_time = time.time()
        get_feat(lnc_fasta,f_out,islncs=1)
        get_feat(pc_fasta,f_out,islncs=0)
        f_out.close()

        f_out = open(outputfile + ".feature.txt", "w")
        model_create(outputfile, f_out, r_para)
        f_out.close()
        print("Cost Time :%ds\n" % (time.time()-start_time))
    elif traintype == "predict":
        f_out = open(outputfile + ".feature.txt", "w")
        get_feat(predictfile,f_out,islncs=-1)
        model_predict(modelfile, outputfile+".feature.txt")
        f_out.close()

if __name__ == '__main__':
    main()
