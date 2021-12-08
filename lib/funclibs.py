from Bio import SeqIO
from Bio.SeqUtils import GC,ProtParam
from Bio.Seq import Seq
import numpy
import os,re
from lib.findcds import FindCDS
from lib.fickett import Fickett
from lib.eiip import *
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import pickle
## (1) farward frame translations
## result: A translated sequence will be generated
def frame_translation(seq, shift, genetic_code=1):
    from Bio.Seq import translate
    length = len(seq) 
    frame = {} 
    fragment_length = 3 * ((length-shift) // 3) 
    frame = translate(seq[shift:shift+fragment_length], genetic_code) 
    return frame


def mRNA_translate(mRNA):
    return Seq(mRNA).translate()

def protein_param(putative_seqprot):
    return putative_seqprot.isoelectric_point()

def txcds(FAfile):
   # ### run txCdsPredict on the input FASTA file ####
    cmd = "txCdsPredict " + FAfile + " -anyStart tmp.cds"
    os.system(cmd)
    cds_len = dict()
    cds_score = dict()
    temp = open("tmp.cds")
    for line in temp:
        line_array = line.split()
        id_array = line_array[0].split('|')
        trans_id = id_array[0]
        pred_start = int(line_array[1])
        pred_end = int(line_array[2])
        pred_len = pred_end-pred_start
        cds_len[trans_id] = pred_len
        cds_score[trans_id] = float(line_array[5])
    temp.close()
    os.system("rm tmp.cds")
    return cds_len, cds_score
    ############ end of running txCdsPredict ################

def getStop(seq):
    # 5
    frame_0 = frame_translation(seq, 0, genetic_code=1)
    frame_1 = frame_translation(seq, 1, genetic_code=1)
    frame_2 = frame_translation(seq, 2, genetic_code=1)
    stop = (frame_0.count("*"), frame_1.count("*"), frame_2.count("*"))
    std_stop = numpy.std(stop)
    return std_stop

def getex(seqRNA,strinfoAmbiguous):
    '''seqCDS:longest ORF'''
    seqCDS, _ , _ , orf_fullness = FindCDS(seqRNA).longest_orf("+")
    seqprot = mRNA_translate(seqCDS)
    pep_len = len(seqprot)  # pep_len = len(seqprot.strip("*"))
    newseqprot = strinfoAmbiguous.sub("", str(seqprot))
    '''exclude ambiguous amio acid X, B, Z, J, Y in peptide sequence'''
    # 7 6
    protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot.strip("*")))
    if pep_len > 0:
        isoelectric_point = protein_param(protparam_obj)
    else:
        orf_fullness = -1
        isoelectric_point = 0.0
    return isoelectric_point, orf_fullness

def getlenperc(cds_len,first,seq_len):
    cds = cds_len[first[0]]
    len_perc = float(cds)/seq_len
    return len_perc

def getcdss(cds_score,first):
    return cds_score[first[0]]

def write_files(seq,cds_len,cds_score,first,seqRNA,strinfoAmbiguous,fickett_obj, f_out,class_type=0):
    seq_len = len(seq)  # 1
    cdsscore = getcdss(cds_score,first) # 2
    fickett_score = fickett_obj.fickett_value(seqRNA) # 3
    gc = GC(seq) # 4
    std_stop = getStop(seq) # 5
    isoelectric_point,orf_fullness = getex(seqRNA,strinfoAmbiguous) # 6 7
    len_perc = getlenperc(cds_len,first,seq_len) # 8
    cwt = get_Ewtf(get_numseq(seq, 1), wavlet)  # 9

    # print features
    print("%d" % seq_len, end='\t', file=f_out)  # 1
    print("%d" % cdsscore, end='\t', file=f_out)  # 2
    print("%s" % str(fickett_score), end='\t', file=f_out)  # 3
    print("%0.2f" % gc, end='\t', file=f_out)  # 4
    print("%.6f" % std_stop, end='\t', file=f_out)  # 5
    print("%s" % str(orf_fullness), end='\t', file=f_out)  # 6
    print("%s" % str(isoelectric_point), end='\t', file=f_out)  # 7
    print("%0.2f" % len_perc, end='\t', file=f_out)  # 8
    for o in cwt:
        print("%0.8f" % o, end='\t', file=f_out)  # 9
    if class_type != -1:
        print("%d" % class_type, file=f_out)


def get_feat(FAfile,outFileLink,islncs = 0):
    strinfoAmbiguous = re.compile("X|B|Z|J|U", re.I)
    ptU = re.compile("U", re.I)
    fickett_obj = Fickett()
    # ------
    
    f_in = open(FAfile)
    cds_len, cds_score = txcds(FAfile) # 2 # ### run txCdsPredict on the input FASTA file ####
    print("Extra Informations Start:")
    proIndex = 0
    for record in SeqIO.parse(f_in, "fasta"):
        proIndex = proIndex + 1
        print("%s\r" % str(proIndex).zfill(7),end='')
        ID = record.id
        first = ID.split("|")
        seq = record.seq
        seqRNA = ptU.sub("T", str(seq).strip())
        seqRNA = seqRNA.upper()
        write_files(seq,cds_len,cds_score,first,seqRNA,strinfoAmbiguous,fickett_obj, outFileLink,islncs)
    print("".zfill(7))
    print("Extra Informations Finished")
    f_in.close()

help_docs = '''
Usage: 
1 . \033[31;1mIf you want to train your model\033[0m
    \033[31;1mExample:\033[0m prelnc2.py \033[31;1m-l\033[0m ath_lnc.fa \033[31;1m-p\033[0m ath_pc.fa \033[31;1m-t\033[0m train \033[31;1m-o\033[0m ath_model.pkl \033[31;1m0 100 0.3 10\033[0m
    -l lncRNAfile    -p pcRNAfile    -train for trainning   -o outputfile
                        \033[31;1m 0 100 0.3 10 \033[0m        
    Default max_depth=\033[31;1mNone\033[0m, min_samples_leaf=\033[31;1m100\033[0m, min_weight_fraction_leaf=\033[31;1m0.3\033[0m,n_estimators=\033[31;1m10\033[0m,
2 . \033[31;1mIf you want to predict your RNA\033[0m
    \033[31;1mExample:\033[0mprelnc2.py \033[31;1m-f\033[0m ath_lnc.fa \033[31;1m-m\033[0m ath_model.pkl \033[31;1m-t\033[0m predict \033[31;1m-o\033[0m result.txt
    -f  file_in for predict  -m modelfile   -predict   -o outputfile
'''

def model_create(model_path, train_f, r_para):
    # create model
    data = pd.read_table(train_f, header=None)
    x = numpy.array(data[:-1])
    y = numpy.array(data[-1])
    rf = RandomForestClassifier(n_jobs=-1,
                                max_depth=r_para[0],	 
                                min_samples_leaf=r_para[1],	 
                                min_weight_fraction_leaf=r_para[2],
                                n_estimators = r_para[3]
                                # max_depth=None, min_samples_leaf=100, min_weight_fraction_leaf=0.3,n_estimators=10,
                                )

    model = rf.fit(x, y)
    with open(model_path, 'wb') as f:
        pickle.dump(model, f)


def model_predict(model_path, test_f):
    with open(model_path, 'rb') as m:
        model = pickle.load(m)

    # predict model
    data = pd.read_table(test_f, header=None)
    x = numpy.array(data)
    y_pre = model.predict(x)
    with open('temp', 'w') as t:
        numpy.savetxt(t, y_pre, fmt='%0.8f')
    cl = "paste " + test_f + " temp > " + test_f + ".predict.txt"
    os.system(cl)
    os.system("rm temp")
