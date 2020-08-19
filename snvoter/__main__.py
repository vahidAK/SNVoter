#! /usr/bin/env python3
# coding=utf-8

# Copyright (C) 2020  Vahid Akbari

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
SNVoter: A top up tool to enhance SNV calling from Nanopore sequencing data.
"""

__author__ = "Vahid Akbari"
__email__ = "vakbari@bcgsc.ca"
__copyright__ = "Copyright (C) 2020, " + __author__
__license__ = "GPLv3"
__collaborator__ = "Jean-Michel Garant"

import os
os.environ['PYTHONHASHSEED'] = '0'
import numpy as np
import random as rn
import tensorflow as tf
import warnings

rn.seed(1)
np.random.seed(1)
tf.random.set_seed(1)

#The tensorflow version is 2.1.0
#The scikit-learn version is 0.21.3
import sys
import multiprocessing as mp
import argparse
import os
import pysam
from matplotlib import pyplot as plt
from tqdm import tqdm
from itertools import repeat
from pickle import dump, load
import statistics
import pandas as pd
import gzip
import bz2
from collections import defaultdict
try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_curve
    from sklearn.metrics import auc, confusion_matrix
except:
    raise ImportError("It seems that sklearn has not been installed. "
                      "Please insall it first.")
try:
    import keras
    from keras.callbacks import ModelCheckpoint
    from keras.models import Sequential, load_model
    from keras.layers import Dense
    from keras import backend as K
except:
    raise ImportError("It seems that keras and/or tensorflow have/has not "
                      "been installed. Please make sure they are installed.")
tf.keras.backend.set_floatx('float64')

def recall(y_true, y_pred):
    """
    Function to calculate model recall during training
    """
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall_m = true_positives / (possible_positives + K.epsilon())
    return recall_m

def precision(y_true, y_pred):
    """
    Function to calculate model precision during training
    """
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision_m = true_positives / (predicted_positives + K.epsilon())
    return precision_m

def F1(y_true, y_pred):
    """
    Function to calculate model F1 during training
    """
    precision_m = precision(y_true, y_pred)
    recall_m = recall(y_true, y_pred)
    return 2*((precision_m*recall_m)/(precision_m+recall_m+K.epsilon()))

def model_plot(fpr_keras, tpr_keras,output,history):
    """
    Function to plot model parameters during training
    """
    auc_keras = auc(fpr_keras, tpr_keras)
    out_roc= output + "ROC_curve.png"
    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_keras, tpr_keras, label="Keras (area = {:.3f})"
             "".format(auc_keras))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')
    plt.savefig(out_roc)
    plt.close('all')
    
    out_acc= output + "accuracy_plot.png"
    plt.plot(history['accuracy'])
    plt.plot(history['val_accuracy'])
    plt.title('Model accuracy')
    plt.ylabel('Accuracy')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
    plt.savefig(out_acc)
    plt.close('all')
    
    out_loss= output + "loss_plot.png"
    plt.plot(history['loss'])
    plt.plot(history['val_loss'])
    plt.title('Model loss')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
    plt.savefig(out_loss)
    plt.close('all')
    
    out_precision= output + "precision_plot.png"
    plt.plot(history['precision'])
    plt.plot(history['val_precision'])
    plt.title('Model precision')
    plt.ylabel('Precision')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
    plt.savefig(out_precision)
    plt.close('all')
    
    out_recall= output + "recall_plot.png"
    plt.plot(history['recall'])
    plt.plot(history['val_recall'])
    plt.title('Model recall')
    plt.ylabel('Recall')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
    plt.savefig(out_recall)
    plt.close('all')


def openfile(file):
    '''
    Opens a file
    '''
    if file.endswith('.gz'):
        opened_file = gzip.open(file,'rt')
    elif file.endswith('bz') or file.endswith('bz2'):
        opened_file = bz2.open(file,'rt')
    else:
        opened_file = open(file,'rt')
    return opened_file

def window_mutation(window_list, bam_file,mq,reference,depth,nine_mer):
    """
    Function to calculate mutation frequenies and qualities (features). The
    features then used by prediction or extraction mudules to be used for
    prediction or training a new model.
    """
    encoder= {'A':[1,0,0,0], 'T':[0,1,0,0], 'C':[0,0,1,0], 'G':[0,0,0,1],
              'N':[0,0,0,0],'U':[0,1,0,0]}
    bam= pysam.AlignmentFile(bam_file, 'rb')
    freq_dict= dict()
    for window in window_list:
        chrom,windowstart,windowend,base_pos = window# thsese are 9-mers
        try:
            fasta= pysam.FastaFile(reference)
        except:
            raise Exception("Cannot load reference file.")
        try:
            window_seq= fasta.fetch(reference=chrom,
                                    start=windowstart,
                                    end=windowend)
        except: # Window is not in the reference
            warnings.warn("Reference sequence for {}:{}-{} "
                          "was not found in the reference file"
                          "".format(chrom,windowstart,windowend))
            continue
        window_seq=str(window_seq.upper())
        pileupcolumns= bam.pileup(chrom, windowstart, windowend,
                                  truncate= True,min_base_quality = 0,
                                  min_mapping_quality = mq)
        refbase_index= 0
        del_list= []
        mis_list= []
        ins_list= []
        qual_list= []
        cov_list= []
        for pileupcolumn in pileupcolumns:
            pileupcolumn.set_min_base_quality(0)
            ref_base= window_seq[refbase_index]
            refbase_index += 1
            coverage= pileupcolumn.get_num_aligned()
            if coverage >= depth:
                deletion= 0
                insertion= 0
                mismatch= 0
                quality= statistics.mean(pileupcolumn.get_query_qualities())
                bases= pileupcolumn.get_query_sequences(mark_matches=True,
                                                        mark_ends=True,
                                                        add_indels=True)
                for base in bases:
                    base= base.upper()
                    if base == ref_base:
                        continue
                    elif base == '*':
                        deletion += 1
                    elif base in ['A','T','C','G','U']:
                        mismatch += 1
                    elif '+' in base:
                        insertion += 1
                    else:
                        continue
                del_list.append(deletion/coverage)
                mis_list.append(mismatch/coverage)
                ins_list.append(insertion/coverage)
                qual_list.append(quality)
                cov_list.append(coverage)
        if len(cov_list) == 9:# to check is list is empty and ignor partial in case no reads were mapped
            #converting to 5-mer
            if not nine_mer:
                for index in range(0,5):
                    kmer_dumy= []
                    for base in window_seq[index: index + 5]:
                        kmer_dumy.extend(encoder[base])
                    mean_quals= statistics.mean(qual_list[index: index + 5])
                    median_quals= statistics.median(qual_list[index: index + 5])
                    std_quals= statistics.stdev(qual_list[index: index + 5])
                    mean_mis= statistics.mean(mis_list[index: index + 5])
                    std_mis= statistics.stdev(mis_list[index: index + 5])
                    mean_dels= statistics.mean(del_list[index: index + 5])
                    std_del= statistics.stdev(del_list[index: index + 5])
                    mean_ins= statistics.mean(ins_list[index: index + 5])
                    std_ins= statistics.stdev(ins_list[index: index + 5])
                    out_frequencies= (cov_list[index: index + 5]+
                                      kmer_dumy+
                                      qual_list[index: index + 5]+
                                      mis_list[index: index + 5]+
                                      del_list[index: index + 5]+
                                      ins_list[index: index + 5]+
                                      [mean_quals, median_quals, std_quals,
                                       mean_mis, std_mis, mean_dels, std_del,
                                       mean_ins,std_ins])
                    key= (chrom, windowstart+index,windowstart+index + 5,
                          base_pos, window_seq[index: index + 5])
                    freq_dict[key]= (window , out_frequencies)
            else:
                kmer_dumy= []
                for base in window_seq:
                    kmer_dumy.extend(encoder[base])
                mean_quals= statistics.mean(qual_list)
                median_quals= statistics.median(qual_list)
                std_quals= statistics.stdev(qual_list)
                mean_mis= statistics.mean(mis_list)
                std_mis= statistics.stdev(mis_list)
                mean_dels= statistics.mean(del_list)
                std_del= statistics.stdev(del_list)
                mean_ins= statistics.mean(ins_list)
                std_ins= statistics.stdev(ins_list)
                out_frequencies= (cov_list+
                                  kmer_dumy+
                                  qual_list+
                                  mis_list+
                                  del_list+
                                  ins_list+
                                  [mean_quals, median_quals, std_quals,
                                   mean_mis, std_mis, mean_dels, std_del,
                                   mean_ins, std_ins])
                key= (chrom, windowstart, windowend,
                      base_pos, window_seq)
                freq_dict[key]= (window , out_frequencies)
        else:
            continue
    return freq_dict

def main_extraction(args):
    """
    Extraction mudule to extract features for training a new model
    """
    bam_file= os.path.abspath(args.bam)
    mq= args.mappingQuality
    threads= args.threads
    chunk= args.chunk_size
    reference= os.path.abspath(args.reference)
    mod_status= args.mod_status
    window_list=list()
    feed_list= list()
    vcf_file= os.path.abspath(args.input)
    vcf= openfile(vcf_file)
    all_lines= []
    for line in vcf:
        if line.startswith('#'):
            continue
        all_lines.append(1)
    all_lines = [all_lines[x:x+chunk]
                         for x in range(0, len(all_lines), chunk)]
    all_lines = [all_lines[x:x+threads]
                         for x in range(0, len(all_lines), threads)]
    vcf.close()
    vcf= openfile(vcf_file)
    with tqdm(total=len(all_lines),desc="Processing: ",
                      bar_format="{l_bar}{bar} [ Estimated time left:"
                      " {remaining} ]") as pbar:
        for line in vcf:
            if line.startswith('#'):
                continue
            line=line.rstrip().split('\t')
            chrom= line[0]
            windowinfo= (chrom, 
                         int(line[1])-4-1,
                         int(line[1])+5-1,
                         int(line[1])-1)
            if windowinfo not in window_list:
                window_list.append(windowinfo)
                if len(window_list) > chunk:
                    feed_list.append(window_list)
                    window_list= list()
                if len(feed_list)== threads:
                    p = mp.Pool(threads)
                    results = p.starmap(window_mutation,
                            list(zip(feed_list,
                                     repeat(bam_file),
                                     repeat(mq),
                                     repeat(reference),
                                     repeat(args.depth),
                                     repeat(args.nine_mer))))
                    p.close()
                    p.join()
                    for freq_dict in results:
                        if freq_dict is not None:
                            for key,val in freq_dict.items():
                                sys.stdout.write(','.join(map(str,key)) + ',' + 
                                      ','.join(map(str,val[1])) +
                                      ',' + str(mod_status)+'\n')
                    feed_list = []
                    pbar.update(1)
        else:
            if feed_list or window_list:
                feed_list.append(window_list)
                p = mp.Pool(len(feed_list))
                results = p.starmap(window_mutation,
                        list(zip(feed_list,
                                 repeat(bam_file),
                                 repeat(mq),
                                 repeat(reference),
                                 repeat(args.depth),
                                 repeat(args.nine_mer))))
                p.close()
                p.join()
                for freq_dict in results:
                    if freq_dict is not None:
                        for key,val in freq_dict.items():
                            sys.stdout.write(','.join(map(str,key)) + ',' + 
                                  ','.join(map(str,val[1])) + 
                                  ',' + str(mod_status)+'\n')
                pbar.update(1)                    
                feed_list = []
    vcf.close()
    sys.stderr.write("Job Finished\n")

def main_prediction(args):
    """
    Prediction module to vote for SNVs using a pretrained model
    """
    bam_file= os.path.abspath(args.bam)
    model_file= os.path.abspath(args.model_file)
    model = load_model(model_file, custom_objects={'F1': F1,
                                                   'precision':precision,
                                                   'recall':recall})
    scaler= load(open(model_file+'.pkl', 'rb'))
    mq= args.mappingQuality
    threads= args.threads
    chunk= args.chunk_size
    reference= os.path.abspath(args.reference)
    vcf_file= os.path.abspath(args.input)
    output= os.path.abspath(args.output)
    if (not os.path.isfile(output+"_Predictions.vcf") and not 
        os.path.isfile(output+"_Weighted_Qualities.vcf")):
        out_preds= open(output+"_Predictions.vcf", 'w')
        out_ready= open(output+"_Weighted_Qualities.vcf",'w')
    else:
        raise FileExistsError("Selected output files already exists.")
    vcf= openfile(vcf_file)
    all_lines= []
    for line in vcf:
        if line.startswith('#'):
            out_ready.write(line)
        all_lines.append(1)
    vcf.close()
    all_lines= [all_lines[x:x+chunk]
                         for x in range(0, len(all_lines), chunk)]
    all_lines= [all_lines[x:x+threads]
                         for x in range(0, len(all_lines), threads)]
    vcf= openfile(vcf_file)
    feed_list= list()
    windowinfo= dict()
    model_input= []
    info_input= []
    with tqdm(total=len(all_lines),desc="Processing: ",
                      bar_format="{l_bar}{bar} [ Estimated time left: {remaining} ]"
                      ) as pbar:
        for line in vcf:
            if line.startswith('#'):
                continue
            line_list= line.rstrip().split('\t')
            chrom= line_list[0]
            windowinfo[(chrom,
                        int(line_list[1])-4-1,
                        int(line_list[1])+5-1,
                        int(line_list[1])-1)]= line.rstrip()
            if len(list(windowinfo.keys())) == (chunk * threads):
                feed_list= [list(windowinfo.keys())[x:x+chunk]
                            for x in range(0, len(list(windowinfo.keys())),
                                           chunk)]
                p = mp.Pool(threads)
                results = p.starmap(window_mutation,
                        list(zip(feed_list,
                                 repeat(bam_file),
                                 repeat(mq),
                                 repeat(reference),
                                 repeat(args.depth),
                                 repeat(args.nine_mer))))
                p.close()
                p.join()
                for freq_dict in results:
                    if freq_dict is not None:
                        for key,val in freq_dict.items():
                            if not args.nine_mer:
                                coverage= int(statistics.mean(val[1][0:5]))
                                info_input.append((key,coverage,val[0]))
                                model_input.append(val[1][5:])
                            else:
                                coverage= int(statistics.mean(val[1][0:9]))
                                info_input.append((key,coverage,val[0]))
                                model_input.append(val[1][9:])
                model_input= np.asarray(model_input, dtype= np.float64)
                model_input= scaler.transform(model_input)
                predictions= model.predict(model_input)
                for i in range(len(predictions)):
                    key= info_input[i][0]
                    cov= info_input[i][1]
                    key_windoinfo= info_input[i][2]
                    pred= float(predictions[i])
                    out_preds.write(windowinfo[key_windoinfo] + '\t' +
                          '\t'.join(map(str,key)) + '\t' +
                          str(cov) + '\t' + str(pred)+'\n')
                model_input= []
                info_input= []
                windowinfo= dict()
                feed_list = []
                pbar.update(1)
        else:
            if windowinfo:
                feed_list= [list(windowinfo.keys())[x:x+chunk]
                             for x in range(0, len(list(windowinfo.keys())),
                                            chunk)]
                p = mp.Pool(threads)
                results = p.starmap(window_mutation,
                        list(zip(feed_list,
                                 repeat(bam_file),
                                 repeat(mq),
                                 repeat(reference),
                                 repeat(args.depth),
                                 repeat(args.nine_mer))))
                p.close()
                p.join()
                for freq_dict in results:
                    if freq_dict is not None:
                        for key,val in freq_dict.items():
                            if not args.nine_mer:
                                coverage= int(statistics.mean(val[1][0:5]))
                                info_input.append((key,coverage,val[0]))
                                model_input.append(val[1][5:])
                            else:
                                coverage= int(statistics.mean(val[1][0:9]))
                                info_input.append((key,coverage,val[0]))
                                model_input.append(val[1][9:])
                model_input= np.asarray(model_input, dtype= np.float64)
                model_input= scaler.transform(model_input)
                predictions= model.predict(model_input)
                for i in range(len(predictions)):
                    key= info_input[i][0]
                    cov= info_input[i][1]
                    key_windoinfo= info_input[i][2]
                    pred= float(predictions[i])
                    out_preds.write(windowinfo[key_windoinfo] + '\t' +
                          '\t'.join(map(str,key)) + '\t' +
                          str(cov)+'\t'+str(pred)+'\n')
                windowinfo= dict()
                feed_list = []
                pbar.update(1)
    vcf.close()
    out_preds.close()
    ready_dict= defaultdict(list)
    with openfile(output+"_Predictions.vcf") as pred:
        for line in pred:
            line=line.rstrip().split('\t')
            ready_dict[tuple(line[0:-7])].append(float(line[-1]))
    for key,val in ready_dict.items():
        if key[5].replace('.','',1).isdigit():
            if not args.nine_mer:
                weighted_qual= float(key[5]) * statistics.mean(val)
                out_write= list(key[0:5])+[round(weighted_qual,4)]+list(key[6:])
                out_ready.write('\t'.join(map(str,out_write))+'\n')
            else:
                weighted_qual= float(key[5]) * val[0]
                out_write= list(key[0:5])+[round(weighted_qual,4)]+list(key[6:])
                out_ready.write('\t'.join(map(str,out_write))+'\n')
    out_ready.close()
    sys.stderr.write("Job Finished\n")

def main_train(args):
    """
    Train module to train a new model
    """
    train_data= os.path.abspath(args.train)
    test_data= os.path.abspath(args.test)
    output= os.path.abspath(args.out_dir)
    nine_mer= args.nine_mer
    epoch_num= args.epochs
    batch_num= args.batch_size

    train = pd.read_csv(train_data, header= None)
    test= pd.read_csv(test_data, header= None)
    if not nine_mer:
        X_train= train.iloc[:, 10:-1].values
        y_train = train.iloc[:, -1].values
        test_val= test.iloc[:, 10:-1].values
        test_label= test.iloc[:, -1].values
    else:
        X_train= train.iloc[:, 14:-1].values
        y_train = train.iloc[:, -1].values
        test_val= test.iloc[:, 14:-1].values
        test_label= test.iloc[:, -1].values
    inputdim= X_train.shape[1]
    # Feature Scaling
    sc = StandardScaler()
    X_train = sc.fit_transform(X_train)
    test_val= sc.transform(test_val)
    outmodel= output + "_trained_model.h5"
    dump(sc, open(outmodel+'.pkl', 'wb'))
    # Initialising the ANN
    initializer = keras.initializers.glorot_uniform(seed=1)
    classifier = Sequential()
    classifier.add(Dense(activation="relu",
                         input_dim=inputdim,
                         units=inputdim,
                         kernel_initializer=initializer))
    classifier.add(Dense(activation="relu",
                         units=inputdim * 2,
                         kernel_initializer=initializer))
    classifier.add(Dense(activation="relu",
                         units=inputdim,
                         kernel_initializer=initializer))
    classifier.add(Dense(activation="relu",
                         units=round(inputdim/2),
                         kernel_initializer=initializer))
    classifier.add(Dense(activation="sigmoid",
                         units=1,
                         kernel_initializer=initializer))
    classifier.compile(optimizer = 'adam',
                       loss = 'binary_crossentropy',
                       metrics = ['accuracy',F1,precision, recall])
    # Fitting the ANN to the Training set
    checkpointer = ModelCheckpoint(filepath=outmodel,
                                   monitor='val_loss',
                                   save_best_only=True,
                                   verbose=1)

    history= classifier.fit(X_train,
                            y_train,
                            validation_data= (test_val, test_label),
                            batch_size = batch_num,
                            epochs = epoch_num,
                            callbacks=[checkpointer],
                            shuffle=False,verbose=1)

    y_predictions = classifier.predict(test_val)
    y_pred = (y_predictions > 0.5)
    out_cm= output + "_Test_prediction_ConfusionMatrix.tsv"
    cm = confusion_matrix(test_label, y_pred)
    print(cm,file=open(out_cm,'w'))

    y_pred_keras = classifier.predict(test_val).ravel()
    fpr_keras, tpr_keras, thresholds_keras = roc_curve(test_label,
                                                       y_pred_keras)
    if args.plot:
        model_plot(fpr_keras,
                   tpr_keras,
                   output,
                   history.history)


def prediction_parser(subparsers):
    """
    Specific argument parser for prediction module.
    """
    sub_prediction = subparsers.add_parser(
        "prediction",
        help="predict using a model",
        description=("Predict based on a model."))
    sp_input = sub_prediction.add_argument_group("required arguments")
    sp_input.add_argument("--input", "-i", action="store", type=str,
                          required=True,help="The path to the input "
                          "vcf or bed file. NOTE. Files must end with "
                          ".bed or .vcf. vcf files are 1-based and beds "
                          "are zero-based")
    sp_input.add_argument("--bam", "-b", action="store", type=str,
                          required=True,help="The path to the alignment bam "
                          "file")
    sp_input.add_argument("--reference", "-r", action="store", type=str,
                          required=True, default=None, help="The path to the "
                          "reference file. File must be indexed by samtools "
                          "faidx.")
    sp_input.add_argument("--output", "-o", action="store", type=str,
                          required=True, default=None, help="The path to the "
                          "output directory and prefix for output file.")
    sp_input = sub_prediction.add_argument_group("optional arguments")
    sp_input.add_argument("--model_file", "-mf", type=str, action="store",
                          required=False,
                          default=os.path.join(os.path.dirname(
                                                   os.path.realpath(__file__)),
                                               "model",
                                               "NA12878_20FC_model.h5"),
                          help="Path to the trained model. Default is "
                               "NA12878_20FC_model.h5")
    sp_input.add_argument("--mappingQuality", "-mq", action="store", type=int,
                          default= 0,required=False,help="Cutt off for filtering "
                          "out low quality mapped reads from bam. Default is 0")
    sp_input.add_argument("--depth", "-d", action="store", type=int,
                          default= 1, required=False, help="Cutt off for "
                          "filtering out regions with low depth to have "
                          "frequencies. Default >= 1")
    sp_input.add_argument("--window_bam", "-w", action="store", type=str,
                          required=False, help="if you want to only do for a "
                          "region or chromosom You must insert region like "
                          "this chr1 or chr1:1000-100000.")
    sp_input.add_argument("--nine_mer", "-nm", action="store_true",
                          required=False, help="Prediction for 9-mer. "
                          "Default is five mer")
    sp_input.add_argument("--threads", "-t", action="store", type=int,
                          required=False, default=4, help="Number of threads. "
                          "Default is 4.")
    sp_input.add_argument("--chunk_size", "-cs", action="store", type=int,
                          required=False, default=100, help="Chunk size. "
                          "Default is 100.")
    sub_prediction.set_defaults(func=main_prediction)

def extraction_parser(subparsers):
    """
    Specific argument parser for extraction module.
    """
    sub_extraction = subparsers.add_parser(
        "extraction",
        help="extract features",
        description=("Extract mutation frequencicies in 5-mer window."))
    se_input = sub_extraction.add_argument_group("required arguments")
    se_input.add_argument("--input", "-i", action="store", type=str,
                          required=True, help="The path to the input vcf or "
                          "bed file. NOTE. Files must end with .bed or .vcf. "
                          "vcf files are 1-based and beds are zero-based")
    se_input.add_argument("--mod_status", "-ms", type=int,action="store",
                          required=True, help= "0 or 1. If you are extracting "
                          "frequencies to train a model, give the"
                        " modification status for your bed file either it is "
                        "modified (1) or unmodified (0) regions.")
    se_input.add_argument("--bam", "-b", action="store", type=str,
                          required= True, help= "The path to the alignment "
                          "bam file")
    se_input.add_argument("--reference", "-r", action="store", type=str,
                          required=True, default=None, help="The path to the "
                          "reference file. File must be indexed by samtools "
                          "faidx")
    se_input = sub_extraction.add_argument_group("optional arguments")
    se_input.add_argument("--mappingQuality", "-mq", action="store", type=int,
                          default= 0,required=False,help="Cutt off for filtering "
                          "out low quality mapped reads from bam. Default is 0")
    se_input.add_argument("--depth", "-d", action="store", type=int,default= 1,
                          required=False,help="Cutt off for filtering out regions "
                          "with low depth to have frequencies. Default >=1")
    se_input.add_argument("--window_bam", "-w", action="store", type=str,
                          required=False, help= "if you want to only do for a "
                          "region or chromosom, you must insert region like "
                          "this chr1 or chr1:1000-100000.")
    se_input.add_argument("--nine_mer", "-nm", action="store_true",
                          required=False, help="Extraction for 9-mer. "
                          "Default is five mer")
    se_input.add_argument("--threads", "-t", action="store", type=int,
                          required=False, default= 4, help="Number of threads")
    se_input.add_argument("--chunk_size", "-cs", action="store", type=int,
                          required= False, default= 400, help= "Number of "
                          "sites send to each processes for parrallel "
                          "processing. Default is 50.")
    sub_extraction.set_defaults(func=main_extraction)

def train_parser(subparsers):
    """
    Specific argument parser for train module.
    """
    sub_train = subparsers.add_parser(
        "train",
        help="train a new model",
        description=("train a new model"))
    st_input = sub_train.add_argument_group("required arguments")
    st_input.add_argument("--train", "-tr", action="store", type=str,
                          required= True, help= "The path to the shuffled and "
                          "ready file for training")
    st_input.add_argument("--test", "-te", action="store", type=str,
                          required= True, help= "The path to the shuffled and "
                          "ready file for testing.")
    st_input.add_argument("--out_dir", "-o", action="store", type=str,
                          required= True, help= "Output directory and prefix "
                          "for saving the model and figures")
    st_input = sub_train.add_argument_group("optional arguments")
    st_input.add_argument("--epochs", "-e", action="store", type=int,
                          default= 100, required= False, help= "Number of "
                          "epochs. Default is 100")
    st_input.add_argument("--batch_size", "-batch", action="store", type=int,
                          required= False, default= 400, help= "batch size for "
                          "model training. Default is 400.")
    st_input.add_argument("--plot", "-plt", action= "store_true",
                          required= False, help= "Select this option if you "
                          "wish to output training plots.")
    st_input.add_argument("--nine_mer", "-nm", action= "store_true",
                          required= False, help="Training for 9-mer. "
                          "Default is five mer")
    sub_train.set_defaults(func=main_train)

def main():
    """
    Docstring placeholder.
    """
    parser = argparse.ArgumentParser(
#        formatter_class=argparse.RawTextHelpFormatter, #  It's bad habit
        prog="snvoter",
        description=("SNVoter: Voting for SNVs"))
    subparsers = parser.add_subparsers(title="Modules")
    extraction_parser(subparsers)
    prediction_parser(subparsers)
    train_parser(subparsers)
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
