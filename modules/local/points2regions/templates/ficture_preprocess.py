#!/usr/bin/env python

import sys
import os
import re
import logging
import pandas as pd
import gzip

def format_xenium():

    print("[START]")

    negctrl_regex = "BLANK\\|NegCon"
    if ( "${params.negative_control_regex}" ):
        negctrl_regex = "${params.negative_control_regex}"

    transripts = "${transcripts}"
    features = "${features}"

    unit_info=['X','Y','gene','cell_id','overlaps_nucleus']
    oheader = unit_info + ['Count']

    feature=pd.DataFrame()
    xmin=sys.maxsize
    xmax=0
    ymin=sys.maxsize
    ymax=0

    output = "processed_transcripts.tsv.gz"
    feature_file = "feature.clean.tsv.gz"
    min_phred_score = 15

    with gzip.open(output, 'wt') as wf:
        wf.write('\\t'.join(oheader) + '\\n')

    for chunk in pd.read_csv(transripts,header=0,chunksize=500000):
        chunk = chunk.loc[(chunk.qv > min_phred_score)]
        chunk.rename(columns = {'feature_name':'gene'}, inplace=True)
        if negctrl_regex != '':
            chunk = chunk[~chunk.gene.str.contains(negctrl_regex, flags=re.IGNORECASE, regex=True)]
        chunk.rename(columns = {'x_location':'X', 'y_location':'Y'}, inplace=True)
        chunk['Count'] = 1
        chunk[oheader].to_csv(output,sep='\\t',mode='a',index=False,header=False,float_format="%.2f")
        logging.info(f"{chunk.shape[0]}")
        feature = pd.concat([feature, chunk.groupby(by='gene').agg({'Count':"sum"}).reset_index()])
        x0 = chunk.X.min()
        x1 = chunk.X.max()
        y0 = chunk.Y.min()
        y1 = chunk.Y.max()
        xmin = min(int(xmin), int(x0))
        xmax = max(int(xmax), int(x1))
        ymin = min(int(ymin), int(y0))
        ymax = max(int(ymax), int(y1))

    if os.path.exists(features):
        feature_list = []
        with open(features, 'r') as ff:
            for line in ff:
                feature_list.append(line.strip('\\n'))
        feature = feature.groupby(by='gene').agg({'Count':"sum"}).reset_index()
        feature = feature[[x in feature_list for x in feature['gene']]]
        feature.to_csv(feature_file,sep='\\t',index=False)

    f = os.path.join( os.path.dirname(output), "coordinate_minmax.tsv" )
    with open(f, 'w') as wf:
        wf.write(f"xmin\\t{xmin}\\n")
        wf.write(f"xmax\\t{xmax}\\n")
        wf.write(f"ymin\\t{ymin}\\n")
        wf.write(f"ymax\\t{ymax}\\n")

    # Output version information
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'ficture_preprocess: v.1.0.0"\\n')

    print("[FINISH]")

if __name__ == '__main__':
    format_xenium()
