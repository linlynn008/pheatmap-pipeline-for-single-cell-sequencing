import pandas as pd
import numpy as np
from functools import reduce

#https://thispointer.com/how-to-change-current-working-directory-in-python/

import os

dirpath = os.getcwd()
   print("current directory is : " + dirpath)
#current directory is : C:\Users\lguo\AppData\Local\Programs\Python\Python37-32
   foldername = os.path.basename(dirpath)
   print("Directory name is : " + foldername)

for file in os.listdir("C:/Users/lguo/Documents/PACT_Pharma/Smartseq/EXP19_1502/fpkms"):
    if file.endswith("bam.tab"):
        print(os.path.join("C:/Users/lguo/Documents/PACT_Pharma/Smartseq/EXP19_1502/fpkms", file))

fn = ['EXP1502A5-right_S12_L001.bam.tab', 'EXP1502A7-right_S20_L001.bam.tab', 'EXP1502A9-right_S26_L001.bam.tab',
      'EXP1502B1_S1_L001.bam.tab', 'EXP1502B7-right_S21_L001.bam.tab', 'EXP1502C1_S2_L001.bam.tab',
      'EXP1502C2_S8_L001.bam.tab', 'EXP1502C6-right_S17_L001.bam.tab', 'EXP1502D1_S3_L001.bam.tab',
      'EXP1502D5-right_S13_L001.bam.tab', 'EXP1502E3_S7_L001.bam.tab', 'EXP1502E5-right_S14_L001.bam.tab',
      'EXP1502E7-right_S22_L001.bam.tab', 'EXP1502F2_S5_L001.bam.tab', 'EXP1502F4-right_S10_L001.bam.tab',
      'EXP1502F5-right_S15_L001.bam.tab', 'EXP1502F6-right_S18_L001.bam.tab', 'EXP1502F8-right_S24_L001.bam.tab',
      'EXP1502F9-right_S27_L001.bam.tab']

def filtered_data (f1):
   f2 = pd.read_csv (f1, sep='\t', index_col=0)
   f3 = f2[f2.apply(lambda z: (z['Reference'].find('_') < 0) and (z.TPM > 0), axis=1)].drop_duplicates(
               subset='Gene Name', keep='first')[['TPM']]
   f3.columns = [ f1[0:f1.find('-')] ];
   return f3

fd = list(map(filtered_data, fn))

mm = reduce(lambda x, y: x.merge(y, how='outer', left_index=True, right_index=True), fd)
mm = mm.fillna(0)
mm.to_csv('unified.csv')







