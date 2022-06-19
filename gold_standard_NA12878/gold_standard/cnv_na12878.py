#Import required modules
import os, argparse, pybedtools, csv
import pandas as pd
import numpy as np
from natsort import index_natsorted
from collections import Counter
from sys import stderr

#----------------------------------------------------------------------------------------------------------------
### SUBFUNCTIONS REQUIRED FOR KEY FUNCTIONS ###
#----------------------------------------------------------------------------------------------------------------
def sort_chrom(Chr):
    """ Sort BED file in alpha-numeric order by chromosome.

    Arguments
    ---------
        Chr: Chromosome (i.e. chr1)

    Returns
    -------
        New: New variable to order autosomal chromosomes as well as sex and mitochondrial chromosomes.
    """

    if Chr:
        New = Chr[3:]
        if New == 'X':
            New = 23
        elif New == 'Y':
            New = 24
        elif New == 'M':
            New = 25
        elif New == 'N':
            New = 26
        else:
            New = int(New)
    else:
        New = 0
    return New

#----------------------------------------------------------------------------------------------------------------
### KEY FUNCTIONS ###
#----------------------------------------------------------------------------------------------------------------
def new_version_bedfile(file, study, out):
    '''
    Generate a new version of the initial BED file to start the analysis.
    A tag with the dbVar study ID is added and the size of each CNV call is calculated.

    Arguments
    ---------
        file: Input BED file
        study: dbVar study ID (the tag)
        out: Initial output path

    Returns
    -------
        The new BED file into the specified output path
    '''

    print("Processing %s..." %(file))

    # Load the input file with the CNVs characterized
    df = pd.read_csv(file, sep='\t', comment='#', header=None)
    header = ['#chrom', 'start', 'end', 'cnv_type', 'copy_num', 'variant_call_id', 'variant_region_id', 
              'variant_region_access', 'validation', 'placements_per_assembly', 'num_probe_CNV']
    df.columns = header[:len(df.columns)]
    df = df.drop("variant_region_access", axis=1)
        
    #Obtain the size in bp and kb of the CNV interval
    df['length_bp'] = df["end"].apply(lambda x: int(x)) - df["start"].apply(lambda x: int(x))
    df['length_kb'] = df['length_bp'].apply(lambda x: round((int(x)/1000),1))
    
    # The numenclature of the CNV type is standardized for further analysis
    df = df.replace({'cnv_type' : {"copy_number_gain": "duplication", "copy_number_loss": "deletion"}})

    #Add the tag belonging to the study
    df.insert(loc=9, column='dbVar_study_id', value=study)
    
    # Sort by chromosome and then by start and end position
    df.sort_values(by=["#chrom", "start", "end"], key=lambda x: np.argsort(index_natsorted(df["#chrom"])),
                   inplace=True, ignore_index=True)

    # Write the BED file (name: 'inputfile_v2.bed')
    fields = ['#chrom', 'start', 'end', 'cnv_type', 'copy_num', 'variant_call_id', 'variant_region_id','validation',
              'placements_per_assembly', 'dbVar_study_id', 'num_probe_CNV', 'length_bp', 'length_kb']
    df.to_csv(os.path.join(out, "%s_v2.bed" % (file)), sep='\t', columns=fields, index=False)

    print("'%s_v2.bed' created in %s... \n" %(file, out))

def consensus_cnv(study1, id1, arrayid1, study2, id2, arrayid2, study3, id3, arrayid3, dst):
    ''' 
    Determine the consensus target regions of CNVs that will set the ground truth, from the CNV calls 
    identified by initial studies; using pybedtools (a Python extension of Aaron Quinlan’s BEDtools suite).
    Assign a confidence score and a classification level to each CNV interval.

    Arguments
    ---------
        study1: new BED file generated previously, corresponding to one of the initial studies.
        id1: dbVar study ID (or a label), referring to 'study1' 
        arrayid1: List of one or more files with the experimental designs (genomic coordinates of the array probes 
                used), corresponding to the 'study1'

        study2: new BED file generated previously, corresponding to another of the initial studies
        id2: dbVar study ID (or a label), referring to 'study2'
        arrayid2: List of one or more files with the experimental designs (genomic coordinates of the array probes 
                used), corresponding to the 'study2'

        study3: new BED file generated previously, corresponding to another of the initial studies
        id3: dbVar study ID (or a label), referring to 'study3'
        arrayid3: List of one or more files with the experimental designs (genomic coordinates of the array probes 
                used), corresponding to the 'study3'
                
        dst: Initial output path where save files generated

    Returns
    -------
        - The 'ID1_ID2_ID3_merge.bed' file with the overlapping features.
        - The 'consensus_target_regions.bed' file with the consensus CNV intervals determined from the initial 
        studies; and the corresponding confidence score and classification level assigned to each region.
        
        Both files are saved in the specified output path.
    '''

    print("Determining the consensus CNVs intervals that set the ground truth...")

    ## 1. Generate the merge file with overlapping features
    # Create a BedTool object from the study1 
    st1 = pybedtools.BedTool((os.path.join(dst, "%s_v2.bed" % (study1))))
    # Create a BedTool object from the study2 
    st2 = pybedtools.BedTool((os.path.join(dst, "%s_v2.bed" % (study2))))
    # Create a BedTool object from the study3 
    st3 = pybedtools.BedTool((os.path.join(dst, "%s_v2.bed" % (study3))))

    # Use BedTools B.coverage(A) with A and B files. After each interval in A, reports:
    # 1. The number of features in B that overlapped (by at least one base pair) the A interval.
    # 2. The number of bases in A that had non-zero coverage from features in B.
    # 3. The length of the entry in A.
    # 4. The fraction of bases in A that had non-zero coverage from features in B.
    st1_cov = st1.coverage([st2, st3])
    st2_cov = st2.coverage([st1, st3])
    st3_cov = st3.coverage([st1, st2])
 
    # Concatenate previous BedTools objects, and sort by chromosome and then by start position 
    cat_sort = (st1_cov.cat(st2_cov, postmerge=False)).cat(st3_cov, postmerge=False).sort()
            
    # BedTools merge() combines overlapping features in an interval file into a single feature which spans all of
    # the combined features
    mergeCNV = cat_sort.merge(c=[2,3,4,5,6,7,8,9,10,11,14,15,16,17], o="collapse",
                              delim="|").saveas(os.path.join(dst, "%s_%s_%s_merge.bed" % (id1, id2, id3)))
                                        # Write the merge BED file (name: 'id1_id2_id3_merge.bed')
        
        
    ##2. Calculate the number of probes matching the consensus CNV interval, taking into account the experimental
    ##design of each initial study (this will be necessary to calculate the confidence score later on) using 
    # BedTools A.intersect(B)
    
    #merge_probestudy1 = mergeCNV.intersect(arrayid1, c=True) #No me deja ejecutar este comando porque  
        #el archivo de 42M sondas ocupa mucho. Por tanto, lo divido en cuatro archivos de ~10M sondas.
    merge_probestudy1A = mergeCNV.intersect(arrayid1[0], c=True)
    merge_probestudy1B = mergeCNV.intersect(arrayid1[1], c=True)
    merge_probestudy1C = mergeCNV.intersect(arrayid1[2], c=True)
    merge_probestudy1D = mergeCNV.intersect(arrayid1[3], c=True)
    cat_probesid1 = merge_probestudy1A.cat(merge_probestudy1B, postmerge=False).cat(
        merge_probestudy1C, postmerge=False).cat(merge_probestudy1D, postmerge=False).sort()
    merge_probestudy1 = cat_probesid1.groupby(g=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17], c=18, o=["sum"])

    merge_probestudy2 = mergeCNV.intersect(arrayid2, c=True)
    merge_probestudy3 = mergeCNV.intersect(arrayid3, c=True)

    # Keep only the field containing the number of probes
    num_probe_study1_merge = [x[17] for x in merge_probestudy1]
    num_probe_study2_merge = [x[17] for x in merge_probestudy2]
    num_probe_study3_merge = [x[17] for x in merge_probestudy3]

    # Create a dataframe that stores all information obtained for each CNV interval.
    data = pybedtools.BedTool.to_dataframe(mergeCNV, names=['#chrom', 'start', 'end', 'startx', 'endx', 'cnv_type', 
                                                            'copy_num', 'variant_call_id', 'variant_region_id', 
                                                            'validation', 'num_map', 'dbVar_study_id', 
                                                            'num_probe_CNV', 'num_features_overlap', 'bp_overlap',
                                                            'feature_size_bp', 'fraction_overlap'])      
        
    data.insert(loc = 13, column = 'probes_id1_mergeCNV', value = num_probe_study1_merge)
    data.insert(loc = 14, column = 'probes_id2_mergeCNV', value = num_probe_study2_merge)
    data.insert(loc = 15, column = 'probes_id3_mergeCNV', value = num_probe_study3_merge)
    
    #Obtain the size in bp and kb of the CNV interval
    data['interval_length_bp'] = data["end"].apply(lambda x: int(x)) - data["start"].apply(lambda x: int(x))
    data['interval_length_kb'] = data['interval_length_bp'].apply(lambda x: round((int(x)/1000),1))
    
    
    ## 3. Calculate the confidence score for each CNV interval.        
    info = []
    for index, row in data.iterrows():
        chrom = row['#chrom']
        start = int(row['start'])
        end = int(row['end'])
        startx = row['startx' ].split("|")
        endx = row['endx'].split("|")
        cnv_type = row['cnv_type'].split("|")
        copy_num = row['copy_num']
        variant_call_id = row['variant_call_id']
        variant_region_id = row['variant_region_id']
        validation = row['validation']
        num_map = row['num_map']
        dbVar_id = row['dbVar_study_id'].split("|")
        num_probe_CNV = row['num_probe_CNV'].split("|")
        probes_id1_mergeCNV = int(row['probes_id1_mergeCNV'])
        probes_id2_mergeCNV = int(row['probes_id2_mergeCNV'])
        probes_id3_mergeCNV = int(row['probes_id3_mergeCNV'])
        bp_overlap = row['bp_overlap'].split("|")
        feature_size_bp = row['feature_size_bp']
        fraction_overlap = row['fraction_overlap'].split("|")
        interval_length_bp = int(row['interval_length_bp'])
        interval_length_kb = float(row['interval_length_kb'])
        
        global confidence_score
        
        ## 3.1. CNVs with an overlap = 0%, i.e. CNVs characterized by a single study.
        if (all(int(x) == 0 for x in bp_overlap)) == True:
            # Assign a CNV classification of 'Low Confidence', because there is no consensus among different studies 
            #(only identified by one study).
            classification = 'LowConfidence'
            # If probes are available to detect CNV interval:
            if (all(int(x) != 0 for x in num_probe_CNV)) == True:
                # First, identify the study id corresponding to the input CNV.
                # Determine if probes from the array belonging to the other studies match in that region (larger CNV 
                #interval). If that CNV could have been detected by other studies (presents probes in that region), 
                #but has not, it is penalized with a 0.0 (due to lack of congruence between studies and experimental 
                #design). Otherwise, a score of 0.5 is assigned. 
                # The minimum number of probes considered to be able to characterize a CNV is 6 (which is the highest 
                #value of the minimum number of probes, in the comparison of the three studies).
                if dbVar_id[0] == id1: 
                    if probes_id2_mergeCNV > 6 or probes_id3_mergeCNV > 6:
                        confidence_score = 0.0
                    else:
                        confidence_score = 0.5
                    info.append((chrom, start, end, '|'.join(cnv_type), copy_num, variant_call_id, variant_region_id, 
                                 validation,num_map,'|'.join(dbVar_id), '|'.join(num_probe_CNV), '|'.join(bp_overlap), 
                                 '|'.join(fraction_overlap), interval_length_bp, interval_length_kb, confidence_score, 
                                 classification))
                elif dbVar_id[0] == id2:  
                    if probes_id1_mergeCNV > 6 or probes_id3_mergeCNV > 6:
                        confidence_score = 0.0
                    else:
                        confidence_score = 0.5
                    info.append((chrom, start, end, '|'.join(cnv_type), copy_num, variant_call_id, variant_region_id, 
                                 validation,num_map,'|'.join(dbVar_id), '|'.join(num_probe_CNV), '|'.join(bp_overlap), 
                                 '|'.join(fraction_overlap), interval_length_bp, interval_length_kb, confidence_score, 
                                 classification))
                elif dbVar_id[0] == id3:  
                    if probes_id1_mergeCNV > 6 or probes_id2_mergeCNV > 6:
                        confidence_score = 0.0
                    else:
                        confidence_score = 0.5    
                    info.append((chrom, start, end, '|'.join(cnv_type), copy_num, variant_call_id, variant_region_id, 
                                 validation,num_map,'|'.join(dbVar_id), '|'.join(num_probe_CNV), '|'.join(bp_overlap), 
                                 '|'.join(fraction_overlap), interval_length_bp, interval_length_kb, confidence_score, 
                                 classification))
            else: # 0 probes
                confidence_score = 0.0
                info.append((chrom, start, end, '|'.join(cnv_type), copy_num, variant_call_id, variant_region_id, 
                                 validation,num_map,'|'.join(dbVar_id), '|'.join(num_probe_CNV), '|'.join(bp_overlap), 
                                 '|'.join(fraction_overlap), interval_length_bp, interval_length_kb, confidence_score, 
                                 classification))
                
        ## 3.2. CNVs with reciprocal overlap greater than or equal to 75%:
        #This is similar to using -f 0.75 -r from BEDtools intersect. In other words, CNVs that overlap at least 75% 
        #with each other receive a high confidence score; since it is assumed that the breakpoints are similar.
        elif all(float(x)>=0.75 for x in fraction_overlap) == True:
            # A) If the overlapping CNVs are characterized as the same type of CNV by the corresponding studies, 
            #keep the information of the larger CNV interval, assign a confidence score of 1.0, and a classification
            #level of HighConfidence.
            classification = 'HighConfidence'
            if len(set(cnv_type)) == 1:
                confidence_score = 1.0
                info.append((chrom, start, end, '|'.join(set(cnv_type)), copy_num, variant_call_id, 
                             variant_region_id, validation,num_map,'|'.join(sorted(set(dbVar_id))), 
                             '|'.join(num_probe_CNV), '|'.join(bp_overlap), '|'.join(fraction_overlap), 
                             interval_length_bp, interval_length_kb, confidence_score, classification))

            # B) If there are multiple overlapping CNVs calls, characterized as different types of CNV: 
            #determine which CNV type is the most supported and maintain that information with the largest CNV 
            #interval, assigning a high classification level (since at least two studies support the CNV); and 
            #assigning a confidence score of 0.5 (penalize it because there is uncertainty between experimental 
            #designs, although the breakpoints are similar). The least supported CNV type is removed.        
            else:
                confidence_score = 0.5
                fields = startx, endx, cnv_type, copy_num.split("|"), variant_call_id.split("|"),\
                variant_region_id.split("|"), validation.split("|"), num_map.split("|"),dbVar_id, num_probe_CNV,\
                bp_overlap, fraction_overlap

                main_type = Counter(cnv_type)
                if main_type['deletion'] > main_type['duplication']:
                    index_dup = [n for n, x in enumerate(cnv_type) if x == "duplication"]
                    list(map(lambda x: [x.pop(i) for i in index_dup[::-1]], fields))

                    info.append((chrom, int(min(fields[0])), int(max(fields[1])), fields[2][0], 
                                 '|'.join(fields[3]),'|'.join(fields[4]), '|'.join(fields[5]), 
                                 '|'.join(fields[6]), '|'.join(fields[7]),'|'.join(sorted(set(fields[8]))), 
                                 '|'.join(fields[9]),'|'.join(fields[10]), '|'.join(fields[11]), 
                                 (int(max(fields[1])) - int(min(fields[0]))), 
                                 round(((int(max(fields[1]))- int(min(fields[0]))) / 1000), 1),
                                 confidence_score, classification))

                elif main_type['deletion'] < main_type['duplication']:
                    index_del = [n for n, x in enumerate(cnv_type) if x == "deletion"]
                    list(map(lambda x: [x.pop(i) for i in index_del[::-1]], fields))

                    info.append((chrom, int(min(fields[0])), int(max(fields[1])), fields[2][0], 
                                 '|'.join(fields[3]),'|'.join(fields[4]), '|'.join(fields[5]), 
                                 '|'.join(fields[6]), '|'.join(fields[7]),'|'.join(sorted(set(fields[8]))), 
                                 '|'.join(fields[9]), '|'.join(fields[10]), '|'.join(fields[11]), 
                                 (int(max(fields[1])) - int(min(fields[0]))), 
                                 round(((int(max(fields[1]))- int(min(fields[0]))) / 1000), 1),
                                 confidence_score, classification))

                else: #If main_type['deletion'] = main_type['duplication']
                    confidence_score = 0.0
                    classification = 'Indeterminate'
                    info.append((chrom, start, end, '|'.join(set(cnv_type)), copy_num, variant_call_id, 
                                 variant_region_id, validation,num_map,'|'.join(sorted(set(dbVar_id))), 
                                 '|'.join(num_probe_CNV), '|'.join(bp_overlap), '|'.join(fraction_overlap), 
                                 interval_length_bp, interval_length_kb, confidence_score, classification))
                
                        
        ## 3.3. CNVs with reciprocal overlap less than 75% (except 0%):
        #This intervals receive a low confidence score, since it is assumed that the breakpoints are different
        elif (any(float(x) < 0.75 and float(x) != 0 for x in fraction_overlap)) == True:
            # A) If the overlapping CNVs are characterized as the same type of CNV by the corresponding studies, 
            #keep the information of the larger CNV interval, assign a confidence score of 0.5 or 0.0 (depending
            #on the power of detection (number of probes)), and a classification level of HighConfidence.
            classification = 'HighConfidence'
            if len(set(cnv_type)) == 1:
                # Create a dictionary that associates to each study (key) the set of matching probe values
                # {'dbvar_id': (number of CNV calls, probes that match the larger CNV interval, probes from
                # the array of the first study that match the larger CNV interval, probes from the array of the
                # second study, probes from the third study)}
                # (example: {'estd20': [(1, 824, 824, 106, 1)], 'nstd46': [(4, 100, 824, 106, 1)]})
                c_s=[]
                dict_study = {}  
                counts = Counter(dbVar_id)
                for i in range(len(dbVar_id)):
                    studyid = dbVar_id[i]
                    probes_CNV = int(num_probe_CNV[i])
                    dict_study.setdefault(studyid, []).append((counts[studyid], probes_CNV, probes_id1_mergeCNV,
                                                               probes_id2_mergeCNV, probes_id3_mergeCNV))
                new_dict = {}
                for key, value in dict_study.items():
                    if len(value) == 1:  
                        new_dict[key] = value
                    else: 
                        probes = []
                        for i in range(len(value)):
                            probes.append(value[i][1])
                        tot_probes_CNV = sum(probes)
                        new_dict.setdefault(key, []).append((value[0][0],tot_probes_CNV,value[0][2],
                                                             value[0][3],value[0][4]))

                for key, value in new_dict.items():
                    if value[0][1] != 0: # If probes are available to detect CNV interval:
                    # First, identify the study id corresponding to the input CNV.
                    #Compare if the number of probes that have detected that CNV call is equal to the number of 
                    #probes in the corresponding array that match in that region (larger CNV interval) (with an 
                    #error of 6 probes (the minimum number of probes considered to be able to identified a CNV)). 
                    #If the number of probes is different (with an error of 6 probes) (i.e. the array could 
                    #potentially have detected the entire CNV interval), penalize the interval with a 0.0 score (due
                    #to lack of congruence between estudies). Otherwise, assign a score of 0.5 
                        if key == id1: 
                            if (value[0][2]-6) <= value[0][1] <= (value[0][2]+6):
                                c_s.append(0.5)
                            else:
                                c_s.append(0.0)
                        elif key == id2: 
                            if (value[0][3]-6) <= value[0][1] <= (value[0][3]+6):
                                c_s.append(0.5)
                            else:
                                c_s.append(0.0)
                        elif key == id3:  
                            if (value[0][4]-6) <= value[0][1] <= (value[0][4]+6):
                                c_s.append(0.5)
                            else:
                                c_s.append(0.0)
                    else: #CNV characterized without probes: confidence score is 0.0
                        c_s.append(0.0)
                #If interval has the same confidence score for each CNV call that matches in that region:
                if len(set(c_s)) == 1: 
                    confidence_score = c_s[0]  # keep one value
                else:  # else, penalize with the minimum score
                    confidence_score = min(c_s)

                info.append((chrom, start, end, '|'.join(set(cnv_type)), copy_num, variant_call_id, 
                             variant_region_id, validation,num_map,'|'.join(sorted(set(dbVar_id))), 
                             '|'.join(num_probe_CNV),'|'.join(bp_overlap), '|'.join(fraction_overlap), 
                             interval_length_bp, interval_length_kb, confidence_score, classification))
                
            # B) If there are multiple overlapping CNVs calls, characterized as different types of CNV: 
            #determine which CNV type is the most supported and maintain that information with the largest CNV 
            #interval, assigning a high classification level (since at least two studies support the CNV); and 
            #assigning a confidence score of 0.0 (penalize it due to the uncertainty between experimental 
            #designs and breakpoints). The least supported CNV type is removed.
            else:
                confidence_score = 0.0
                fields = startx, endx, cnv_type, copy_num.split("|"), variant_call_id.split("|"),\
                variant_region_id.split("|"), validation.split("|"), num_map.split("|"),dbVar_id, num_probe_CNV,\
                bp_overlap, fraction_overlap

                main_type = Counter(cnv_type)
                if main_type['deletion'] > main_type['duplication']:
                    index_dup = [n for n, x in enumerate(cnv_type) if x == "duplication"]
                    list(map(lambda x: [x.pop(i) for i in index_dup[::-1]], fields))

                    info.append((chrom, int(min(fields[0])), int(max(fields[1])), fields[2][0], 
                                 '|'.join(fields[3]),'|'.join(fields[4]), '|'.join(fields[5]), 
                                 '|'.join(fields[6]), '|'.join(fields[7]),'|'.join(sorted(set(fields[8]))), 
                                 '|'.join(fields[9]), '|'.join(fields[10]), '|'.join(fields[11]), 
                                 (int(max(fields[1])) - int(min(fields[0]))), 
                                 round(((int(max(fields[1]))- int(min(fields[0]))) / 1000), 1),
                                 confidence_score, classification))

                elif main_type['deletion'] < main_type['duplication']:
                    index_del = [n for n, x in enumerate(cnv_type) if x == "deletion"]
                    list(map(lambda x: [x.pop(i) for i in index_del[::-1]], fields))

                    info.append((chrom, int(min(fields[0])), int(max(fields[1])), fields[2][0], 
                                 '|'.join(fields[3]),'|'.join(fields[4]), '|'.join(fields[5]), 
                                 '|'.join(fields[6]), '|'.join(fields[7]),'|'.join(sorted(set(fields[8]))), 
                                 '|'.join(fields[9]), '|'.join(fields[10]), '|'.join(fields[11]), 
                                 (int(max(fields[1])) - int(min(fields[0]))), 
                                 round(((int(max(fields[1]))- int(min(fields[0]))) / 1000), 1),
                                 confidence_score, classification))

                else: #If main_type['deletion'] = main_type['duplication']
                    classification = 'Indeterminate'
                    info.append((chrom, start, end, '|'.join(set(cnv_type)), copy_num, variant_call_id, 
                                 variant_region_id, validation,num_map,'|'.join(sorted(set(dbVar_id))), 
                                 '|'.join(num_probe_CNV), '|'.join(bp_overlap), '|'.join(fraction_overlap), 
                                 interval_length_bp, interval_length_kb, confidence_score, classification))
                
                
    info_sorted = sorted(info, key=lambda x: (sort_chrom(x[0]), x[1], x[2]))

    # 4. Write the BED file (name: "consensus_target_regions.bed")
    with open(os.path.join(dst, "consensus_target_regions.bed"), 'w') as ffo:
        ffo.write('#chrom\tstart\tend\tcnv_type\tcopy_num\tvariant_call_id\tvariant_region_id\t\
        validation\tnum_map\tdbVar_study_id\tnum_probe_CNV\tbp_overlap\tfraction_overlap\tinterval_length_bp\t\
        interval_length_kb\tconfidence_score\tclassification\n')
        ffo.write("\n".join(map(lambda x: "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
                                (x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], 
                                 x[12],x[13], x[14], x[15], x[16]), info_sorted)))
        
    print("'%s_%s_%s_merge.bed' file with the overlapping features created in %s..." %(id1, id2, id3, dst))
    print("'consensus_target_regions.bed' file with consensus CNV intervals created in %s... \n" %(dst))


def annotate(cnvfile, exonfile, dest):
    ''' 
    Annotate the BED consensus target regions file with gene names, based on overlap between the BED file CNV 
    intervals and coding exons (based on RefSeq hg19 exons (default setting) unless another exon/gene file is 
    indicated).

    Arguments
    ---------
        cnvfile: BED file with consensus CNV intervals.
        exonfile: File with exon coordinates for annotation (default: NCBI RefSeq exon list based on the
                GRCh37/hg19 assembly).
        dest: Initial output path 

    Returns
    -------
        - The 'ground_truth_global.bed' file with all the consensus CNV intervals, and their corresponding annotation
        - The 'ground_truth_exome.bed' file with the consensus CNV intervals that match in coding exon regions, and
        thir corresponding annotation.
        - The 'ground_truth_indeterminate.bed' file with the CNV intervals that cannot be assigned to a particular 
        region or type of alteration due to their high uncertainty and high inconsistency between initial analysis 
        studies.
        
        All files are saved in the specified output path.

    '''

    print("Annotating the consensus CNV intervals... ")
    
    # Load the file with the consensus CNVs and the file with exon coordinates to determine the number of fields
    # of each initial file
    with open(os.path.join(dest, cnvfile)) as cnv, open(exonfile) as exon:
        num_cols_cnv = len(next(csv.reader(cnv, delimiter='\t')))
        num_cols_exon = len(next(csv.reader(exon, delimiter='\t')))

    # If the exon/gene list is the default file, continue with the script.
    if exonfile == 'exome_INGEMM.no_pseudo.bed':
        fields = ['#chrom', 'start', 'end', 'cnv_type', 'copy_num', 'variant_call_id', 'variant_region_id', 
                  'validation', 'placements_per_assembly', 'dbVar_study_id', 'num_probe_CNV', 'interval_bp_overlap',
                  'fraction_overlap','interval_length_bp', 'interval_length_kb', 'confidence_score', 
                  'classification', 'gene_name', 'gene_ID','gene_number', 'bp_overlap_exon', 'count_exon']
        pass
    #Otherwise, check if the file contains at least 3 fields (chromosome, start and end) with chromosome coordinates.
    else:
        if num_cols_exon < 3:
            raise ValueError('''Error: The number of fields in the exon/gene list file is not correct. 
                             File must contain  at least chromosome, start and end separated by tabs. ''')
        else:
            fields = False

    # Create a BedTool object from the consensus CNV file
    ctr = pybedtools.BedTool(os.path.join(dest, cnvfile))

    # Use BedTools intersect to calculate the number of base pairs of overlap between the analysis features and
    # the matching exons (-wao); and the number of exons that match CNV intervals (-c).
    # Use BedTools groupby to summarized the dataset based upon common column groupings ('gr'). The 'co' list 
    # indicates the columns that should be summarized. The operation applied to 'co' is 'freqdesc', that prints
    # a separated list of values observed and the number of times they were observed, reported in 
    # descending order of frequency. 
    gr = list(range(1, num_cols_cnv+1))
    co = list(range(num_cols_cnv+4, num_cols_cnv+(num_cols_exon+3)))
    result = ctr.intersect(exonfile, wao=True).intersect(exonfile, c=True).groupby(g=gr, c=co, o=['freqdesc'])
    
    # Create a dataframe from the bedtool object 
    data = pybedtools.BedTool.to_dataframe(result, disable_auto_names=True, header=None)
    
    # If all placements per assembly are the same for a CNV interval, keep only one value
    data[8] = data[8].apply(lambda x: x.split('|')[0] if len(set(x.split('|'))) == 1 else x) #num_map
    # Keep only the number of exons that match the CNV interval
    data.iloc[:,-1] = data.iloc[:,-1].apply(lambda x: int(x.split(':')[0])) #exon count
    # If there is no overlap between CNV interval and exon, substitute the value by '0'
    data.iloc[:,-2] = data.iloc[:,-2].apply(lambda x: '0' if x == '0:1' else x) #bp_overlap_exon
    # If file is 'RefSeq_exon_pc.bed', when there is no overlap substitute gene_biotype, gene_id and gene_name by '.'
    data.iloc[:,-3] = data.iloc[:,-3].apply(lambda x: '.' if x == '.:1' else x) #if gene_number
    data.iloc[:,-4] = data.iloc[:,-4].apply(lambda x: '.' if x == '.:1' else x) #if gene_id
    data.iloc[:,-5] = data.iloc[:,-5].apply(lambda x: '.' if x == '.:1' else x) #if gene_name

    #1. Identify the consensus CNV interval classify as 'indeterminate' due to their uncertainty 
    indeterminates = (data[16] == 'Indeterminate')
    indet = data[indeterminates].copy(deep=True)
    
    #2. Identify the consensus CNV interval classify as 'HighConfidence' or 'LowConfidence'
    final = data[~indeterminates].copy(deep=True)
    
    #3. Identify the consensus CNV interval classify as 'HighConfidence' or 'LowConfidence', that match in 
    #coding exon regions
    exon = (data[16] != 'Indeterminate') & (data.iloc[:,-1] != 0)
    exome = data[exon].copy(deep=True)

    # Sort by chromosome and then by start and end position
    indet.sort_values(by=[0,1,2], key=lambda x: np.argsort(index_natsorted(indet[0])),inplace=True,ignore_index=True)
    final.sort_values(by=[0,1,2], key=lambda x: np.argsort(index_natsorted(final[0])),inplace=True,ignore_index=True)
    exome.sort_values(by=[0,1,2], key=lambda x: np.argsort(index_natsorted(exome[0])),inplace=True,ignore_index=True)

    # Write the BED files
    indet.to_csv(os.path.join(dest, "ground_truth_indeterminate.bed"), sep='\t', header=fields, index=False)
    final.to_csv(os.path.join(dest, "ground_truth_global.bed"), sep='\t', header=fields, index=False)
    exome.to_csv(os.path.join(dest, "ground_truth_exome.bed"), sep='\t', header=fields, index=False)
    
    
    print(''' 'ground_truth_global.bed' annotated file with consensus CNV intervals created in %s... ''' %(dest))
    print(''' 'ground_truth_exome.bed' annotated file with consensus CNV intervals matching exon regions \
    created in %s...''' %(dest))
    print(''' 'ground_truth_indeterminate.bed' annotated file with consensus CNV intervals with high uncertainty \
    created in %s... \n''' %(dest))

#----------------------------------------------------------------------------------------------------------------
### MAIN FUNCTION ###
#----------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='This script allows to obtain a set of consensus CNV intervals \
    which will determine the ground truth of CNVs, from the CNV calls characterized by three major CNV determination \
    benchmark studies (i.e. dbVar studies nstd46, estd20 and estd195 based on NA12878/HG001 sample). The resulting \
    CNV intervals have an associated confidence score, classification level and annotation.')
  
    parser.add_argument("-i1", "--infile1", required=True, help="Input file in BED format with the CNVs characterized \
                        by one of the benchmark studies after probe density analysis per CNV.\
                        This file should contain the following fields: 'chromosome', 'start', 'end', 'cnv type', \
                        'copy number', 'variant call ID', 'variant region ID', 'variant region accession', \
                        'validation', 'placements per assembly', 'number of probes'.",
                        dest="first_input")
    parser.add_argument("-id1", "--idinfile1", required=True, type=str, help="Study identifier or label referring to \
                        the study determined in the parameter '-i1'.")
    parser.add_argument("-exp1", "--experimental1", required=True, action="extend", nargs="+", type=str,
                        help="Experimental details referring to the study determined in the parameter '-i1' \
                        (i.e. BED file with the genomic coordinates of the probes used in the corresponding array).")

    parser.add_argument("-i2", "--infile2", required=True, help="Input file in BED format with the CNVs characterized \
                        by another of the benchmark studies after probe density analysis per CNV.\
                        This file should contain the following fields (the same as the first study provided): \
                        'chromosome', 'start', 'end', 'cnv type', 'copy number', 'variant call ID', 'variant region ID',\
                        'variant region accession', 'validation', 'placements per assembly', 'number of probes'.",
                        dest="second_input")
    parser.add_argument("-id2", "--idinfile2", required=True, type=str, help="Study identifier or label referring to \
                        the study determined in the parameter '-i2'.")
    parser.add_argument("-exp2", "--experimental2", required=True, action="extend", nargs="+", type=str,
                        help="Experimental details referring to the study determined in the parameter '-i2' \
                        (i.e. BED file with the genomic coordinates of the probes used in the corresponding array)")

    parser.add_argument("-i3", "--infile3", required=True, help="Input file in BED format with the CNVs characterized \
                        by another of the benchmark studies after probe density analysis per CNV.\
                        This file should contain the following fields (the same as the first study provided): \
                        'chromosome', 'start', 'end', 'cnv type', 'copy number', 'variant call ID', 'variant region ID',\
                        'variant region accession', 'validation', 'placements per assembly', 'number of probes'.", 
                        dest="third_input")
    parser.add_argument("-id3", "--idinfile3", required=True, type=str, help="Study identifier or label referring to \
                        the study determined in the parameter '-i3'.")
    parser.add_argument("-exp3", "--experimental3", required=True, action="extend", nargs="+", type=str,
                        help="Experimental details referring to the study determined in the parameter '-i3' \
                        (i.e. BED file with the genomic coordinates of the probes used in the corresponding array)")

    parser.add_argument("-o", "--outpath", required=True, help="Output path to save the files generated during \
                        script execution", dest="outputpath")

    parser.add_argument("-e", "--exon", default="exome_INGEMM.no_pseudo.bed", required=False, help="File with exon \
                        coordinates for annotation. \
                        Default file: NCBI RefSeq exon list based on GRCh37/hg19 assembly. \
                        Another exon list may be provided.")

    args = parser.parse_args()

    print('''\n#-----------------------------------------------------------------------------------------------------
        Starting python script and checking parameters... 
#----------------------------------------------------------------------------------------------------- \n''')

    #Check first study file exists
    first_study = args.first_input
    if not os.path.exists(first_study):
        raise IOError('File %s does not exist.' % (first_study))
    #Check experimental design file corresponding to first study exists
    for f in args.experimental1:
        if not os.path.exists(f):
            raise IOError('The experimental design file %s does not exist.' % (f))
    
    #Check second study file exists
    second_study = args.second_input
    if not os.path.exists(second_study):
        raise IOError('File %s does not exist.' % (second_study))
    #Check experimental design file corresponding to second study exists
    for f in args.experimental2:
        if not os.path.exists(f):
            raise IOError('The experimental design file %s does not exist.' % (f))
            
    #Check third study file exists
    third_study = args.third_input
    if not os.path.exists(third_study):
        raise IOError('File %s does not exist.' % (third_study))
    #Check experimental design file corresponding to third study exists
    for f in args.experimental3:
        if not os.path.exists(f):
            raise IOError('The experimental design file %s does not exist.' % (f))

    #Check that the output path exists, if not create it
    outgroundtruth = args.outputpath
    if not os.path.exists(outgroundtruth):
        os.makedirs(outgroundtruth)

    ##Calling the required functions:
    new_version_bedfile(first_study, args.idinfile1, outgroundtruth)
    new_version_bedfile(second_study, args.idinfile2, outgroundtruth)
    new_version_bedfile(third_study, args.idinfile3, outgroundtruth)

    consensus_cnv(first_study, args.idinfile1, args.experimental1,
                  second_study, args.idinfile2, args.experimental2,
                  third_study, args.idinfile3, args.experimental3,
                  outgroundtruth)

    annotate('consensus_target_regions.bed', args.exon, outgroundtruth)

##################################################################################################################

if __name__ == '__main__':
    try:
        main()
        print('''\n     Script executed \n''')
    except Exception as e:
        stderr.write(f"Error: {e} \n")
        exit()

