import allel
import os
import pandas as pd

# Function to read in the vcf exactly as it is into a dataframe.
def vcf_body_to_df(file_path):
    with open(file_path, 'r') as vcf_text_wrapper:
        vcf_lines = vcf_text_wrapper.readlines()

    vcf_labels_preprocessed = [line for line in vcf_lines if line.startswith("#") and not(line.startswith("##"))]

    for line in vcf_labels_preprocessed:
        line_list = line.split('\n')
        vcf_label_split = [l.split('\t') for l in line_list if len(l) > 0]

    # flatten list...
    vcf_labels = []
    for label in vcf_label_split:
        for item in label:
            vcf_labels.append(item)

    vcf_body_preprocessed = [line for line in vcf_lines if not line.startswith("#")]

    vcf_body_split = []
    for line in vcf_body_preprocessed:
        line_list = line.split('\n')
        vcf_body_split.append([l.split('\t') for l in line_list if len(l) > 0])

# flatten list...
    vcf_body = []
    for entry in vcf_body_split:
        for item in entry:
            vcf_body.append(item)

    df = pd.DataFrame.from_records(vcf_body, columns=vcf_labels)

    return df

def split_multigene_vcf(file_path):
    file_list = os.listdir(file_path)

    for vcf_name in file_list:
        name_root, name_ext = os.path.splitext(vcf_name)
        if name_ext == ".vcf":
            patient_id = name_root.split('-')[0]
            #TODO: temporary string to replace build number...
            genome_build = 'b37'
            path_to_vcf = file_path+vcf_name
            header_lines = allel.read_vcf_headers(path_to_vcf).headers
            '''
            # Fix to make single gene VCF files compatible with elimu tool
            header_lines_ps = header_lines[:]
            header_lines_ps[len(header_lines)-1] = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
            header_lines_ps.append(header_lines[len(header_lines)-1])
            '''
            gi_df = allel.vcf_to_dataframe(path_to_vcf, fields='GI', alt_number=1)
            unique_GI = set(gi_df['GI'])
            vcf_body_complete_df = vcf_body_to_df(path_to_vcf)
            vcf_info_field = vcf_body_complete_df.INFO
            vcf_info_split = [el.split(';') for el in vcf_info_field]
            vcf_info_gi = []
            for row in vcf_info_split:
                for item in row:
                    if item.startswith('GI'):
                        vcf_gi = item.split('=')[1]
                        vcf_info_gi.append(vcf_gi)
                        break
            for gi in unique_GI:
                unique_GI_vars_df = vcf_body_complete_df[[x == gi for x in vcf_info_gi]]
                new_file_name = patient_id+'.'+genome_build+'.'+gi+name_ext
                f = open(file_path+'single_gene_vcf/'+new_file_name, 'w+')
                f.writelines(header_lines)
                vcf_file_values = unique_GI_vars_df.values
                '''
                # Fix for charite vcf files to be compatible with elimu vcf2fhir tool!!
                vcf_file_values[:,8] = vcf_file_values[:,8] + ':PS'
                vcf_file_values[:,9] = vcf_file_values[:,9] + ':0'
                '''
                for item in vcf_file_values:
                    line = '\t'.join(item)+'\n'
                    f.write(line)
                f.close()

if __name__ == '__main__':
    file_path = "./charite_example_files/"
    split_multigene_vcf(file_path)