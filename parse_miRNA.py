import pandas as pd

dict_miRNA_target_mirTarbase = {}
dict_miRNA_target_msigDB = {}
dict_miRNA_target = {}

mirTarbase_file = "/Users/guangrong/Documents/GitHub_project/fastqpi_BigGIM/KGs/miRTarBase/hsa_MTI.csv"

mirTarbase = pd.read_csv(mirTarbase_file)

def format(miRNA):
    miRNA = miRNA.replace("hsa-miR-", "MIR")
    miRNA = miRNA.replace("hsa-let-", "LET")
    #print(miRNA)
    miRNA = miRNA.replace("-", "_")
    miRNA = miRNA.upper()
    return(miRNA)

new_id = []
for miRNA in list(set(mirTarbase['miRNA'])):
    new_id.append(format(miRNA))

print("Number of unique miRNA in miRtarBase: " + str(len(new_id)))

for i in range(0,mirTarbase.shape[0]):
    miRNA = format(mirTarbase.iloc[i,1])
    if miRNA not in dict_miRNA_target_mirTarbase:
        dict_miRNA_target_mirTarbase[miRNA] = set([mirTarbase.iloc[i,3]])
    else:
        dict_miRNA_target_mirTarbase[miRNA].add(mirTarbase.iloc[i,3])

msigDB_miRNA_file = "/Users/guangrong/Documents/GitHub_project/fastqpi_BigGIM/KGs/MsigDB/v7.4/c3.mir.v7.4.symbols.gmt"
fin = open(msigDB_miRNA_file,'r')
mirna_list = []
target_list = []
for line in fin.readlines():
    words = line.strip().split('\t')
    
    for i in range(2, len(words)):
        mirna_list.append(words[0])
        target_list.append(words[i])
        if words[0] not in dict_miRNA_target_msigDB:
            dict_miRNA_target_msigDB[words[0]] = set(words[i])
        else:
            dict_miRNA_target_msigDB[words[0]].add(words[i])

fin.close()


print("Number of unique miRNA in MsigDB: " + str(len(set(mirna_list))))
print("Number of unique target genes in MsigDB: " + str( len(set(target_list))))
print( "Number of intersection miRNAs between two resources: " + str(len(set(new_id).intersection(set(mirna_list)))))

print(set(mirna_list) - set(new_id))

ALL_miRNA = list(set(list(dict_miRNA_target_msigDB.keys()) + list(dict_miRNA_target_mirTarbase.keys())))

Evidence_list = []
mirna_list_all = []
target_list_all = []
for mi in ALL_miRNA:
    #print(mi)
    if mi in dict_miRNA_target_msigDB and mi in dict_miRNA_target_mirTarbase:
        targets = set(list(dict_miRNA_target_msigDB[mi]) + list(dict_miRNA_target_mirTarbase[mi])) 
        for target in targets:
            mirna_list_all.append(mi)
            target_list_all.append(target)
            
            Evidence_list.append('MSigDBv7.4_and_miRTarbase')
            
    elif mi in dict_miRNA_target_msigDB and mi not in dict_miRNA_target_mirTarbase:
        for target in dict_miRNA_target_msigDB[mi]:
            Evidence_list.append('MSigDBv7.4')
            mirna_list_all.append(mi)
            target_list_all.append(target)

    else:
        for target in dict_miRNA_target_mirTarbase[mi]:
            
            mirna_list_all.append(mi)
            target_list_all.append(target)
            Evidence_list.append('miRTarbase')
print(len(mirna_list_all))
print(len(target_list_all))
result = pd.DataFrame({"Subject": mirna_list_all, 
                        "Subject_category": ['miRNA'] * len(mirna_list_all),
                        "Object":target_list_all, 
                        "Object_category": ['Gene'] * len(target_list_all), 
                        "predicates":['miRNA_target']*len(mirna_list_all),
                        "Evidence": Evidence_list})
result.to_csv("KG_miRNA_target_edges.csv", index = False)

print("Number of unique miRNA inm mirna_list_all: " + str(len(set(mirna_list_all))))
print("Number of unique miRNA inm target_list_all: " + str(len(set(target_list_all))))

