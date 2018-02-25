import sys
import getopt
import glob
import re
import time

def main(argv):
    inFile = ''
    koFullIn = ''
    annoIn = ''
    try:
        opts, args = getopt.getopt(argv,'i:k:a:h',['in','koFull','anno','help'])
    except getopt.GetoptError:
        print('python analyzeHamfasOnly.py -i input.list -k list.KO.FULL -a annoDir/pfam.xml')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('h','--help'):
            print('python analyzeHamfasOnly.py -i input.list -k list.KO.FULL -a annoDir/pfam.xml')
            sys.exit()
        elif opt in ('-i','--in'):
            inFile = arg
        elif opt in ('-k','--koFull'):
            koFullIn = arg
        elif opt in ('-a','--anno'):
            annoIn = arg

    anno = ''
    with open(annoIn, 'r') as annoFile:
        anno = annoFile.read().replace('\n', '')

    ### get type of reference species where KO annotation come from
    type = {}
    type['ape']='3'
    type['mja']='3'

    type['aae']='4'
    type['bsu']='4'
    type['eco']='4'
    type['hpy']='4'
    type['lla']='4'
    type['mge']='4'
    type['mtu']='4'
    type['nme']='4'
    type['syn']='4'

    type['sacce']='1'
    type['ago']='1'
    type['aspni']='1'
    type['canal']='1'
    type['neucr']='1'
    type['schpo']='1'

    type['cel']='2'
    type['dme']='2'
    type['dre']='2'
    type['hsa']='2'
    type['mmu']='2'
    type['nemve']='2'
    type['rno']='2'
    type['monbr']='2'
    type['cho']='2'
    type['pfa']='2'
    type['trybr']='2'
    type['ath']='2'
    type['ehi']='2'

    koType = {}     # koType{geneID:koID:type}
    koTypeMin = {}  # koTypeMin{genID:minType}

    fungiOrtho = {} # fungiOrtho{geneID:Y/N}
    fungalList = ['ago','aspni','canal','neucr','schpo']

    with open(koFullIn, 'r') as koFullFile:
        ### split into multiple orthogoups
        koFullGroup = re.split(r'### ',koFullFile.read())
        for i in range(1,len(koFullGroup)):
            ## split group into lines
            orthologs = re.split(r'\n',koFullGroup[i])

            sacceID = orthologs[1].rstrip('\t')
            koType[sacceID] = {}

            ## check if this gene has fungal ortholog
            fungiOrtho[sacceID] = "N"
            for fungus in fungalList:
                if re.search(fungus,koFullGroup[i]):
                    fungiOrtho[sacceID] = "Y"
                    break

            ## process only group that have KO annotated
            if re.search('K\d{5};',orthologs[0]):
                for ortho in orthologs:
                    ## get orthologs that have KO
                    if re.search('K\d{5}\t',ortho):
                        # get KO id
                        orthoKOid = re.split('\t',ortho)[1]
                        # check if this KO ID was transfered to the whole group
                        if re.search(orthoKOid,orthologs[0]):
                            # get reference species type
                            refSpec = re.split(':',ortho)[0]
                            if not orthoKOid in koType[sacceID]:
                                koType[sacceID][orthoKOid] = type[refSpec]
                            else:
                                if koType[sacceID][orthoKOid] > type[refSpec]:
                                    koType[sacceID][orthoKOid] = type[refSpec]

                            if not sacceID in koTypeMin:
                                koTypeMin[sacceID] = type[refSpec]
                            else:
                                if koTypeMin[sacceID] > type[refSpec]:
                                    koTypeMin[sacceID] = type[refSpec]

    ### read input list
    with open(inFile) as f:
        for line in f:
            ## get gene ID
            geneID = re.split(r'\t+',line.rstrip())[0]
            #print(geneID)

            ## get number of pfam domains annotated for this protein
            annoPattern = re.compile('protein id=\"'+geneID+'\".+?<\/protein')
            matchAnno = annoPattern.search(anno)
            count = 0
            if(matchAnno):
                count = len(re.findall('feature type',matchAnno.group(0)))
                print(count)
            else:
                print(geneID + " noAnno")

            ### print OUTPUT
            output = [geneID,str(count),"t"+str(koTypeMin[geneID]),fungiOrtho[geneID]]
#            for ko in koType[geneID].keys():
#                output.append(ko+"#"+koType[geneID][ko])
            # print(output)
#            print('\t'.join(output))
            time.sleep(1)

if __name__ == "__main__":
    if len(sys.argv[1:]) < 6:
        print('python analyzeHamfasOnly.py -i input.list -k list.KO.FULL -a annoDir/pfam.xml')
        sys.exit(2)
    else:
        main(sys.argv[1:])
