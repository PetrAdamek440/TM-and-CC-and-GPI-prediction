
import time

t = time.time()

import os
#import gc
import glob

#
filesPNG = glob.glob('figures/*.png') # , recursive=True)

for fpng in filesPNG:
    try:
        os.remove(fpng)
    except OSEError as epng:
        print("Error: %s : %s" % (fpng,epng.strerror))

filesSVG = glob.glob('figures/*.svg') # , recursive=True)

for fsvg in filesSVG:
    try:
        os.remove(fsvg)
    except OSEError as esvg:
        print("Error: %s : %s" % (fsvg,esvg.strerror))
#
#
filesPNGg = glob.glob('figures/given/*.png') # , recursive=True)

for fpngg in filesPNGg:
    try:
        os.remove(fpngg)
    except OSEError as epngg:
        print("Error: %s : %s" % (fpng,epngg.strerror))

filesSVGg = glob.glob('figures/given/*.svg') # , recursive=True)
#
#
for fsvgg in filesSVGg:
    try:
        os.remove(fsvgg)
    except OSEError as esvgg:
        print("Error: %s : %s" % (fsvg,esvgg.strerror))

filesPNGc = glob.glob('figures/complementary/*.png') # , recursive=True)

for fpngc in filesPNGc:
    try:
        os.remove(fpngc)
    except OSEError as epngc:
        print("Error: %s : %s" % (fpng,epngc.strerror))

filesSVGc = glob.glob('figures/complementary/*.svg') # , recursive=True)
#
#
for fsvgc in filesSVGc:
    try:
        os.remove(fsvgc)
    except OSEError as esvgc:
        print("Error: %s : %s" % (fsvg,esvgc.strerror))

filesCCp = glob.glob('CC_plots/plotsP/*.png') # , recursive=True)

for fCCp in filesCCp:
    try:
        os.remove(fCCp)
    except OSEError as e:
        print("Error: %s : %s" % (f,e.strerror))

filesCCn = glob.glob('CC_plots/plotsN/*.png') # , recursive=True)
#
#
for fCCn in filesCCn:
    try:
        os.remove(fCCn)
    except OSEError as e:
        print("Error: %s : %s" % (f,e.strerror))

filesCCpR = glob.glob('CC_plots/plotsP-realCC/*.png') # , recursive=True)

for fCCpR in filesCCpR:
    try:
        os.remove(fCCpR)
    except OSEError as e:
        print("Error: %s : %s" % (f,e.strerror))

filesCCnR = glob.glob('CC_plots/plotsN-realCC/*.png') # , recursive=True)

for fCCnR in filesCCnR:
    try:
        os.remove(fCCnR)
    except OSEError as e:
        print("Error: %s : %s" % (f,e.strerror))
#
#
filesCCpL = glob.glob('CC_plots/plotsP-realCC/line-picture-CC-P/*.png') # , recursive=True)

for fCCpL in filesCCpL:
    try:
        os.remove(fCCpL)
    except OSEError as e:
        print("Error: %s : %s" % (f,e.strerror))

filesCCnL = glob.glob('CC_plots/plotsN-realCC/line-picture-CC-N/*.png') # , recursive=True)

for fCCnL in filesCCnL:
    try:
        os.remove(fCCnL)
    except OSEError as e:
        print("Error: %s : %s" % (f,e.strerror))
#
#
filesFASTAp = glob.glob('FASTA_output/peptidesP/*.*') # , recursive=True)

for fFASTAp in filesFASTAp:
    try:
        os.remove(fFASTAp)
    except OSEError as e:
        print("Error: %s : %s" % (f,e.strerror))

filesFASTAn = glob.glob('FASTA_output/peptidesN/*.*') # , recursive=True)

for fFASTAn in filesFASTAn:
    try:
        os.remove(fFASTAn)
    except OSEError as e:
        print("Error: %s : %s" % (f,e.strerror))
#
#

filesGPIp = glob.glob('figures/GPI/P/*.png')

for fGPIp in filesGPIp:
    try:
        os.remove(fGPIp)
    except OSEError as eg:
        print("Error: %s : %s" % (f,eg.strerror))

filesGPIn = glob.glob('figures/GPI/N/*.png')

for fGPIn in filesGPIn:
    try:
        os.remove(fGPIn)
    except OSEError as eg:
        print("Error: %s : %s" % (f,eg.strerror))

#

vstupni_sekvence = input('Zadejte cestu k souboru s DNA sekvenci pro prohlednuti: ')
nucleotides = open(vstupni_sekvence)

#nucleotides = open('sekvence_zadavane/rand_string.txt') # test random sekvence

#nucleotides = open('sekvence_zadavane/Human_tetherin_locus (1500bp).txt') # test delky lokusu
#nucleotides = open('sekvence_zadavane/Human_tetherin_locus (3000bp).txt') # test delky lokusu
#nucleotides = open('sekvence_zadavane/Human_tetherin_locus (30kbp).txt') # test delky lokusu

#nucleotides = open('sekvence_zadavane/Human_tetherin_locus.txt') # posledni

#nucleotides = open('sekvence_zadavane/DP/Human_tetherin_locus-vDP.txt') # posledni - k DP
#nucleotides = open('sekvence_zadavane/DP/Chicken_tetherin_locus_v2.txt') # posledni - k DP
#nucleotides = open('sekvence_zadavane/DP/Mouse_tetherin_locus.txt') # posledni - k DP
#nucleotides = open('sekvence_zadavane/DP/hg38_dna-chr16-NONannotated-vDEvPAv1.txt') # posledni - k DP
#nucleotides = open('sekvence_zadavane/DP/rand_string-2021-04-21.txt') # posledni - k DP
#nucleotides = open('sekvence_zadavane/DP/rand/rand_string-v04.txt') # posledni - k DP
#nucleotides = open('sekvence_zadavane/DP/Manis_javanica_tetherin_locus-ORIZNUTO.txt') # posledni - k DP

#nucleotides = open('sekvence_zadavane/DP/Pteropus_alecto_tetherin.txt') # posledni - k DP

#nucleotides = open('sekvence_zadavane/Chicken_tetherin_locus.txt') # posledni
#nucleotides = open('sekvence_zadavane/zaba_Xenopus.txt') # pred-posledni

dolni_limit = input('Zadejte dolni limit povolene delky vstupni sekvence v bp: ') # 320 # bazi
dolni_limit = int(dolni_limit)
horni_limit = input('Zadejte horni limit povolene delky vstupni sekvence v bp: ') # 32000 # bazi
horni_limit = int(horni_limit)

CC_threshold_P = input('Zadejte prah pro predikci coiled-coil struktury v zadanem vlakne (napr. 0.10 (tj 10 %)): ')
CC_threshold_P = float(CC_threshold_P)

CC_threshold_N = input('Zadejte prah pro predikci coiled-coil struktury v komplementarnim vlakne (napr. 0.10 (tj 10 %)): ')
CC_threshold_N = float(CC_threshold_N)

bases = nucleotides.read()

basesRP = bases.replace('\n','')

basesP = basesRP

lenbP = len(basesP)
#print("Nukleotidu (delka NA) je (po redukci o pripadny znak noveho radku) = ",lenbP)

basesN1 = bases.replace('A','%temp%').replace('T','A').replace('%temp%','T')
basesN2 = basesN1.replace('G','%temp%').replace('C','G').replace('%temp%','C')
basesN35 = basesN2
basesN53 = basesN35[len(basesN35)::-1]
basesNop = basesN53

basesRN = basesNop.replace('\n','')

basesN = basesRN

lenbN = len(basesN)

if (lenbP < dolni_limit):
    print('Delka vstupni sekvence je pod dolnim povolenym limitem ',dolni_limit,' paru bazi.')
elif (lenbP > horni_limit):
    print('Delka vstupni sekvence je nad hornim povolenym limitem, ',horni_limit,' paru bazi.')
else:

    limitORFzadany = input('Zadejte minimalni delku ORFu v bp (napr. 150): ') # napr. 150 # DULEZITY nastavitelny parametr # Je vcetne start a stop kodonu.
    limitORFzadany = float(limitORFzadany)
    
    limitORF = abs(round(limitORFzadany)) # DULEZITY nastavitelny parametr # Je vcetne start a stop kodonu.

    print('\n')

    print('Pouzita minimalni delka ORFu je: ', limitORF,' bp (vstupni udaj v absolutni hodnote zaokrouhleny na cele cislo).')

    if ((limitORF < 30) | (limitORF > lenbP) | (limitORF > lenbN)) == True:
        print('Delka ORFu musi byt minimalne 30 bp a maximalne delka zadaneho vlakna (v bp).')
    else:        
        limitORFP = limitORF # DULEZITY nastavitelny parametr
    
        print('\n')

        #print('Vypis zadane sekvence:')
        #print(bases)
        #print('\n')
        #print('Vypis komplementarni sekvence:')
        #print(basesN)
        #print('\n')

        #lenb = len(bases)
        #print("Nukleotidu (delka NA) je = ",lenb)
        #print('\n')

        import re

        # vyhledani ORFu (v zadanem vlakne)
        # def.2 podle :
        '''
        Def. 2: ORF je usek DNA mezi dvemi STOP kodony a pocet jeho nukleotidu je deliteny 3 beze zbytku
        '''
        stop_codon_findP = [mPstop.start() for mPstop in re.finditer('(TAA|TAG|TGA)', basesP)]
        #print('Stop kodon v zadanem vlakne je na techto pozicich: \n',stop_codon_findP)
        #print('Stop kodon je tolikrat:',len(stop_codon_findP))
        ## najde opravdu, kde zacina dana sekvence, tedy kde je "T"
        ## vyhodi cislo, kde je prave to pismeno, napr. je-li "T" 6. zleva pak: "5"

        varianty_ambiguitni_kod = input('Osetreni ambiguitniho kodu, nezname znaky - pro nahodne nahrazeni aminokyselinou: 1, pro vystrizeni znaku: 2, pro vyrazeni celeho ORFu s nezm. zn.: 3--> zadejte volbu: ')

        proces123 = int(varianty_ambiguitni_kod)

        if (proces123 == 1):

            stop_kodon_1_ramec_P = []
            stop_kodon_2_ramec_P = []
            stop_kodon_3_ramec_P = []

            listORFsP = []
            #listORFsN = [] ##### - je nize ... 

            listORFsP_1 = []
            listORFsP_2 = []
            listORFsP_3 = []

            listORFsP_1_pozice = []
            listORFsP_2_pozice = []
            listORFsP_3_pozice = []

            for ixyzP in stop_codon_findP:
                if ((ixyzP)%3 == 0):
                    stop_kodon_1_ramec_P += [ixyzP]
                elif ((ixyzP+2)%3 == 0):
                    stop_kodon_2_ramec_P += [ixyzP]
                elif ((ixyzP+1)%3 == 0):
                    stop_kodon_3_ramec_P += [ixyzP]

            #print('Zadane vlakno - stop_kodon_1_ramec: \n',stop_kodon_1_ramec_P)
            #print('Zadane vlakno - stop_kodon_2_ramec: \n',stop_kodon_2_ramec_P)
            #print('Zadane vlakno - stop_kodon_3_ramec: \n',stop_kodon_3_ramec_P)


            if (len(stop_kodon_1_ramec_P) < 2):
                print('Nedostecny pocet STOP kodonu v 1. ramci v zadanem vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix1P in range(0,len(stop_kodon_1_ramec_P)-1,1):
                    if ((stop_kodon_1_ramec_P[ix1P+1] - stop_kodon_1_ramec_P[ix1P]) < limitORFP):
                        continue
                    elif ((stop_kodon_1_ramec_P[ix1P+1] - stop_kodon_1_ramec_P[ix1P]) >= limitORFP):
                        listORFsP_1 += [basesP[(stop_kodon_1_ramec_P[ix1P]+3):(stop_kodon_1_ramec_P[ix1P+1]+3)]]
                        listORFsP_1_pozice += [[stop_kodon_1_ramec_P[ix1P], stop_kodon_1_ramec_P[ix1P+1]]]
                #print(listORFsP_1)
                #print('Pozice ORFu v zadanem vlakne v 1. ramci jsou (v kodu: "listORFsP_1_pozice"): \n',listORFsP_1_pozice)

            if (len(stop_kodon_2_ramec_P) < 2):
                print('Nedostecny pocet STOP kodonu v 2. ramci v zadanem vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix2P in range(0,len(stop_kodon_2_ramec_P)-1,1):
                    if ((stop_kodon_2_ramec_P[ix2P+1] - stop_kodon_2_ramec_P[ix2P]) < limitORFP):
                        continue
                    elif ((stop_kodon_2_ramec_P[ix2P+1] - stop_kodon_2_ramec_P[ix2P]) >= limitORFP):
                        listORFsP_2 += [basesP[(stop_kodon_2_ramec_P[ix2P]+3):(stop_kodon_2_ramec_P[ix2P+1]+3)]]
                        listORFsP_2_pozice += [[stop_kodon_2_ramec_P[ix2P], stop_kodon_2_ramec_P[ix2P+1]]]
                #print(listORFsP_2)
                #print('Pozice ORFu v zadanem vlakne ve 2. ramci jsou (v kodu: "listORFsP_2_pozice"): \n',listORFsP_2_pozice)            

            if (len(stop_kodon_3_ramec_P) < 2):
                print('Nedostecny pocet STOP kodonu v 3. ramci v zadanem vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix3P in range(0,len(stop_kodon_3_ramec_P)-1,1):
                    if ((stop_kodon_3_ramec_P[ix3P+1] - stop_kodon_3_ramec_P[ix3P]) < limitORFP):
                        continue
                    elif ((stop_kodon_3_ramec_P[ix3P+1] - stop_kodon_3_ramec_P[ix3P]) >= limitORFP):
                        listORFsP_3 += [basesP[(stop_kodon_3_ramec_P[ix3P]+3):(stop_kodon_3_ramec_P[ix3P+1]+3)]]
                        listORFsP_3_pozice += [[stop_kodon_3_ramec_P[ix3P], stop_kodon_3_ramec_P[ix3P+1]]]
                #print(listORFsP_3)
                #print('Pozice ORFu v zadanem vlakne ve 3. ramci jsou (v kodu: "listORFsP_3_pozice"): \n',listORFsP_3_pozice)     

            listORFsP = listORFsP_1 + listORFsP_2 + listORFsP_3
            #print('Pocet ORFu pro min. delku ORFu ',limitORFP,'je: ',len(listORFsP))
            #print('Seznam ORFu (pro "P"), tj. zadane vlakno: ',listORFsP)

            listORFsP_pozice = listORFsP_1_pozice + listORFsP_2_pozice + listORFsP_3_pozice ## !!! ve skutecnosti se hodi <-- toto
            #sorted_listORFsP_pozice = sorted(listORFsP_pozice) # bude se hodit pro vyhledavani pozice pouziteho ORFu pro vysetrovany peptid/protein

            '''
            for iP in range(LP):
                #listBP = list(refallBasesP[iP][:])
                listBP2x = list(refallBasesP[iP][:])
                listBP = list(listBP2x[:][1:])
                sP = ''
                seqP = listBP
                JseqP = sP.join(seqP)
                if (len(JseqP))%3 != 0:
                 #try:
                    #JseqP = ''
                    JseqP = None
                elif (((len(JseqP))%3 == 0 )) & ((len(JseqP)) >= limitORFP):
                    #print(JseqP)
                    countP += 1
                    listORFsP += [JseqP]
                    continue
                else:
                    JseqP = None
                    #break

            #if countP:
            #    print('pocet ORFu v zadanem vlaknu je: ', countP)
            '''
            '''
            basesN1 = bases.replace('A','%temp%').replace('T','A').replace('%temp%','T')
            basesN2 = basesN1.replace('G','%temp%').replace('C','G').replace('%temp%','C')
            basesN35 = basesN2
            basesN53 = basesN35[len(basesN35)::-1]
            basesN = basesN53

            basesRN = basesN.replace('\n','')

            basesN = basesRN

            lenbN = len(basesN)
            '''

            ####print("Nukleotidu (delka NA) je po redukci o znaky pro novy radek: ",lenbN)

            #####
            #print('Vypis komplementarni sekvence k zadane sekvenci:')
            #print(basesN)
            #####

            # vyhledani ORFu (v komplementarnim vlakne)
            # def.2 podle:
            '''
            Def. 2: ORF je usek DNA mezi dvemi STOP kodony a pocet jeho nukleotidu je deliteny 3 beze zbytku
            '''
            stop_codon_findN = [mNstop.start() for mNstop in re.finditer('(TAA|TAG|TGA)', basesN)]
            #print('Stop kodon v komplementarnim vlakne je na techto pozicich: \n',stop_codon_findN)
            #print('Stop kodon je tolikrat:',len(stop_codon_findN))
            ## najde opravdu, kde zacina dana sekvence, tedy kde je "T"
            ## vyhodi cislo, kde je prave to pismeno, napr. je-li "T" 6. zleva pak: "5"

            limitORFN = limitORF # DULEZITY nastavitelny parametr

            stop_kodon_1_ramec_N = []
            stop_kodon_2_ramec_N = []
            stop_kodon_3_ramec_N = []

            listORFsN = []

            listORFsN_1 = []
            listORFsN_2 = []
            listORFsN_3 = []

            listORFsN_1_pozice = []
            listORFsN_2_pozice = []
            listORFsN_3_pozice = []

            for ixyzN in stop_codon_findN:
                if ((ixyzN)%3 == 0):
                    stop_kodon_1_ramec_N += [ixyzN]
                elif ((ixyzN+2)%3 == 0):
                    stop_kodon_2_ramec_N += [ixyzN]
                elif ((ixyzN+1)%3 == 0):
                    stop_kodon_3_ramec_N += [ixyzN]

            #print('stop_kodon_1_ramec: \n',stop_kodon_1_ramec_N)
            #print('stop_kodon_2_ramec: \n',stop_kodon_2_ramec_N)
            #print('stop_kodon_3_ramec: \n',stop_kodon_3_ramec_N)


            if (len(stop_kodon_1_ramec_N) < 2):
                print('Nedostecny pocet STOP kodonu v 1. ramci v komplemtarnim vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix1N in range(0,len(stop_kodon_1_ramec_N)-1,1):
                    if ((stop_kodon_1_ramec_N[ix1N+1] - stop_kodon_1_ramec_N[ix1N]) < limitORFN):
                        continue
                    elif ((stop_kodon_1_ramec_N[ix1N+1] - stop_kodon_1_ramec_N[ix1N]) >= limitORFN):
                        listORFsN_1 += [basesN[(stop_kodon_1_ramec_N[ix1N]+3):(stop_kodon_1_ramec_N[ix1N+1]+3)]]
                        listORFsN_1_pozice += [[stop_kodon_1_ramec_N[ix1N], stop_kodon_1_ramec_N[ix1N+1]]]
                #print(listORFsN_1)
                #print('Pozice ORFu v komplementarnim vlakne v 1. ramci jsou (v kodu: "listORFsN_1_pozice"): \n',listORFsN_1_pozice)

            if (len(stop_kodon_2_ramec_N) < 2):
                print('Nedostecny pocet STOP kodonu v 2. ramci v komplemtarnim vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix2N in range(0,len(stop_kodon_2_ramec_N)-1,1):
                    if ((stop_kodon_2_ramec_N[ix2N+1] - stop_kodon_2_ramec_N[ix2N]) < limitORFN):
                        continue
                    elif ((stop_kodon_2_ramec_N[ix2N+1] - stop_kodon_2_ramec_N[ix2N]) >= limitORFN):
                        listORFsN_2 += [basesN[(stop_kodon_2_ramec_N[ix2N]+3):(stop_kodon_2_ramec_N[ix2N+1]+3)]]
                        listORFsN_2_pozice += [[stop_kodon_2_ramec_N[ix2N], stop_kodon_2_ramec_N[ix2N+1]]]
                #print(listORFsN_2)
                #print('Pozice ORFu v komplementarnim vlakne ve 2. ramci jsou (v kodu: "listORFsN_2_pozice"): \n',listORFsN_2_pozice)            

            if (len(stop_kodon_3_ramec_N) < 2):
                print('Nedostecny pocet STOP kodonu v 3. ramci v komplemtarnim vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix3N in range(0,len(stop_kodon_3_ramec_N)-1,1):
                    if ((stop_kodon_3_ramec_N[ix3N+1] - stop_kodon_3_ramec_N[ix3N]) < limitORFN):
                        continue
                    elif ((stop_kodon_3_ramec_N[ix3N+1] - stop_kodon_3_ramec_N[ix3N]) >= limitORFN):
                        listORFsN_3 += [basesN[(stop_kodon_3_ramec_N[ix3N]+3):(stop_kodon_3_ramec_N[ix3N+1]+3)]]
                        listORFsN_3_pozice += [[stop_kodon_3_ramec_N[ix3N], stop_kodon_3_ramec_N[ix3N+1]]]
                #print(listORFsN_3)
                #print('Pozice ORFu v komplementarnim vlakne ve 3. ramci jsou (v kodu: "listORFsN_3_pozice"): \n',listORFsN_3_pozice)    

            listORFsN = listORFsN_1 + listORFsN_2 + listORFsN_3
            #print('Pocet ORFu pro min. delku ORFu ',limitORFN,'je: ',len(listORFsN))
            #print('Seznam ORFu (pro "N"), tj. komplemtarni vlakno: ',listORFsN)

            listORFsN_pozice = listORFsN_1_pozice + listORFsN_2_pozice + listORFsN_3_pozice ## !!! ve skutecnosti se hodi <-- toto
            #sorted_listORFsP_pozice = sorted(listORFsP_pozice) # bude se hodit pro vyhledavani pozice pouziteho ORFu pro vysetrovany peptid/protein

            # Nyni preklad do peptidu:

            # Standardni geneticky kod:
            Base1 = ['TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG']
            Base2 = ['TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG']
            Base3 = ['TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG']
            Standard_Code = ['FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']

            Standard_Code_nonSTOP = ['FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            # je Standard_Code bez STOP kodonu, pro preklenuti neznamych znaku v DNA sekvenci
            # je pouzit napr. v radku s: "AMKP_c00 = random.choice(Standard_Code_nonSTOP[0])"

            # Alternativni geneticke kody:
            Vertebrate_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG']
            Yeast_Mitochondrial_Code = ['FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Mold_Protozoan_Coelenterate_Mitochondrial_and_Mycoplasma_Spiroplasma_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Invertebrate_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG']
            Ciliate_Dasycladacean_and_Hexamita_Nuclear_Code = ['FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Echinoderm_and_Flatworm_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG']
            Euplotid_Nuclear_Code = ['FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Bacterial_Archaeal_and_Plant_Plastid_Code = ['FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Alternative_Yeast_Nuclear_Code = ['FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Ascidian_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG']
            Alternative_Flatworm_Mitochondrial_Code = ['FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG']
            Blepharisma_Nuclear_Code = ['FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Chlorophycean_Mitochondrial_Code = ['FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Trematode_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG']
            Scenedesmus_obliquus_mitochondrial_Code = ['FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Thraustochytrium_Mitochondrial_Code = ['FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
            Pterobranchia_mitochondrial_code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG']
            Candidate_Division_SR1_and_Gracilibacteria_Code = ['FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']

            listBasesAll = Base1 + Base2 + Base3

            ##### Nejprve pro zadane vlakno ("P"):

            peptidP_c00 = []
            peptidP_c01 = []
            peptidP_c02 = []
            peptidP_c03 = []
            peptidP_c04 = []
            peptidP_c05 = []
            peptidP_c06 = []
            peptidP_c07 = []
            peptidP_c08 = []
            peptidP_c09 = []
            peptidP_c10 = []
            peptidP_c11 = []
            peptidP_c12 = []
            peptidP_c13 = []
            peptidP_c14 = []
            peptidP_c15 = []
            peptidP_c16 = []
            peptidP_c17 = []
            peptidP_c18 = []

            #peptidP = []
            #PeptidyP = []

            if (len(listORFsP) == 0):
                print('Zadny ORF a peptid v zadanem vlaknu.')
                print('\n')
            else:


                for iorfxP in range(len(listORFsP)):
                 for kodonP in range(int(len(listORFsP[iorfxP][:])/3)):

                # print('kodon je: ',kodon+1)
                  res1P = [i1P for i1P in range(len(listBasesAll[0])) if listBasesAll[0].startswith(listORFsP[iorfxP][kodonP+0+2*kodonP], i1P)]
                  res2P = [i2P for i2P in range(len(listBasesAll[1])) if listBasesAll[1].startswith(listORFsP[iorfxP][kodonP+1+2*kodonP], i2P)]
                  res3P = [i3P for i3P in range(len(listBasesAll[2])) if listBasesAll[2].startswith(listORFsP[iorfxP][kodonP+2+2*kodonP], i3P)]

                  if ((len(res1P)==0)|(len(res2P)==0)|(len(res3P)==0)):
                      import random
                      AMKP_c00 = random.choice(Standard_Code_nonSTOP[0])
                      #AMKP_c00 = Standard_Code[0][7]
                  else:
                      set_res1P = set(res1P)
                      set_res2P = set(res2P)
                      set_res3P = set(res3P)

                      #print('set_res1P je: ',set_res1P)
                      #print('set_res2P je: ',set_res2P)
                      #print('set_res3P je: ',set_res3P)

                      prunik12P = set_res1P.intersection(set_res2P)
                      prunik123P = prunik12P.intersection(set_res3P)

                      # print('prunik je: ',prunik123P)

                      indexAAsP = int(list(prunik123P)[0])

                      AMKP_c00 = Standard_Code[0][indexAAsP]

                      #print('AMK_c00 je: ',AMKP_c00)

                  peptidP_c00 += [AMKP_c00]

                  sPepP_c00 = ''
                  seqPepP_c00 = peptidP_c00
                  JseqPepP_c00 = sPepP_c00.join(seqPepP_c00)
                  strPepP_c00 = JseqPepP_c00

                  countP_STOP_c00 = strPepP_c00.count('*')
                  PepsP_c00 = strPepP_c00.split('*',countP_STOP_c00)

                  #PepsP = strPepP.split('*',len(listORFsP))
                  #PepsP = PepsP[0:-1]
                  
                  #if AMKP:
                  #print('peptidy jsou: ', PepsP)
                  #PeptidyP += [PepsP]

                  #preramec = []
                  
                  OpeptidesP_c00 = open("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa","w+")
                  for iOpP_c00 in range(len(PepsP_c00)-1):
                        preramecP_c00 = listORFsP_pozice[iorfxP][0]
                        if ((preramecP_c00)%3 == 0):
                            frameP_c00 = 1
                        elif ((preramecP_c00 + 1)%3 == 0):
                            frameP_c00 = 3
                        elif ((preramecP_c00 + 2)%3 == 0):
                            frameP_c00 = 2
                        else:
                            print('Nenasel jsem cteci ramec...')
                        #strORF = str(ORF)
                        #OpeptidesP_c00.write(">" + "frame: " + str(frameP_c00) + ". ," + "ORF-position: " + str(listORFsP_pozice[iOpP_c00]) + "\n" + PepsP_c00[iOpP_c00] + "\n")
                        OpeptidesP_c00.write(">"+str(listORFsP_pozice[iOpP_c00][0])+" "+str(listORFsP_pozice[iOpP_c00][1]) + "\n" + PepsP_c00[iOpP_c00] + "\n")
                  OpeptidesP_c00.close()
                  
                  '''
                  with open('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.fasta','w+') as fPepPfasta_c00:
                    for Pf_c00 in range(len(PepsP_c00)-1):
                        fPepPfasta_c00.write("> ORF: " + str(listORFsP_pozice[Pf_c00]) + "\n" + PepsP_c00[Pf_c00] + "\n")
                  fPepPfasta_c00.close()
                  '''
                  
                  #
                  AMKP_c08 = Bacterial_Archaeal_and_Plant_Plastid_Code[0][indexAAsP]

                  peptidP_c08 += [AMKP_c08]
                  sPepP_c08 = ''
                  seqPepP_c08 = peptidP_c08
                  JseqPepP_c08 = sPepP_c08.join(seqPepP_c08)
                  strPepP_c08 = JseqPepP_c08

                  countP_STOP_c08 = strPepP_c08.count('*')
                  PepsP_c08 = strPepP_c08.split('*',countP_STOP_c08)
                  
                  OpeptidesP_c08 = open("FASTA_output/peptidesP/peptidesP_c08_Bacterial_Archaeal_and_Plant_Plastid_Code.faa","w+")
                  for iOpP_c08 in range(len(PepsP_c08)-1):
                        preramecP_c08 = listORFsP_pozice[iorfxP][0]
                        if ((preramecP_c08)%3 == 0):
                            frameP_c08 = 1
                        elif ((preramecP_c08 + 1)%3 == 0):
                            frameP_c08 = 3
                        elif ((preramecP_c08 + 2)%3 == 0):
                            frameP_c08 = 2
                        else:
                            print('Nenasel jsem cteci ramec...')

                        #OpeptidesP_c08.write(">" + "frame: " + str(frameP_c08) + ". ," + "ORF-position: " + str(listORFsP_pozice[iOpP_c08]) + "\n" + PepsP_c08[iOpP_c08] + "\n")
                        #OpeptidesP_c08.write(">" + "orf"+str(listORFsP_pozice[iOpP_c08][0]) + "\n" + str(PepsP_c08[iOpP_c08]) + "\n")
                        OpeptidesP_c08.write(">"+str(listORFsP_pozice[iOpP_c08][0])+" "+str(listORFsP_pozice[iOpP_c08][1]) + "\n" + PepsP_c08[iOpP_c08] + "\n")
                  OpeptidesP_c08.close()

                  #
                  AMKP_c09 = Alternative_Yeast_Nuclear_Code[0][indexAAsP]

                  peptidP_c09 += [AMKP_c09]
                  sPepP_c09 = ''
                  seqPepP_c09 = peptidP_c09
                  JseqPepP_c09 = sPepP_c09.join(seqPepP_c09)
                  strPepP_c09 = JseqPepP_c09

                  countP_STOP_c09 = strPepP_c09.count('*')
                  PepsP_c09 = strPepP_c09.split('*',countP_STOP_c09)
                  
                  OpeptidesP_c09 = open("FASTA_output/peptidesP/peptidesP_c09_Alternative_Yeast_Nuclear_Code.faa","w+")
                  for iOpP_c09 in range(len(PepsP_c09)-1):
                        preramecP_c09 = listORFsP_pozice[iorfxP][0]
                        if ((preramecP_c09)%3 == 0):
                            frameP_c09 = 1
                        elif ((preramecP_c09 + 1)%3 == 0):
                            frameP_c09 = 3
                        elif ((preramecP_c09 + 2)%3 == 0):
                            frameP_c09 = 2
                        else:
                            print('Nenasel jsem cteci ramec...')
                            
                        #OpeptidesP_c09.write(">" + "frame: " + str(frameP_c09) + ". ," + "ORF-position: " + str(listORFsP_pozice[iOpP_c09]) + "\n" + PepsP_c09[iOpP_c09] + "\n")
                        #OpeptidesP_c09.write(">" + "orf"+str(listORFsP_pozice[iOpP_c09][0]) + "\n" + str(PepsP_c09[iOpP_c09]) + "\n")
                        OpeptidesP_c09.write(">"+str(listORFsP_pozice[iOpP_c09][0])+" "+str(listORFsP_pozice[iOpP_c09][1]) + "\n" + PepsP_c09[iOpP_c09] + "\n")
                  OpeptidesP_c09.close()


                  # ukladam sekvence peptiduu pro zadane vlakno ("P"):
                  # pri kazdem behu programu se soubor cely prepise...
                  # ukladam jako txt a "fasta"

                  with open('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.txt','w+') as fPepP_c00:
                      for P_c00 in PepsP_c00:
                          fPepP_c00.write(str(P_c00) + '\n')
                      fPepP_c00.close()
                  '''
                  #zaloha prechoziho postupu pro "fasta" format:

                  with open('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.fasta','w+') as fPepPfasta:
                      for Pf in PepsP:
                          fPepPfasta.write(str(Pf) + '\n')
                      fPepPfasta.close()
                  '''
                  OpeptidesP_c00 = open("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.fasta", "w")
                  for iOpP_c00 in range(len(PepsP_c00)):
                      OpeptidesP_c00.write(">" + "\n" +PepsP_c00[iOpP_c00] + "\n")
                  #napoveda pro hlavicku:
                  #mame 2 seznamy:
                  # list_seq = [sequence1, sequence2, sequence3, sequence4]
                  # list_name = [name1, name2, name3, name4]
                  #udelame dle:
                  # OpeptidesP.write(">" + list_name[i] + "\n" +list_seq[i] + "\n")
                  OpeptidesP_c00.close()

                  #print('Seznam peptidu pro zadane vlakno je v TXT souboru /home/oem/Documents/PacesDr/peptidesP.txt')
                  #print('Seznam peptidu pro zadane vlakno je ve FASTA souboru /home/oem/Documents/PacesDr/peptidesP.fasta')


            ### Nyni preklad do peptidu pro komplementarni vlakno (zn. "N"):

            print('\n')

            peptidN_c00 = []
            peptidN_c01 = []
            peptidN_c02 = []
            peptidN_c03 = []
            peptidN_c04 = []
            peptidN_c05 = []
            peptidN_c06 = []
            peptidN_c07 = []
            peptidN_c08 = []
            peptidN_c09 = []
            peptidN_c10 = []
            peptidN_c11 = []
            peptidN_c12 = []
            peptidN_c13 = []
            peptidN_c14 = []
            peptidN_c15 = []
            peptidN_c16 = []
            peptidN_c17 = []
            peptidN_c18 = []

            #peptidN = []
            #PeptidyN = []

            if (len(listORFsN) == 0):
                print('Zadny ORF a peptid v komplementarnim vlaknu.')
                print('\n')
            else:

                for iorfxN in range(len(listORFsN)):
                 for kodonN in range(int(len(listORFsN[iorfxN][:])/3)):
                     
                # print('kodon je: ',kodon+1)
                  res1N = [i1N for i1N in range(len(listBasesAll[0])) if listBasesAll[0].startswith(listORFsN[iorfxN][kodonN+0+2*kodonN], i1N)]
                  res2N = [i2N for i2N in range(len(listBasesAll[1])) if listBasesAll[1].startswith(listORFsN[iorfxN][kodonN+1+2*kodonN], i2N)]
                  res3N = [i3N for i3N in range(len(listBasesAll[2])) if listBasesAll[2].startswith(listORFsN[iorfxN][kodonN+2+2*kodonN], i3N)]

                  if ((len(res1N)==0)|(len(res2N)==0)|(len(res3N)==0)):
                      #AMKN_c00 = Standard_Code[0][7]
                      import random
                      AMKN_c00 = random.choice(Standard_Code_nonSTOP[0])
                  else:
                      set_res1N = set(res1N)
                      set_res2N = set(res2N)
                      set_res3N = set(res3N)

                      #print('set_res1N je: ',set_res1N)
                      #print('set_res2N je: ',set_res2N)
                      #print('set_res3N je: ',set_res3N)

                      prunik12N = set_res1N.intersection(set_res2N)
                      prunik123N = prunik12N.intersection(set_res3N)

                    # print('prunik je: ',prunik123N)

                      indexAAsN = int(list(prunik123N)[0])

                      # pro standardni kod (komplementarni, "N" vlakno):

                    #
                      AMKN_c00 = Standard_Code[0][indexAAsN]
                      # print('AMK je: ',AMKN_c00)

                  peptidN_c00 += [AMKN_c00]

                  sPepN_c00 = ''
                  seqPepN_c00 = peptidN_c00
                  JseqPepN_c00 = sPepN_c00.join(seqPepN_c00)
                  strPepN_c00 = JseqPepN_c00

                  countN_STOP_c00 = strPepN_c00.count('*')
                  PepsN_c00 = strPepN_c00.split('*',countN_STOP_c00)
                  
                  OpeptidesN_c00 = open("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa","w+")
                  for iOpN_c00 in range(len(PepsN_c00)-1):
                        preramecN_c00 = listORFsN_pozice[iorfxN][0]
                        if ((preramecN_c00)%3 == 0):
                            frameN_c00 = 1
                        elif ((preramecN_c00 + 1)%3 == 0):
                            frameN_c00 = 3
                        elif ((preramecN_c00 + 2)%3 == 0):
                            frameN_c00 = 2
                        else:
                            print('Nenasel jsem cteci ramec...')

                        #OpeptidesN_c00.write(">" + "frame: " + str(frameN_c00) + ". ," + "ORF-position: " + str(listORFsN_pozice[iOpN_c00]) + "\n" + PepsN_c00[iOpN_c00] + "\n")
                        OpeptidesN_c00.write(">"+str(listORFsN_pozice[iOpN_c00][0])+" "+str(listORFsN_pozice[iOpN_c00][1]) + "\n" + str(PepsN_c00[iOpN_c00]) + "\n")
                  OpeptidesN_c00.close()

                  '''
                  sPepN_c00 = ''
                  seqPepN_c00 = peptidN_c00
                  JseqPepN_c00 = sPepN_c00.join(seqPepN_c00)
                  strPepN_c00 = JseqPepN_c00

                  countN_STOP_c00 = strPepN_c00.count('*')
                  PepsN_c00 = strPepN_c00.split('*',countN_STOP_c00)
                  '''

                  '''
                  while (len(PepsN_c00) < len(listORFsN)):
                      continue
                  if (len(PepsN_c00) == len(listORFsN)):
                      print('Peptidy v komplementarnim vlaknu jsou: \n', PepsN_c00)

                  '''
                  #
                  AMKN_c08 = Bacterial_Archaeal_and_Plant_Plastid_Code[0][indexAAsN]

                  peptidN_c08 += [AMKN_c08]
                  sPepN_c08 = ''
                  seqPepN_c08 = peptidN_c08
                  JseqPepN_c08 = sPepN_c08.join(seqPepN_c08)
                  strPepN_c08 = JseqPepN_c08

                  countN_STOP_c08 = strPepN_c08.count('*')
                  PepsN_c08 = strPepN_c08.split('*',countN_STOP_c08)
                  
                  OpeptidesN_c08 = open("FASTA_output/peptidesN/peptidesN_c08_Bacterial_Archaeal_and_Plant_Plastid_Code.faa","w+")
                  for iOpN_c08 in range(len(PepsN_c08)-1):
                        preramecN_c08 = listORFsN_pozice[iorfxN][0]
                        if ((preramecN_c08)%3 == 0):
                            frameN_c08 = 1
                        elif ((preramecN_c08 + 1)%3 == 0):
                            frameN_c08 = 3
                        elif ((preramecN_c08 + 2)%3 == 0):
                            frameN_c08 = 2
                        else:
                            print('Nenasel jsem cteci ramec...')
                            
                        #OpeptidesN_c08.write(">" + "frame: " + str(frameN_c08) + ". ," + "ORF-position: " + str(listORFsN_pozice[iOpN_c08]) + "\n" + PepsN_c08[iOpN_c08] + "\n")
                        #OpeptidesN_c08.write(">" + "orf"+str(listORFsN_pozice[iOpN_c08][0]) + "\n" + str(PepsN_c08[iOpN_c08]) + "\n")
                        OpeptidesN_c08.write(">"+str(listORFsN_pozice[iOpN_c08][0])+" "+str(listORFsN_pozice[iOpN_c08][1]) + "\n" + PepsN_c08[iOpN_c08] + "\n")
                  OpeptidesN_c08.close()

                  #
                  AMKN_c09 = Alternative_Yeast_Nuclear_Code[0][indexAAsN]

                  peptidN_c09 += [AMKN_c09]
                  sPepN_c09 = ''
                  seqPepN_c09 = peptidN_c09
                  JseqPepN_c09 = sPepN_c09.join(seqPepN_c09)
                  strPepN_c09 = JseqPepN_c09

                  countN_STOP_c09 = strPepN_c09.count('*')
                  PepsN_c09 = strPepN_c09.split('*',countN_STOP_c09)
                  
                  OpeptidesN_c09 = open("FASTA_output/peptidesN/peptidesN_c09_Alternative_Yeast_Nuclear_Code.faa","w+")
                  for iOpN_c09 in range(len(PepsN_c09)-1):
                        preramecN_c09 = listORFsN_pozice[iorfxN][0]
                        if ((preramecN_c09)%3 == 0):
                            frameN_c09 = 1
                        elif ((preramecN_c09 + 1)%3 == 0):
                            frameN_c09 = 3
                        elif ((preramecN_c09 + 2)%3 == 0):
                            frameN_c09 = 2
                        else:
                            print('Nenasel jsem cteci ramec...')
                            
                        #OpeptidesN_c09.write(">" + "frame: " + str(frameN_c09) + ". ," + "ORF-position: " + str(listORFsN_pozice[iOpN_c09]) + "\n" + PepsN_c09[iOpN_c09] + "\n")
                        #OpeptidesN_c09.write(">" + "orf"+str(listORFsN_pozice[iOpN_c09][0]) + "\n" + str(PepsN_c09[iOpN_c09]) + "\n")
                        OpeptidesN_c09.write(">"+str(listORFsN_pozice[iOpN_c09][0])+" "+str(listORFsN_pozice[iOpN_c09][1]) + "\n" + PepsN_c09[iOpN_c09] + "\n")
                  OpeptidesN_c09.close()

                  #print('\n')

                  '''
                  if (len(PepsN_c00) == 0):
                      print('Zadny ORF a peptid v komplementarnim vlaknu.')
                  #elif (len(PepsN_c00) < len(listORFsN)):
                      #continue
                  elif (len(PepsN_c00) == len(listORFsN)):
                      print('Peptidy v komplementarnim vlaknu jsou: \n', PepsN_c00)
                  else:
                      None
                      #continue
                  '''
                  
                  # ukladam sekvence peptiduu pro komplementarni vlakno ("N"):
                  # pri kazdem behu programu se soubor cely prepise...
                  # ukladam jako txt a "fasta"

                  with open('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.txt','w+') as fPepN_c00:
                      for N_c00 in PepsN_c00:
                          fPepN_c00.write(str(N_c00) + '\n')
                      fPepN_c00.close()
                  '''
                  #zaloha prechoziho postupu pro "fasta" format:

                  with open('/home/oem/Documents/PacesDr/peptidesN.fasta','w+') as fPepNfasta:
                      for Nf in PepsN:
                          fPepNfasta.write(str(Nf) + '\n')
                      fPepNfasta.close()
                  '''
                  OpeptidesN_c00 = open("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.fasta", "w")
                  for iOpN_c00 in range(len(PepsN_c00)):
                      OpeptidesN_c00.write(">" + "\n" +PepsN_c00[iOpN_c00] + "\n")
                  #napoveda pro hlavicku:
                  #mame 2 seznamy:
                  # list_seq = [sequence1, sequence2, sequence3, sequence4]
                  # list_name = [name1, name2, name3, name4]
                  #udelame dle:# OpeptidesP.write(">" + list_name[i] + "\n" +list_seq[i] + "\n")
                  OpeptidesN_c00.close()

                  #print('\n')

                  #print('Seznam peptidu pro komplementarni vlakno je v TXT souboru /home/oem/Documents/PacesDr/peptidesN.txt')
                  #print('Seznam peptidu pro komplementarni vlakno je ve FASTA souboru /home/oem/Documents/PacesDr/peptidesN.fasta')

            count_TM_dvoj_obrazku_P = 0
            count_TM_dvoj_obrazku_N = 0

            PepsP_c00_4finalA = []
            
            if (len(listORFsP) == 0):
                print('Zadny ORF a peptid v zadanem vlaknu.')
                print('\n')
            else:
                for irpP_c00 in range(len(PepsP_c00)):
                  if (len(PepsP_c00[irpP_c00]) <= 0):
                      #print(irpP_c00,'je pod limitem')
                      continue
                  elif (len(PepsP_c00[irpP_c00]) > 0):
                      PepsP_c00_4finalA += [PepsP_c00[irpP_c00]]
                      #print(irpP_c00,'je nad limitem')
                  else:
                      break
            
            PepsN_c00_4finalA = []

            
            if (len(listORFsN) == 0):
                print('Zadny ORF a peptid v komplementarnim vlaknu.')
                print('\n')
            else:
                for irpN_c00 in range(len(PepsN_c00)):
                  if (len(PepsN_c00[irpN_c00]) <= 0):
                      #print(irpN_c00,'je pod limitem')
                      continue
                  elif (len(PepsN_c00[irpN_c00]) > 0):
                      PepsN_c00_4finalA += [PepsN_c00[irpN_c00]]
                      #print(irpN_c00,'je nad limitem')
                  else:
                      break
            
            import os.path

            import tmhmm

            import matplotlib.pyplot as plt
            import matplotlib.text
            import matplotlib.cm as cm
            import numpy as np
            import matplotlib.mlab as mlab
            from matplotlib.artist import Artist

            from matplotlib.lines import Line2D
            from matplotlib.patches import Rectangle
            from matplotlib.text import Text
            from matplotlib.image import AxesImage

            import tkinter

            if (len(listORFsP) == 0):
                print('Zadny ORF a peptid v zadanem vlaknu.')
                listTMsP_pozice = []
                print('\n')
            else:

                listTMsP_pozice = []

                count_M_P = []

                apM_P_coords = []

                for irpP_c00_4fA in range(len(PepsP_c00_4finalA)):

                    annotationP, posteriorP = tmhmm.predict(PepsP_c00_4finalA[irpP_c00_4fA],'3','TM_model/TMHMM20.model')

                    print(annotationP, posteriorP)

                    #print(listORFsP_pozice[irpP_c00_4fA])

                    apP = annotationP, posteriorP

                    #import matplotlib.pyplot as plt
                    #import numpy as np

                    data1P = []
                    data2P = []
                    data3P = []

                    for iapP in range(len((annotationP, posteriorP)[0])):
                        data1P += [apP[1][iapP][0]]
                        data2P += [apP[1][iapP][1]]
                        data3P += [apP[1][iapP][2]]

                        if (((len(data1P) & len(data2P) & len(data3P)) == len((annotationP, posteriorP)[0])) & (str('M') in ((annotationP, posteriorP)[0]))):

                          apP0 = (annotationP, posteriorP)[0]
                          apM_P = [apMP for apMP in range(len(apP0)) if apP0.startswith('M', apMP)]
                          min_apM_P = min(apM_P)
                          max_apM_P = max(apM_P)
                          apM_P_coords += [[min_apM_P,max_apM_P]]

                          count_TM_dvoj_obrazku_P += 1 # pocet ORFu s TM sek. strukturou, pocita oba obr. (tento a "linkovy") pro jeden OFR s TM jako jeden
                          #print('pocet ORFu s TM sekundarni strukturou (zadane vlakno): ', count_TM_dvoj_obrazku_P)

                          count_M_P_i = ((annotationP, posteriorP)[0]).count('M')

                          XaxP = np.arange(0,len(PepsP_c00_4finalA[irpP_c00_4fA]),1.0)

                          #plt.plot(XaxP,data1P,'g--',label='intracellular',data2P,'b--',label='in membrane',data3P,'y--',label='extracellular')
                          plt.plot(XaxP,data1P,'g--', label='intracellular')
                          plt.plot(XaxP,data2P,'b--', label='in membrane')
                          plt.plot(XaxP,data3P,'y--', label='extracellular')
                          plt.xlabel('poradi aminokyseliny (N-C konec)')
                          plt.ylabel('pravdepodobnost umisteni')
                          plt.title('TM predikce' + ' - pozice ORFu: ' + str(listORFsP_pozice[irpP_c00_4fA]) + ' bp.')
                          plt.legend()
                          #plt.show()
                          plt.savefig('figures/given/Fig No. ' + str(listORFsP_pozice[irpP_c00_4fA]) + ' provided_string.png')
                          #plt.savefig('/home/oem/Documents/_html/figures/Fig No. ' + str([irpP_c00_4fA]) + ' provided_string.png')
                          plt.close()#

                          # "linkovy obr.":
                          xsP = np.linspace(1,1,lenbP)
                          plt.figure(figsize=(lenbP/120,lenbP/4800))
                          plt.hlines(y=1,xmin=1,xmax=lenbP,color='blue',linestyle='-',lw=20)
                          plt.hlines(y=1,xmin=(listORFsP_pozice[irpP_c00_4fA][0]),xmax=(listORFsP_pozice[irpP_c00_4fA][1]),color='r',linestyle='-',lw=30)
                          plt.yticks([])
                          plt.xticks([])
                          '''
                          if ((listORFsP_pozice[irpP_c00_4fA][0])%3 == 0):
                              frameP_ORFall_lin = 1
                          elif ((listORFsP_pozice[irpP_c00_4fA][0] + 1)%3 == 0):
                              frameP_ORFall_lin = 3
                          elif ((listORFsP_pozice[irpP_c00_4fA][0] + 2)%3 == 0):
                              frameP_ORFall_lin = 2
                          else:
                              continue
                          '''
                          #plt.annotate(str(listORFsP_pozice[irpP_c00_4fA][0]),(listORFsP_pozice[irpP_c00_4fA][0],1),fontsize=16)
                          #plt.annotate(str(listORFsP_pozice[irpP_c00_4fA][1]),(listORFsP_pozice[irpP_c00_4fA][1],1),fontsize=16)
                          #plt.annotate(str(lenbP),(lenbP,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
                          #plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
                          '''
                          plt.annotate('Fr.:'+str(frameP_ORFall_lin),((listORFsP_pozice[irpP_c00_4fA][0]+listORFsP_pozice[irpP_c00_4fA][1])/2,1),fontsize=16,color='green',verticalalignment='top')
                          '''
                          #plt.savefig('/home/oem/Documents/_html/figures/NA string - fig No. ' + str([irpP_c00_4fA]) + ' provided_string.png')
                          plt.savefig('figures/given/NA string - fig No. ' + str(listORFsP_pozice[irpP_c00_4fA]) + ' provided_string.png')
                          #plt.show()
                          plt.close()#

                          listTMsP_pozice += [listORFsP_pozice[irpP_c00_4fA]]

                          count_M_P += [[listORFsP_pozice[irpP_c00_4fA][0], count_M_P_i]]
                          
                          ## listTMsP_pozice += listORFsP_pozice[irpP_c00_4fA]
                          ## pro souhrn typu [2310, 2460, 3150, 3306, 8037, 8202, 8304, 8460, 2368, 2581, 5236, 5443, 6997, 7456, 10210, 10441, 1847, 2000, 3035, 3260, 4175, 4394, 7958, 8330]
                          ## tj. bez vzeti do dvojic
                #print('ORFy s TM strukturou, zadane vlakno:\n',listTMsP_pozice)

            if (len(listORFsN) == 0):
                print('Zadny ORF a peptid v komplementarnim vlaknu.')
                listTMsN_pozice = []
                print('\n')
            else:

                listTMsN_pozice = []

                count_M_N = []

                apM_N_coords = []

                for irpN_c00_4fA in range(len(PepsN_c00_4finalA)):

                    annotationN, posteriorN = tmhmm.predict(PepsN_c00_4finalA[irpN_c00_4fA],'4','TM_model/TMHMM20.model')

                    print(annotationN, posteriorN)

                    apN = annotationN, posteriorN

                    #import matplotlib.pyplot as plt
                    #import numpy as np

                    data1N = []
                    data2N = []
                    data3N = []

                    for iapN in range(len((annotationN, posteriorN)[0])):
                        data1N += [apN[1][iapN][0]]
                        data2N += [apN[1][iapN][1]]
                        data3N += [apN[1][iapN][2]]

                        if (((len(data1N) & len(data2N) & len(data3N)) == len((annotationN, posteriorN)[0])) & (str('M') in ((annotationN, posteriorN)[0]))):

                          apN0 = (annotationN, posteriorN)[0]
                          apM_N = [apMN for apMN in range(len(apN0)) if apN0.startswith('M', apMN)]
                          min_apM_N = min(apM_N)
                          max_apM_N = max(apM_N)
                          apM_N_coords += [[min_apM_N,max_apM_N]]

                          count_TM_dvoj_obrazku_N += 1 # pocet ORFu s TM sek. strukturou, pocita oba obr. pro jeden OFR s TM jako jeden
                          #print('pocet ORFu s TM sekundarni strukturou (komplementarni vlakno): ', count_TM_dvoj_obrazku_N)

                          count_M_N_i = ((annotationN, posteriorN)[0]).count('M')

                          XaxN = np.arange(0,len(PepsN_c00_4finalA[irpN_c00_4fA]),1.0)

                          #plt.plot(XaxN,data1N,'g--',label='intracellular',data2N,'b--',label='in membrane',data3N,'y--',label='extracellular')
                          plt.plot(XaxN,data1N,'g--', label='intracellular')
                          plt.plot(XaxN,data2N,'b--', label='in membrane')
                          plt.plot(XaxN,data3N,'y--', label='extracellular')
                          plt.xlabel('poradi aminokyseliny (N-C konec)')
                          plt.ylabel('pravdepodobnost umisteni')
                          plt.title('TM predikce' + ' - pozice ORFu: ' + str(listORFsN_pozice[irpN_c00_4fA]) + ' bp.')
                          plt.legend()
                          #plt.show()
                          #plt.savefig('/home/oem/Documents/_html/figures/Fig No. ' + str([irpN_c00_4fA]) + ' complementary_string.png')
                          plt.savefig('figures/complementary/Fig No. ' + str(listORFsN_pozice[irpN_c00_4fA]) + ' complementary_string.png')
                          plt.close()#

                          # "linkovy obr.":
                          xsN = np.linspace(1,1,lenbN)
                          plt.figure(figsize=(lenbN/120,lenbN/4800))
                          plt.hlines(y=1,xmin=1,xmax=lenbN,color='blue',linestyle='-',lw=20)
                          plt.hlines(y=1,xmin=(listORFsN_pozice[irpN_c00_4fA][0]),xmax=(listORFsN_pozice[irpN_c00_4fA][1]),color='r',linestyle='-',lw=30)
                          plt.yticks([])
                          plt.xticks([])
                          '''
                          if ((listORFsN_pozice[irpN_c00_4fA][0])%3 == 0):
                              frameN_ORFall_lin = 1
                          elif ((listORFsN_pozice[irpN_c00_4fA][0] + 1)%3 == 0):
                              frameN_ORFall_lin = 3
                          elif ((listORFsN_pozice[irpN_c00_4fA][0] + 2)%3 == 0):
                              frameN_ORFall_lin = 2
                          else:
                              continue
                          '''
                          #plt.annotate(str(listORFsN_pozice[irpN_c00_4fA][0]),(listORFsN_pozice[irpN_c00_4fA][0],1),fontsize=16)
                          #plt.annotate(str(listORFsN_pozice[irpN_c00_4fA][1]),(listORFsN_pozice[irpN_c00_4fA][1],1),fontsize=16)
                          #plt.annotate(str(lenbN),(lenbN,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
                          #plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
                          '''
                          plt.annotate('Fr.:'+str(frameN_ORFall_lin),((listORFsN_pozice[irpN_c00_4fA][0]+listORFsN_pozice[irpN_c00_4fA][1])/2,1),fontsize=16,color='green',verticalalignment='top')
                          '''
                          plt.savefig('figures/complementary/NA string - fig No. ' + str(listORFsN_pozice[irpN_c00_4fA]) + ' complementary_string.png')
                          #plt.savefig('/home/oem/Documents/_html/figures/NA string - fig No. ' + str([irpN_c00_4fA]) + ' complementary_string.png')
                          #plt.show()
                          plt.close()#

                          listTMsN_pozice += [listORFsN_pozice[irpN_c00_4fA]]

                          count_M_N += [[listORFsN_pozice[irpN_c00_4fA][0], count_M_N_i]]
                          
                          ## listTMsN_pozice += listORFsN_pozice[irpN_c00_4fA]
                          ## pro souhrn typu [2310, 2460, 3150, 3306, 8037, 8202, 8304, 8460, 2368, 2581, 5236, 5443, 6997, 7456, 10210, 10441, 1847, 2000, 3035, 3260, 4175, 4394, 7958, 8330]
                          ## tj. bez vzeti do dvojic
                #print('ORFy s TM strukturou, komplementarni vlakno:\n',listTMsN_pozice)

            # urcovani CC struktury (zadane i komlementarni vlakno) - zacatek:

            #CC_threshold_P = 0.10
            #CC_threshold_N = 0.10

            CC_count_P = 0
            CC_count_N = 0

            CC_ORFstart_P = []
            CC_ORFstart_N = []

            max_results_CC_P = []
            max_results_CC_N = []

            ind_MrkcP_coords = []
            ind_MrkcN_coords = []
            
            if (len(listORFsP) == 0):
                print('Zadny ORF a peptid v zadanem vlaknu. Nelze hledat CC strukturu.')
                print('\n')
            else:
                from deepcoil import DeepCoil
                from deepcoil.utils import plot_preds
                from Bio import SeqIO

                dc = DeepCoil(use_gpu=False)
                inp = {str(entry.id): str(entry.seq) for entry in SeqIO.parse('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa', 'fasta')}
                results = dc.predict(inp)

                for key in results.keys():
                    plot_preds(results[key], out_file='CC_plots/plotsP/{}.png'.format(key))

                for key in results.keys():
                    #print('Max z ['+key+'] CC (zadane vlakno) je: ',max(results[key]['cc']))
                    if (max(results[key]['cc']) >= CC_threshold_P):
                        plot_preds(results[key], out_file='CC_plots/plotsP-realCC/{}.png'.format(key))
                        CC_count_P += 1
                        CC_ORFstart_P += [int(key)]
                        max_results_CC_P += [[int(key), max(results[key]['cc'])]]

                        MrkcP = max(results[key]['cc'])
                        list_rkcP = (results[key]['cc']).tolist()
                        ind_MrkcP = list_rkcP.index(MrkcP) # udava jen 1 cislo - prvni vyskyt z leva, cislovano od 0 (nuly)
                        ###min_ind_MrkcP = min(ind_MrkcP) # nelze - min jde ze seznamu, ne z jednoho cisla 'int' - neni treba - viz koment. k ind_MrkcP o 1 radek vyse
                        #print(ind_MrkcP)
                        ind_MrkcP_coords += [ind_MrkcP]

                    else:
                        None
            
            if (len(listORFsN) == 0):
                print('Zadny ORF a peptid v komplementarnim vlaknu. Nelze hledat CC strukturu.')
                print('\n')
            else:
                from deepcoil import DeepCoil
                from deepcoil.utils import plot_preds
                from Bio import SeqIO

                dc = DeepCoil(use_gpu=False)
                inp = {str(entry.id): str(entry.seq) for entry in SeqIO.parse('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa', 'fasta')}
                results = dc.predict(inp)

                for key in results.keys():
                    plot_preds(results[key], out_file='CC_plots/plotsN/{}.png'.format(key))

                for key in results.keys():
                    #print('Max z ['+key+'] CC (komplementarni vlakno) je: ',max(results[key]['cc']))
                    if (max(results[key]['cc']) >= CC_threshold_N):
                        plot_preds(results[key], out_file='CC_plots/plotsN-realCC/{}.png'.format(key))
                        CC_count_N += 1
                        CC_ORFstart_N += [int(key)]
                        max_results_CC_N += [[int(key), max(results[key]['cc'])]]

                        MrkcN = max(results[key]['cc'])
                        list_rkcN = (results[key]['cc']).tolist()
                        ind_MrkcN = list_rkcN.index(MrkcN) # udava jen 1 cislo - prvni vyskyt z leva, cislovano od 0 (nuly)
                        ###min_ind_MrkcN = min(ind_MrkcN) # nelze - min jde ze seznamu, ne z jednoho cisla 'int' - neni treba - viz koment. k ind_MrkcN o 1 radek vyse
                        #print(ind_MrkcN)
                        ind_MrkcN_coords += [ind_MrkcN]
                        
                    else:
                        None
            
            # urcovani CC struktury (zadane i komlementarni vlakno) - konec.

            # urcovani obrazku/ORFu s jen CC strukturou (ale muze tam byt i GPI struktura), jen zadane, "P" vlakno - zacatek:

            listTMsP_pozice_0 = []

            if len(listTMsP_pozice) == 0 :
                print('Zadny ORF s TM strukturou v zadanem vlaknu.')

            for itmPp in range(len(listTMsP_pozice)):
                listTMsP_pozice_0 += [listTMsP_pozice[itmPp][0]]

            set_listTMsP_pozice_0 = set(listTMsP_pozice_0)
            set_CC_ORFstart_P = set(CC_ORFstart_P)

            CC_notTM_P_set = set_CC_ORFstart_P.difference(set_listTMsP_pozice_0)
            CC_notTM_P_list = list(CC_notTM_P_set)

            CC_notTM_P_coords = []

            for iOPp in listORFsP_pozice:
                for iCP in CC_notTM_P_list:
                    if (iCP == iOPp[0]):
                        CC_notTM_P_coords += [iOPp]
                    else:
                        None

            # cilem zde je hlavne: CC_notTM_P_coords
            
            TM_notCC_P_set = set_listTMsP_pozice_0.difference(set_CC_ORFstart_P)    

            TM_a_CC_P_set = set_listTMsP_pozice_0.intersection(set_CC_ORFstart_P)    

            # urcovani obrazku/ORFu s jen CC strukturou (ale muze tam byt i GPI struktura), jen zadane, "P" vlakno - konec.

            # urcovani ORFu s CC strukturou (neresi se zda je pritomna dalsi sek. struktura), zadane vlakno - zacatek:

            CC_all_P_coords = []

            for iOPcc in listORFsP_pozice:
                for iPcc in CC_ORFstart_P:
                    if (iPcc == iOPcc[0]):
                        CC_all_P_coords += [iOPcc]
                    else:
                        None

            # urcovani ORFu s CC strukturou (neresi se zda je pritomna dalsi sek. struktura), zadane vlakno - konec.

            # urcovani obrazku/ORFu s jen CC strukturou (ale muze tam byt i GPI struktura), jen komplementarni, "N" vlakno - zacatek:

            listTMsN_pozice_0 = []

            for itmNp in range(len(listTMsN_pozice)):
                listTMsN_pozice_0 += [listTMsN_pozice[itmNp][0]]

            set_listTMsN_pozice_0 = set(listTMsN_pozice_0)
            set_CC_ORFstart_N = set(CC_ORFstart_N)

            CC_notTM_N_set = set_CC_ORFstart_N.difference(set_listTMsN_pozice_0)
            CC_notTM_N_list = list(CC_notTM_N_set)

            CC_notTM_N_coords = []

            for iONp in listORFsN_pozice:
                for iCN in CC_notTM_N_list:
                    if (iCN == iONp[0]):
                        CC_notTM_N_coords += [iONp]
                    else:
                        None

            # cilem zde je hlavne: CC_notTM_N_coords
            
            TM_notCC_N_set = set_listTMsN_pozice_0.difference(set_CC_ORFstart_N)    

            TM_a_CC_N_set = set_listTMsN_pozice_0.intersection(set_CC_ORFstart_N)    

            # urcovani obrazku/ORFu s jen CC strukturou (ale muze tam byt i GPI struktura), jen komplementarni, "N" vlakno - konec.

            # urcovani ORFu s CC strukturou (neresi se zda je pritomna dalsi sek. struktura), komplementarni vlakno - zacatek:

            CC_all_N_coords = []

            for iONcc in listORFsN_pozice:
                for iNcc in CC_ORFstart_N:
                    if (iNcc == iONcc[0]):
                        CC_all_N_coords += [iONcc]
                    else:
                        None

            # urcovani ORFu s CC strukturou (neresi se zda je pritomna dalsi sek. struktura), komplementarni vlakno - konec.    

            # tvorba "linkoveho" obr. pro ORF (mozna) jen s CC (je mozna GPI) strukturou, zadane, "P" vlakno - zacatek:

            for iLccP in range(len(CC_notTM_P_coords)):

                xsP = np.linspace(1,1,lenbP)
                plt.figure(figsize=(lenbP/120,lenbP/4800))
                plt.hlines(y=1,xmin=1,xmax=lenbP,color='blue',linestyle='-',lw=20)
                plt.hlines(y=1,xmin=(CC_notTM_P_coords[iLccP][0]),xmax=(CC_notTM_P_coords[iLccP][1]),color='r',linestyle='-',lw=30)
                plt.yticks([])
                plt.xticks([])
                if ((CC_notTM_P_coords[iLccP][0])%3 == 0):
                  frameP_iLccP = 1
                elif ((CC_notTM_P_coords[iLccP][0] + 1)%3 == 0):
                  frameP_iLccP = 3
                elif ((CC_notTM_P_coords[iLccP][0] + 2)%3 == 0):
                  frameP_iLccP = 2
                else:
                  continue
                plt.annotate(str(CC_notTM_P_coords[iLccP][0]),(CC_notTM_P_coords[iLccP][0],1),fontsize=16)
                plt.annotate(str(CC_notTM_P_coords[iLccP][1]),(CC_notTM_P_coords[iLccP][1],1),fontsize=16)
                plt.annotate(str(lenbP),(lenbP,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
                plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
                #
                plt.annotate('Fr.:'+str(frameP_iLccP),((CC_notTM_P_coords[iLccP][0]+CC_notTM_P_coords[iLccP][1])/2,1),fontsize=16,color='green',verticalalignment='top')
                plt.savefig('CC_plots/plotsP-realCC/line-picture-CC-P/NA string - fig No. ' + str([iLccP]) + ' provided_string.png')
                #plt.show()
                plt.close()#

            # tvorba "linkoveho" obr. pro ORF (mozna) jen s CC (je mozna GPI) strukturou, zadane, "P" vlakno - konec.

            # tvorba "linkoveho" obr. pro ORF (mozna) jen s CC (je mozna GPI) strukturou, komplementarni, "N" vlakno - zacatek:

            for iLccN in range(len(CC_notTM_N_coords)):

                xsN = np.linspace(1,1,lenbN)
                plt.figure(figsize=(lenbN/120,lenbN/4800))
                plt.hlines(y=1,xmin=1,xmax=lenbN,color='blue',linestyle='-',lw=20)
                plt.hlines(y=1,xmin=(CC_notTM_N_coords[iLccN][0]),xmax=(CC_notTM_N_coords[iLccN][1]),color='r',linestyle='-',lw=30)
                plt.yticks([])
                plt.xticks([])
                if ((CC_notTM_N_coords[iLccN][0])%3 == 0):
                  frameN_iLccN = 1
                elif ((CC_notTM_N_coords[iLccN][0] + 1)%3 == 0):
                  frameN_iLccN = 3
                elif ((CC_notTM_N_coords[iLccN][0] + 2)%3 == 0):
                  frameN_iLccN = 2
                else:
                  continue
                plt.annotate(str(CC_notTM_N_coords[iLccN][0]),(CC_notTM_N_coords[iLccN][0],1),fontsize=16)
                plt.annotate(str(CC_notTM_N_coords[iLccN][1]),(CC_notTM_N_coords[iLccN][1],1),fontsize=16)
                plt.annotate(str(lenbN),(lenbN,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
                plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
                #
                plt.annotate('Fr.:'+str(frameN_iLccN),((CC_notTM_N_coords[iLccN][0]+CC_notTM_N_coords[iLccN][1])/2,1),fontsize=16,color='green',verticalalignment='top')
                plt.savefig('CC_plots/plotsN-realCC/line-picture-CC-N/NA string - fig No. ' + str([iLccN]) + ' complementary_string.png')
                #plt.show()
                plt.close()#

            # tvorba "linkoveho" obr. pro ORF (mozna) jen s CC (je mozna GPI) strukturou, komplementarni, "N" vlakno - konec.

            # 2D obrazky pro souhrn ORFu, pro zadane vlakno, po jednotlivych ctecich ramcich - zacatek :

            from urllib import *
            import json
            from matplotlib import *
            from tkinter import *
            import matplotlib.pyplot as plt
            from matplotlib.lines import Line2D
            from matplotlib.patches import Rectangle
            from matplotlib.text import Text
            import matplotlib.image
            import numpy as np 

            import webbrowser

            ###

            if (len(listORFsP_1) == 0):
                print('Zadny ORF a peptid v zadanem vlaknu, v 1. ramci.')
                print('\n')
            else:

                xP1 = []
                yP1 = []

                for iAOP1 in range(len(listORFsP_1_pozice)):

                    xP1 += [listORFsP_1_pozice[iAOP1][1]]
                    yP1 += [listORFsP_1_pozice[iAOP1][0]]

                liste_xP1 = xP1
                liste_yP1 = yP1

                liste_stringP1 = liste_yP1

                fig = plt.figure()

                fig, ax = plt.subplots() #

                ax.set_xlim([0,lenbP])
                ax.set_ylim([0,lenbP])


                plt.figure(figsize=(lenbP/120,lenbP/120))
                plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
                plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
                plt.title('Souhrn ORFu - zadane vlakno \n (1. cteci ramec)', fontsize=round(lenbP/65))

                #ax = plt.gca()
                #fig = plt.gcf()
                    
                for sxP1, dyP1 in zip(liste_xP1, liste_yP1):
                    plt.annotate(' '+str(dyP1), xy = (sxP1, dyP1), fontsize=round(lenbP/120))
                    plt.annotate('            -  '+str(sxP1), xy = (sxP1, dyP1), fontsize=round(lenbP/120))
                    plt.scatter(liste_xP1, liste_stringP1, s = round(lenbP/5), c ='r') # picker=(True & 10000)

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/NA-ORFs-Pframe1-given_string-2D.svg')
                plt.savefig('figures/NA-ORFs-Pframe1-given_string-2D.png', dpi=32) # , dpi=16
                plt.close()

            ###

            if (len(listORFsP_2) == 0):
                print('Zadny ORF a peptid v zadanem vlaknu, v 2. ramci.')
                print('\n')
            else:

                xP2 = []
                yP2 = []

                for iAOP2 in range(len(listORFsP_2_pozice)):

                    xP2 += [listORFsP_2_pozice[iAOP2][1]]
                    yP2 += [listORFsP_2_pozice[iAOP2][0]]

                liste_xP2 = xP2
                liste_yP2 = yP2

                liste_stringP2 = liste_yP2

                fig = plt.figure()

                fig, ax = plt.subplots() ###

                plt.figure(figsize=(lenbP/120,lenbP/120))
                plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
                plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
                plt.title('Souhrn ORFu - zadane vlakno \n (2. cteci ramec)', fontsize=round(lenbP/65))

                ax = plt.gca()
                fig = plt.gcf()

                for sxP2, dyP2 in zip(liste_xP2, liste_yP2):
                    plt.annotate(' '+str(dyP2), xy = (sxP2, dyP2), fontsize=round(lenbP/120))
                    plt.annotate('            -  '+str(sxP2), xy = (sxP2, dyP2), fontsize=round(lenbP/120))
                    plt.scatter(liste_xP2, liste_stringP2, s = round(lenbP/5), c ='b') # picker=(True & 10000)

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/NA-ORFs-Pframe2-given_string-2D.svg')
                plt.savefig('figures/NA-ORFs-Pframe2-given_string-2D.png', dpi=32) # , dpi=16
                plt.close()

            ###

            if (len(listORFsP_3) == 0):
                print('Zadny ORF a peptid v zadanem vlaknu, v 3. ramci.')
                print('\n')
            else:

                xP3 = []
                yP3 = []

                for iAOP3 in range(len(listORFsP_3_pozice)):

                    xP3 += [listORFsP_3_pozice[iAOP3][1]]
                    yP3 += [listORFsP_3_pozice[iAOP3][0]]

                liste_xP3 = xP3
                liste_yP3 = yP3

                liste_stringP3 = liste_yP3

                fig = plt.figure()

                fig, ax = plt.subplots() ###

                plt.figure(figsize=(lenbP/120,lenbP/120))
                plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
                plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
                plt.title('Souhrn ORFu - zadane vlakno \n (3. cteci ramec)', fontsize=round(lenbP/65))

                ax = plt.gca()
                fig = plt.gcf()

                for sxP3, dyP3 in zip(liste_xP3, liste_yP3):
                    plt.annotate(' '+str(dyP3), xy = (sxP3, dyP3), fontsize=round(lenbP/120))
                    plt.annotate('            -  '+str(sxP3), xy = (sxP3, dyP3), fontsize=round(lenbP/120))
                    plt.scatter(liste_xP3, liste_stringP3, s = round(lenbP/5), c ='g') # picker=(True & 10000)


                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                ##plt.scatter(liste_x, liste_string, s = 10000, c ='r', alpha = 0.8).set_urls(['/home/oem/Documents/_html/index_123.html/#:~:text=['+str(s)+', '+str(d)+']'])
                #plt.plot(liste_x, liste_string,'x', markersize=128, color='black')

                #plt.show()
                #plt.savefig('figures/NA-ORFs-Pframe3-given_string-2D.svg')
                plt.savefig('figures/NA-ORFs-Pframe3-given_string-2D.png', dpi=32) # , dpi=16
                plt.close()

            # 2D obrazky pro souhrn ORFu, pro zadane vlakno, po jednotlivych ctecich ramcich - konec.

            # 2D obrazky pro souhrn ORFu, pro komplementarni vlakno, po jednotlivych ctecich ramcich - zacatek:
            # Tj. nyni obrazky s ORFy pro komplementarni vlakno - zacatek:

            if (len(listORFsN_1) == 0):
                print('Zadny ORF a peptid v komplementarnim vlaknu, v 1. ramci.')
                print('\n')
            else:

                xN1 = []
                yN1 = []

                for iAON1 in range(len(listORFsN_1_pozice)):

                    xN1 += [listORFsN_1_pozice[iAON1][1]]
                    yN1 += [listORFsN_1_pozice[iAON1][0]]

                liste_xN1 = xN1
                liste_yN1 = yN1

                liste_stringN1 = liste_yN1

                fig = plt.figure()

                fig, ax = plt.subplots() ###

                plt.figure(figsize=(lenbN/120,lenbN/120))
                plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
                plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
                plt.title('Souhrn ORFu - komplementarni vlakno \n (1. cteci ramec)', fontsize=round(lenbP/65))

                ax = plt.gca()
                fig = plt.gcf()
                    
                for sxN1, dyN1 in zip(liste_xN1, liste_yN1):
                    plt.annotate(' '+str(dyN1), xy = (sxN1, dyN1), fontsize=round(lenbP/120))
                    plt.annotate('            -  '+str(sxN1), xy = (sxN1, dyN1), fontsize=round(lenbP/120))
                    plt.scatter(liste_xN1, liste_stringN1, s = round(lenbP/5), c ='r') # picker=(True & 10000)

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/NA-ORFs-Nframe1-complementary_string-2D.svg')
                plt.savefig('figures/NA-ORFs-Nframe1-complementary_string-2D.png', dpi=32) # , dpi=16
                plt.close()

            ###

            if (len(listORFsN_2) == 0):
                print('Zadny ORF a peptid v komplementarnim vlaknu, v 2. ramci.')
                print('\n')
            else:

                xN2 = []
                yN2 = []

                for iAON2 in range(len(listORFsN_2_pozice)):

                    xN2 += [listORFsN_2_pozice[iAON2][1]]
                    yN2 += [listORFsN_2_pozice[iAON2][0]]

                liste_xN2 = xN2
                liste_yN2 = yN2

                liste_stringN2 = liste_yN2

                fig = plt.figure()

                fig, ax = plt.subplots() ###

                plt.figure(figsize=(lenbN/120,lenbN/120))
                plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
                plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
                plt.title('Souhrn ORFu - komplementarni vlakno \n (2. cteci ramec)', fontsize=round(lenbP/65))

                ax = plt.gca()
                fig = plt.gcf()

                for sxN2, dyN2 in zip(liste_xN2, liste_yN2):
                    plt.annotate(' '+str(dyN2), xy = (sxN2, dyN2), fontsize=round(lenbP/120))
                    plt.annotate('            -  '+str(sxN2), xy = (sxN2, dyN2), fontsize=round(lenbP/120))
                    plt.scatter(liste_xN2, liste_stringN2, s = round(lenbP/5), c ='b') # picker=(True & 10000)

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/NA-ORFs-Nframe2-complementary_string-2D.svg')
                plt.savefig('figures/NA-ORFs-Nframe2-complementary_string-2D.png', dpi=32) # , dpi=16
                plt.close()

            ###

            if (len(listORFsN_3) == 0):
                print('Zadny ORF a peptid v komplementarnim vlaknu, v 3. ramci.')
                print('\n')
            else:

                xN3 = []
                yN3 = []

                for iAON3 in range(len(listORFsN_3_pozice)):

                    xN3 += [listORFsN_3_pozice[iAON3][1]]
                    yN3 += [listORFsN_3_pozice[iAON3][0]]

                liste_xN3 = xN3
                liste_yN3 = yN3

                liste_stringN3 = liste_yN3

                fig = plt.figure()

                fig, ax = plt.subplots() ###

                plt.figure(figsize=(lenbN/120,lenbN/120))
                plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
                plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
                plt.title('Souhrn ORFu - komplementarni vlakno \n (3. cteci ramec)', fontsize=round(lenbP/65))

                ax = plt.gca()
                fig = plt.gcf()

                for sxN3, dyN3 in zip(liste_xN3, liste_yN3):
                    plt.annotate(' '+str(dyN3), xy = (sxN3, dyN3), fontsize=round(lenbP/120))
                    plt.annotate('            -  '+str(sxN3), xy = (sxN3, dyN3), fontsize=round(lenbP/120))
                    plt.scatter(liste_xN3, liste_stringN3, s = round(lenbP/5), c ='g') # picker=(True & 10000)


                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                ##plt.scatter(liste_x, liste_string, s = 10000, c ='r', alpha = 0.8).set_urls(['/home/oem/Documents/_html/index_123.html/#:~:text=['+str(s)+', '+str(d)+']'])
                #plt.plot(liste_x, liste_string,'x', markersize=128, color='black')

                #plt.show()
                #plt.savefig('figures/NA-ORFs-Nframe3-complementary_string-2D.svg')
                plt.savefig('figures/NA-ORFs-Nframe3-complementary_string-2D.png', dpi=32) # , dpi=16
                plt.close()

            # 2D obrazky pro souhrn ORFu, pro komplementarni vlakno, po jednotlivych ctecich ramcich - konec.
            # Tj. nyni pro obrazky s OFRy pro komplementarni vlako - konec.

            '''
            im = plt.imread('/home/oem/Documents/_html/figures/NA-ORFs-Pframe1-given_string-2D.svg')
            implot = ax.imshow(im)

            def onclick(event):
                for sx, dy in zip(liste_x, liste_y):
                    if event.xdata == sx and event.ydata == dy:
                        webbrowser.open('https://www.bbc.com')
                    else:
                        continue
            cid = fig.canvas.mpl_connect('press_event', onclick)
            '''
            #plt.savefig('/home/oem/Documents/_html/figures/NA-ORFs-all-given_string-2D.png')
            
            # puvodni 2D obr. pro souhrn vsech ORFu, nebo jednoho/1. ramce, zadane vlakno - konec.

            # urocvani  ORFu (pozice) s jen TM, s TM a CC, a jen s CC, sek. strukturou (mnozinove operace), zadane vlakno - zacatek:
            # dulezite vstupni promenne:
            # listTMsP_pozice
            # CC_all_P_coords
            # listORFsP_pozice

            listTMsP_pozice2_0 = []
            CC_all_P_coords_0 = []

            for itmPp2 in range(len(listTMsP_pozice)):
                listTMsP_pozice2_0 += [listTMsP_pozice[itmPp2][0]]

            set_listTMsP_pozice_0 = set(listTMsP_pozice2_0)

            for iccPp2 in range(len(CC_all_P_coords)):
                CC_all_P_coords_0 += [CC_all_P_coords[iccPp2][0]]    
            
            set_CC_all_P_coords_0 = set(CC_all_P_coords_0)
                
            set_only_TMsP_pozice_0 = set_listTMsP_pozice_0.difference(set_CC_all_P_coords_0)

            set_CC_only_P_coords_0 = set_CC_all_P_coords_0.difference(set_listTMsP_pozice_0)

            set_intersection_TM_CC_P_0 = set_listTMsP_pozice_0.intersection(set_CC_all_P_coords_0)

            intersection_TM_CC_P_0 = list(set_intersection_TM_CC_P_0)

            only_TMsP_pozice_0 = list(set_only_TMsP_pozice_0)

            CC_only_P_coords_0 = list(set_CC_only_P_coords_0)

            intersection_TM_CC_P = []

            for iOPx1 in listORFsP_pozice:
                for itmccP0 in intersection_TM_CC_P_0:
                    if (itmccP0 == iOPx1[0]):
                        intersection_TM_CC_P += [iOPx1]
                    else:
                        None

            only_TMsP_pozice = []

            for iOPx2 in listORFsP_pozice:
                for itmP0 in only_TMsP_pozice_0:
                    if (itmP0 == iOPx2[0]):
                        only_TMsP_pozice += [iOPx2]
                    else:
                        None

            CC_only_P_coords = []

            for iOPx3 in listORFsP_pozice:
                for iccP0 in CC_only_P_coords_0:
                    if (iccP0 == iOPx3[0]):
                        CC_only_P_coords += [iOPx3]
                    else:
                        None

            # ziskane promenne , hlavne:
            # intersection_TM_CC_P
            # only_TMsP_pozice
            # CC_only_P_coords

            # urocvani ORFu (pozice) s jen TM, s TM a CC, a jen s CC, sek. strukturou (mnozinove operace), zadane vlakno - konec.

            # urocvani ORFu (pozice) s jen TM, s TM a CC, a jen s CC, sek. strukturou (mnozinove operace), komplementarni vlakno - zacatek:
            # dulezite vstupni promenne:
            # listTMsN_pozice
            # CC_all_N_coords
            # listORFsN_pozice

            listTMsN_pozice2_0 = []
            CC_all_N_coords_0 = []
            
            for itmNp2 in range(len(listTMsN_pozice)):
                listTMsN_pozice2_0 += [listTMsN_pozice[itmNp2][0]]

            set_listTMsN_pozice_0 = set(listTMsN_pozice2_0)

            for iccNp2 in range(len(CC_all_N_coords)):
                CC_all_N_coords_0 += [CC_all_N_coords[iccNp2][0]]    
            
            set_CC_all_N_coords_0 = set(CC_all_N_coords_0)
                
            set_only_TMsN_pozice_0 = set_listTMsN_pozice_0.difference(set_CC_all_N_coords_0)

            set_CC_only_N_coords_0 = set_CC_all_N_coords_0.difference(set_listTMsN_pozice_0)

            set_intersection_TM_CC_N_0 = set_listTMsN_pozice_0.intersection(set_CC_all_N_coords_0)

            intersection_TM_CC_N_0 = list(set_intersection_TM_CC_N_0)    

            only_TMsN_pozice_0 = list(set_only_TMsN_pozice_0)

            CC_only_N_coords_0 = list(set_CC_only_N_coords_0)

            intersection_TM_CC_N = []

            for iONx1 in listORFsN_pozice:
                for itmccN0 in intersection_TM_CC_N_0:
                    if (itmccN0 == iONx1[0]):
                        intersection_TM_CC_N += [iONx1]
                    else:
                        None

            only_TMsN_pozice = []

            for iONx2 in listORFsN_pozice:
                for itmN0 in only_TMsN_pozice_0:
                    if (itmN0 == iONx2[0]):
                        only_TMsN_pozice += [iONx2]
                    else:
                        None

            CC_only_N_coords = []

            for iONx3 in listORFsN_pozice:
                for iccN0 in CC_only_N_coords_0:
                    if (iccN0 == iONx3[0]):
                        CC_only_N_coords += [iONx3]
                    else:
                        None

            # ziskane promenne , hlavne:
            # intersection_TM_CC_N
            # only_TMsN_pozice
            # CC_only_N_coords
            
            # urcovani ORFu (pozice) s jen TM, s TM a CC, a jen s CC, sek. strukturou (mnozinove operace), komplementarni vlakno - konec.

            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - zacatek:
            # dava souradnice ORFu a CC struktury
            
            if (len(CC_all_P_coords) == 0):
                print('Zadny ORF s CC sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:
                xPccp = []
                yPccp = []

                for iPccp in range(len(CC_all_P_coords)):

                    xPccp += [max_results_CC_P[iPccp][0] + 3*ind_MrkcP_coords[iPccp] + 3]
                    yPccp += [round(max_results_CC_P[iPccp][1]*1000)/10]

                liste_xPccp = xPccp
                liste_yPccp = yPccp

                liste_xPccp = [0] + liste_xPccp
                liste_xPccp = liste_xPccp + [lenbP]

                liste_yPccp = [0] + liste_yPccp
                liste_yPccp = liste_yPccp + [0]

                liste_stringPccp = liste_yPccp

                fig = plt.figure()
                
                fig, ax = plt.subplots() ###

                plt.figure(figsize=(round(lenbP/120),round(lenbP/75)))
                plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, zadane vlakno', fontsize=round(lenbP/65))
                plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbP/65))
                #plt.title('Souhrn ORFu s CC strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxPccp, dyPccp in zip(liste_xPccp, liste_yPccp):
                    if (sxPccp == 0):
                        plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                        plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                    elif (sxPccp > 0) and (sxPccp < lenbP):
                        for sxPccpX, dyPccpX in zip(liste_xPccp[1:-1], liste_yPccp[1:-1]):
                            plt.annotate(''+str(dyPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifPccX = liste_xPccp[1:-1].index(sxPccpX)
                            plt.annotate('     - ORF:'+str(sxPccpX - 3*ind_MrkcP_coords[ifPccX] - 3)+' - CC:'+str(sxPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxPccp == lenbP):
                        plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                        plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                    else:
                        None
                    #plt.scatter(liste_xPccp, liste_stringPccp, s = lenbP/1.5, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xPccp, liste_stringPccp, width = 10, align='center', color='black') # width = lenbP/120
                        
                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno.svg')
                plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno.png', dpi=32) # , dpi=16
                plt.close()

            # dava souradnice ORFu a CC struktury
            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - zacatek:
            # dava souradnice ORFu a CC struktury
            # pro BIG verzi - roztazeni obr. pres celou str.
            
            if (len(CC_all_P_coords) == 0):
                print('Zadny ORF s CC sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:
                xPccp = []
                yPccp = []

                for iPccp in range(len(CC_all_P_coords)):

                    xPccp += [max_results_CC_P[iPccp][0] + 3*ind_MrkcP_coords[iPccp] + 3]
                    yPccp += [round(max_results_CC_P[iPccp][1]*1000)/10]

                liste_xPccp = xPccp
                liste_yPccp = yPccp

                liste_xPccp = [0] + liste_xPccp
                liste_xPccp = liste_xPccp + [lenbP]

                liste_yPccp = [0] + liste_yPccp
                liste_yPccp = liste_yPccp + [0]

                liste_stringPccp = liste_yPccp

                fig = plt.figure()
                
                fig, ax = plt.subplots() ###

                plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))
                plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, zadane vlakno', fontsize=round(lenbP/65))
                plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbP/65))
                #plt.title('Souhrn ORFu s CC strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxPccp, dyPccp in zip(liste_xPccp, liste_yPccp):
                    if (sxPccp == 0):
                        plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                        plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                    elif (sxPccp > 0) and (sxPccp < lenbP):
                        for sxPccpX, dyPccpX in zip(liste_xPccp[1:-1], liste_yPccp[1:-1]):
                            plt.annotate(''+str(dyPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifPccX = liste_xPccp[1:-1].index(sxPccpX)
                            plt.annotate('     - ORF:'+str(sxPccpX - 3*ind_MrkcP_coords[ifPccX] - 3)+' - CC:'+str(sxPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxPccp == lenbP):
                        plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                        plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                    else:
                        None
                    #plt.scatter(liste_xPccp, liste_stringPccp, s = lenbP/1.5, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xPccp, liste_stringPccp, width = 10, align='center', color='black') # width = lenbP/120
                        
                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG.svg')
                plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG.png', dpi=32) # , dpi=16
                plt.close()

            # pro BIG verzi - roztazeni obr. pres celou str.
            # dava souradnice ORFu a CC struktury
            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - zacatek:
            # dava souradnice ORFu a CC struktury
            # pro BIG verzi - roztazeni obr. pres celou str.
            # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro CC
            
            if (len(CC_all_P_coords) == 0):
                print('Zadny ORF s CC sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:

                ## lidsky lokus
                ## popis: >hg38_dna range=chr19:17389648-17406217 5'pad=0 3'pad=0 strand=+
                #CCs = [[816,924],[1780,1847],[2028,2089],[2393,2445],[12236,12383],[14371,14468],[14874,15007],[15558,15574]]

                ## kureci lokus
                ## popis: >galGal6_dna range=chr28:3531163-3538110 5'pad=0 3'pad=0 strand=-
                #CCs = [[2122,2314],[2387,2484],[2601,2671],[2741,2802],[3384,3597],[4366,4463],[4545,4678],[4760,4806]]

                ## mysi lokus
                ## popis: >mm10_dna range=chr8:71519771-71538962 5'pad=0 3'pad=0 strand=-
                #CCs = [[1712,1835],[2528,2604],[3137,3177],[4170,4204],[13428,13578],[16406,16503],[16991,17124],[17580,17626]]

                ## luskoun ostrovni
                ## popis:
                ## >ref|NW_023435982.1|:679097-704212 Manis javanica isolate MJ74 unplaced genomic scaffold, YNU_ManJav_2.0 scaffold_88, whole genome shotgun sequence
                #CCs = [[5987,6080],[6935,7011],[7258,7319],[7581,7603],[15298,15439],[16489,16586],[16958,17091],[17551,17591]]

                ## kalon vabivy
                ## popis:
                ## >ref|NW_006429864.1|:923553-939337 Pteropus alecto unplaced genomic scaffold, ASM32557v1 scaffold160, whole genome shotgun sequence
                CCs = [[4140,4263],[5090,5166],[5375,5436],[5686,5717],[10693,10843],[11931,12021]]

                xPccp = []
                yPccp = []

                for iPccp in range(len(CC_all_P_coords)):

                    xPccp += [max_results_CC_P[iPccp][0] + 3*ind_MrkcP_coords[iPccp] + 3]
                    yPccp += [round(max_results_CC_P[iPccp][1]*1000)/10]

                liste_xPccp = xPccp
                liste_yPccp = yPccp

                liste_xPccp = [0] + liste_xPccp
                liste_xPccp = liste_xPccp + [lenbP]

                liste_yPccp = [0] + liste_yPccp
                liste_yPccp = liste_yPccp + [0]

                liste_stringPccp = liste_yPccp

                fig = plt.figure()
                
                fig, ax = plt.subplots() ###

                plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))
                plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, zadane vlakno', fontsize=round(lenbP/65))
                plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbP/65))
                #plt.title('Souhrn ORFu s CC strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxPccp, dyPccp in zip(liste_xPccp, liste_yPccp):
                    if (sxPccp == 0):
                        plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                        plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                    elif (sxPccp > 0) and (sxPccp < lenbP):
                        for sxPccpX, dyPccpX in zip(liste_xPccp[1:-1], liste_yPccp[1:-1]):
                            plt.annotate(''+str(dyPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifPccX = liste_xPccp[1:-1].index(sxPccpX)
                            plt.annotate('     - ORF:'+str(sxPccpX - 3*ind_MrkcP_coords[ifPccX] - 3)+' - CC:'+str(sxPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxPccp == lenbP):
                        plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                        plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                    else:
                        None
                    #plt.scatter(liste_xPccp, liste_stringPccp, s = lenbP/1.5, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xPccp, liste_stringPccp, width = 10, align='center', color='black') # width = lenbP/120

                for iCCs in range(len(CCs)):
                    plt.hlines(y=0,xmin=(CCs[iCCs][0]),xmax=(CCs[iCCs][1]),color='blue',linestyle='-',lw=1200)
                        
                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG-L.svg')
                plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG-L.png', dpi=32) # , dpi=16
                plt.close()

            # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro CC
            # pro BIG verzi - roztazeni obr. pres celou str.
            # dava souradnice ORFu a CC struktury
            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,komplementarni vlakno, vsechny ramce - zacatek:
            # dava souradnice ORFu a CC struktury
            
            if (len(CC_all_N_coords) == 0):
                print('Zadny ORF s CC sekundarni strukturou v komplementarnim vlaknu.')
                print('\n')
            else:
                xNccp = []
                yNccp = []

                for iNccp in range(len(CC_all_N_coords)):

                    xNccp += [max_results_CC_N[iNccp][0] + 3*ind_MrkcN_coords[iNccp] + 3]
                    yNccp += [round(max_results_CC_N[iNccp][1]*1000)/10]

                liste_xNccp = xNccp
                liste_yNccp = yNccp

                liste_xNccp = [0] + liste_xNccp
                liste_xNccp = liste_xNccp + [lenbN]

                liste_yNccp = [0] + liste_yNccp
                liste_yNccp = liste_yNccp + [0]

                liste_stringNccp = liste_yNccp

                fig = plt.figure()
                
                fig, ax = plt.subplots() ###

                plt.figure(figsize=(round(lenbN/120),round(lenbN/75)))
                plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, komplementarni vlakno', fontsize=round(lenbN/65))
                plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbN/65))
                #plt.title('Souhrn ORFu s CC strukturou, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxNccp, dyNccp in zip(liste_xNccp, liste_yNccp):
                    if (sxNccp == 0):
                        plt.annotate(' '+str(dyNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95))
                        plt.annotate('     -  '+str(sxNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95), rotation=90)
                    elif (sxNccp > 0) and (sxNccp < lenbN):
                        for sxNccpX, dyNccpX in zip(liste_xNccp[1:-1], liste_yNccp[1:-1]):
                            plt.annotate(''+str(dyNccpX), xy = (sxNccpX, dyNccpX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifNccX = liste_xNccp[1:-1].index(sxNccpX)
                            plt.annotate('     - ORF:'+str(sxNccpX - 3*ind_MrkcN_coords[ifNccX] - 3)+' - CC:'+str(sxNccpX), xy = (sxNccpX, dyNccpX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxNccp == lenbN):
                        plt.annotate(' '+str(dyNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95))
                        plt.annotate('     -  '+str(sxNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95), rotation=90)
                    else:
                        None
                    #plt.scatter(liste_xNccp, liste_stringNccp, s = lenbN/1.5, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xNccp, liste_stringNccp, width = 10, align='center', color='black') # width = lenbN/120
                        
                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbN/80))
                plt.yticks(size=round(lenbN/80))

                #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno.svg')
                plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno.png', dpi=32) # , dpi=16
                plt.close()

            # dava souradnice ORFu a CC struktury
            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,komplementarni vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,komplementarni vlakno, vsechny ramce - zacatek:
            # dava souradnice ORFu a CC struktury
            # pro BIG verzi - roztazeni obr. pres celou str.
            
            if (len(CC_all_N_coords) == 0):
                print('Zadny ORF s CC sekundarni strukturou v komplementarnim vlaknu.')
                print('\n')
            else:
                xNccp = []
                yNccp = []

                for iNccp in range(len(CC_all_N_coords)):

                    xNccp += [max_results_CC_N[iNccp][0] + 3*ind_MrkcN_coords[iNccp] + 3]
                    yNccp += [round(max_results_CC_N[iNccp][1]*1000)/10]

                liste_xNccp = xNccp
                liste_yNccp = yNccp

                liste_xNccp = [0] + liste_xNccp
                liste_xNccp = liste_xNccp + [lenbN]

                liste_yNccp = [0] + liste_yNccp
                liste_yNccp = liste_yNccp + [0]

                liste_stringNccp = liste_yNccp

                fig = plt.figure()
                
                fig, ax = plt.subplots() ###

                plt.figure(figsize=(round(lenbN/40),round(lenbN/75)))
                plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, komplementarni vlakno', fontsize=round(lenbN/65))
                plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbN/65))
                #plt.title('Souhrn ORFu s CC strukturou, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxNccp, dyNccp in zip(liste_xNccp, liste_yNccp):
                    if (sxNccp == 0):
                        plt.annotate(' '+str(dyNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95))
                        plt.annotate('     -  '+str(sxNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95), rotation=90)
                    elif (sxNccp > 0) and (sxNccp < lenbN):
                        for sxNccpX, dyNccpX in zip(liste_xNccp[1:-1], liste_yNccp[1:-1]):
                            plt.annotate(''+str(dyNccpX), xy = (sxNccpX, dyNccpX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifNccX = liste_xNccp[1:-1].index(sxNccpX)
                            plt.annotate('     - ORF:'+str(sxNccpX - 3*ind_MrkcN_coords[ifNccX] - 3)+' - CC:'+str(sxNccpX), xy = (sxNccpX, dyNccpX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxNccp == lenbN):
                        plt.annotate(' '+str(dyNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95))
                        plt.annotate('     -  '+str(sxNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95), rotation=90)
                    else:
                        None
                    #plt.scatter(liste_xNccp, liste_stringNccp, s = lenbN/1.5, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xNccp, liste_stringNccp, width = 10, align='center', color='black') # width = lenbN/120
                        
                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbN/80))
                plt.yticks(size=round(lenbN/80))

                #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno-BIG.svg')
                plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno-BIG.png', dpi=32) # , dpi=16
                plt.close()

            # pro BIG verzi - roztazeni obr. pres celou str.
            # dava souradnice ORFu a CC struktury
            # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,komplementarni vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - zacatek:
            # upravene - ukaze ORFy a pozici zacatku TM struktury

            if (len(listTMsP_pozice) == 0):
                print('Zadny ORF s TM sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:
                xPtmA = []
                yPtmA = []

                for iPtmA in range(len(listTMsP_pozice)):

                    xPtmA += [count_M_P[iPtmA][0] + 3*apM_P_coords[iPtmA][0] + 3]
                    yPtmA += [count_M_P[iPtmA][1]]

                liste_xPtmA = xPtmA
                liste_yPtmA = yPtmA
                
                liste_xPtmA = [0] + liste_xPtmA
                liste_xPtmA = liste_xPtmA + [lenbP]
                
                liste_yPtmA = [0] + liste_yPtmA
                liste_yPtmA = liste_yPtmA + [0]
                
                liste_stringPtmA = liste_yPtmA

                fig = plt.figure()
                fig, ax = plt.subplots()

                #ax.set_xlim([0,lenbP])
                #fig.add_axes([0,0,lenbP,lenbP])
                
                plt.figure(figsize=(round(lenbP/120),round(lenbP/75)))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxPtmA, dyPtmA in zip(liste_xPtmA, liste_yPtmA):
                 
                    #plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                    if (sxPtmA == 0):
                        plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    elif (sxPtmA > 0) and (sxPtmA < lenbP):
                        for sxPtmAX, dyPtmAX in zip(liste_xPtmA[1:-1], liste_yPtmA[1:-1]):
                            plt.annotate(''+str(dyPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifPX = liste_xPtmA[1:-1].index(sxPtmAX)
                            plt.annotate('     - ORF:'+str(sxPtmAX - 3*apM_P_coords[ifPX][0] - 3)+' - TM:'+str(sxPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxPtmA == lenbP):
                        plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                        
                    #plt.scatter(liste_xPtmA, liste_stringPtmA, s = lenbP/1.8, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xPtmA, liste_stringPtmA, width = 10, align='center', color='black') # width = lenbP/250,

                plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, zadane vlakno', fontsize=round(lenbP/65))
                plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbP/65))
                #plt.title('Souhrn ORFu s TM strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno.svg')
                plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno.png', dpi=32) # , dpi=16
                plt.close()

            # upravene - ukaze ORFy a pozici zacatku TM struktury
            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - zacatek:
            # upravene - ukaze ORFy a pozici zacatku TM struktury
            # pro BIG verzi - roztazeni obr. pres celou str.

            if (len(listTMsP_pozice) == 0):
                print('Zadny ORF s TM sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:
                xPtmA = []
                yPtmA = []

                for iPtmA in range(len(listTMsP_pozice)):

                    xPtmA += [count_M_P[iPtmA][0] + 3*apM_P_coords[iPtmA][0] + 3]
                    yPtmA += [count_M_P[iPtmA][1]]

                liste_xPtmA = xPtmA
                liste_yPtmA = yPtmA
                
                liste_xPtmA = [0] + liste_xPtmA
                liste_xPtmA = liste_xPtmA + [lenbP]
                
                liste_yPtmA = [0] + liste_yPtmA
                liste_yPtmA = liste_yPtmA + [0]
                
                liste_stringPtmA = liste_yPtmA

                fig = plt.figure()
                fig, ax = plt.subplots()

                #ax.set_xlim([0,lenbP])
                #fig.add_axes([0,0,lenbP,lenbP])
                
                plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxPtmA, dyPtmA in zip(liste_xPtmA, liste_yPtmA):
                 
                    #plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                    if (sxPtmA == 0):
                        plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    elif (sxPtmA > 0) and (sxPtmA < lenbP):
                        for sxPtmAX, dyPtmAX in zip(liste_xPtmA[1:-1], liste_yPtmA[1:-1]):
                            plt.annotate(''+str(dyPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifPX = liste_xPtmA[1:-1].index(sxPtmAX)
                            plt.annotate('     - ORF:'+str(sxPtmAX - 3*apM_P_coords[ifPX][0] - 3)+' - TM:'+str(sxPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxPtmA == lenbP):
                        plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                        
                    #plt.scatter(liste_xPtmA, liste_stringPtmA, s = lenbP/1.8, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xPtmA, liste_stringPtmA, width = 10, align='center', color='black') # width = lenbP/250,

                plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, zadane vlakno', fontsize=round(lenbP/65))
                plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbP/65))
                #plt.title('Souhrn ORFu s TM strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG.svg')
                plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG.png', dpi=32) # , dpi=16
                plt.close()

            # pro BIG verzi - roztazeni obr. pres celou str.
            # upravene - ukaze ORFy a pozici zacatku TM struktury
            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - konec.
            
            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - zacatek:
            # upravene - ukaze ORFy a pozici zacatku TM struktury
            # pro BIG verzi - roztazeni obr. pres celou str.
            # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro TM

            if (len(listTMsP_pozice) == 0):
                print('Zadny ORF s TM sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:

                ## lidsky lokus
                ## popis: >hg38_dna range=chr19:17389648-17406217 5'pad=0 3'pad=0 strand=+
                #TMs = [[705,774],[12101,12170]]

                ## kureci lokus
                ## popis: >galGal6_dna range=chr28:3531163-3538110 5'pad=0 3'pad=0 strand=-
                #TMs = [[1990,2059],[3303,3372]]

                ## mysi lokus
                ## popis: >mm10_dna range=chr8:71519771-71538962 5'pad=0 3'pad=0 strand=-
                #TMs = [[1619,1688],[13293,13362]]

                ## luskoun ostrovni
                ## popis:
                ## >ref|NW_023435982.1|:679097-704212 Manis javanica isolate MJ74 unplaced genomic scaffold, YNU_ManJav_2.0 scaffold_88, whole genome shotgun sequence
                #TMs = [[5864,5933],[15151,15220]]

                ## kalon vabivy
                ## popis:
                ## >ref|NW_006429864.1|:923553-939337 Pteropus alecto unplaced genomic scaffold, ASM32557v1 scaffold160, whole genome shotgun sequence
                TMs = [[4047,4116],[10555,10624]]                
                
                xPtmA = []
                yPtmA = []

                for iPtmA in range(len(listTMsP_pozice)):

                    xPtmA += [count_M_P[iPtmA][0] + 3*apM_P_coords[iPtmA][0] + 3]
                    yPtmA += [count_M_P[iPtmA][1]]

                liste_xPtmA = xPtmA
                liste_yPtmA = yPtmA
                
                liste_xPtmA = [0] + liste_xPtmA
                liste_xPtmA = liste_xPtmA + [lenbP]
                
                liste_yPtmA = [0] + liste_yPtmA
                liste_yPtmA = liste_yPtmA + [0]
                
                liste_stringPtmA = liste_yPtmA

                fig = plt.figure()
                fig, ax = plt.subplots()

                #ax.set_xlim([0,lenbP])
                #fig.add_axes([0,0,lenbP,lenbP])
                
                plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxPtmA, dyPtmA in zip(liste_xPtmA, liste_yPtmA):
                 
                    #plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                    if (sxPtmA == 0):
                        plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    elif (sxPtmA > 0) and (sxPtmA < lenbP):
                        for sxPtmAX, dyPtmAX in zip(liste_xPtmA[1:-1], liste_yPtmA[1:-1]):
                            plt.annotate(''+str(dyPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifPX = liste_xPtmA[1:-1].index(sxPtmAX)
                            plt.annotate('     - ORF:'+str(sxPtmAX - 3*apM_P_coords[ifPX][0] - 3)+' - TM:'+str(sxPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxPtmA == lenbP):
                        plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                        
                    #plt.scatter(liste_xPtmA, liste_stringPtmA, s = lenbP/1.8, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xPtmA, liste_stringPtmA, width = 10, align='center', color='black') # width = lenbP/250,

                for iTMs in range(len(TMs)):
                    plt.hlines(y=0,xmin=(TMs[iTMs][0]),xmax=(TMs[iTMs][1]),color='red',linestyle='-',lw=1200)

                plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, zadane vlakno', fontsize=round(lenbP/65))
                plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbP/65))
                #plt.title('Souhrn ORFu s TM strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG-L.svg')
                plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG-L.png', dpi=32) # , dpi=16
                plt.close()

            # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro TM
            # pro BIG verzi - roztazeni obr. pres celou str.
            # upravene - ukaze ORFy a pozici zacatku TM struktury
            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,komplementarni vlakno, vsechny ramce - zacatek:
            # upravene - ukaze ORFy a pozici zacatku TM struktury

            if (len(listTMsN_pozice) == 0):
                print('Zadny ORF s TM sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:
                xNtmA = []
                yNtmA = []

                for iNtmA in range(len(listTMsN_pozice)):

                    xNtmA += [count_M_N[iNtmA][0] + 3*apM_N_coords[iNtmA][0] + 3]
                    yNtmA += [count_M_N[iNtmA][1]]

                liste_xNtmA = xNtmA
                liste_yNtmA = yNtmA
                
                liste_xNtmA = [0] + liste_xNtmA
                liste_xNtmA = liste_xNtmA + [lenbN]
                
                liste_yNtmA = [0] + liste_yNtmA
                liste_yNtmA = liste_yNtmA + [0]
                
                liste_stringNtmA = liste_yNtmA

                fig = plt.figure()
                fig, ax = plt.subplots()

                #ax.set_xlim([0,lenbN])
                #fig.add_axes([0,0,lenbN,lenbN])
                
                plt.figure(figsize=(round(lenbN/120),round(lenbN/75)))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxNtmA, dyNtmA in zip(liste_xNtmA, liste_yNtmA):
                 
                    #plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                    if (sxNtmA == 0):
                        plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    elif (sxNtmA > 0) and (sxNtmA < lenbN):
                        for sxNtmAX, dyNtmAX in zip(liste_xNtmA[1:-1], liste_yNtmA[1:-1]):
                            plt.annotate(''+str(dyNtmAX), xy = (sxNtmAX, dyNtmAX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifNX = liste_xNtmA[1:-1].index(sxNtmAX)
                            plt.annotate('     - ORF:'+str(sxNtmAX - 3*apM_N_coords[ifNX][0] - 3)+' - TM:'+str(sxNtmAX), xy = (sxNtmAX, dyNtmAX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxNtmA == lenbN):
                        plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                        
                    #plt.scatter(liste_xNtmA, liste_stringNtmA, s = lenbN/1.8, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xNtmA, liste_stringNtmA, width = 10, align='center', color='black') # width = lenbN/250,

                plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, komplementarni vlakno', fontsize=round(lenbN/65))
                plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbN/65))
                #plt.title('Souhrn ORFu s TM strukturou, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbN/80))
                plt.yticks(size=round(lenbN/80))

                #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno.svg')
                plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno.png', dpi=32) # , dpi=16
                plt.close()

            # upravene - ukaze ORFy a pozici zacatku TM struktury
            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,komplementarni vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,komplementarni vlakno, vsechny ramce - zacatek:
            # upravene - ukaze ORFy a pozici zacatku TM struktury
            # pro BIG verzi - roztazeni obr. pres celou str.
            
            if (len(listTMsN_pozice) == 0):
                print('Zadny ORF s TM sekundarni strukturou v komplementarnim vlaknu.')
                print('\n')
            else:
                xNtmA = []
                yNtmA = []

                for iNtmA in range(len(listTMsN_pozice)):

                    xNtmA += [count_M_N[iNtmA][0] + 3*apM_N_coords[iNtmA][0] + 3]
                    yNtmA += [count_M_N[iNtmA][1]]

                liste_xNtmA = xNtmA
                liste_yNtmA = yNtmA
                
                liste_xNtmA = [0] + liste_xNtmA
                liste_xNtmA = liste_xNtmA + [lenbN]
                
                liste_yNtmA = [0] + liste_yNtmA
                liste_yNtmA = liste_yNtmA + [0]
                
                liste_stringNtmA = liste_yNtmA

                fig = plt.figure()
                fig, ax = plt.subplots()

                #ax.set_xlim([0,lenbN])
                #fig.add_axes([0,0,lenbN,lenbN])
                
                plt.figure(figsize=(round(lenbN/40),round(lenbN/75)))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxNtmA, dyNtmA in zip(liste_xNtmA, liste_yNtmA):
                 
                    #plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                    if (sxNtmA == 0):
                        plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    elif (sxNtmA > 0) and (sxNtmA < lenbN):
                        for sxNtmAX, dyNtmAX in zip(liste_xNtmA[1:-1], liste_yNtmA[1:-1]):
                            plt.annotate(''+str(dyNtmAX), xy = (sxNtmAX, dyNtmAX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifNX = liste_xNtmA[1:-1].index(sxNtmAX)
                            plt.annotate('     - ORF:'+str(sxNtmAX - 3*apM_N_coords[ifNX][0] - 3)+' - TM:'+str(sxNtmAX), xy = (sxNtmAX, dyNtmAX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxNtmA == lenbN):
                        plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                        
                    #plt.scatter(liste_xNtmA, liste_stringNtmA, s = lenbN/1.8, c ='black') # picker=(True & 10000) # yellow
                    plt.bar(liste_xNtmA, liste_stringNtmA, width = 10, align='center', color='black') # width = lenbN/250,

                plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, komplementarni vlakno', fontsize=round(lenbN/65))
                plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbN/65))
                #plt.title('Souhrn ORFu s TM strukturou, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbN/80))
                plt.yticks(size=round(lenbN/80))

                #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno-BIG.svg')
                plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno-BIG.png', dpi=32) # , dpi=16
                plt.close()

            # pro BIG verzi - roztazeni obr. pres celou str.
            # upravene - ukaze ORFy a pozici zacatku TM struktury
            # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,komplementarni vlakno, vsechny ramce - konec.

            # PRIPRAVA pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane i komplementarni vlakno, vsechny ramce - zacatek:

            import os
            import subprocess

            owd = os.getcwd() # owd = original working directory
            
            os.chdir("/home/oem/Documents/HTMLproTMaCCaGPIrelatall/kohgpi-1.5/")
            #os.system("./kohgpi t")              # neni vubec treba
            #os.system("./mapfill2A >kohgpi.map") # neni vubec treba
            os.system("kohgpi i")
            process_GPI_N = subprocess.Popen(['/home/oem/Documents/HTMLproTMaCCaGPIrelatall/kohgpi-1.5/kohgpi','/home/oem/Documents/HTMLproTMaCCaGPIrelatall/FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa'],
                                       stdout = subprocess.PIPE,
                                       stderr = subprocess.PIPE)
            stdout, stderr = process_GPI_N.communicate()
            stdout

            os.chdir("/home/oem/Documents/HTMLproTMaCCaGPIrelatall/kohgpi-1.5/")
            #os.system("./kohgpi t")              # neni vubec treba
            #os.system("./mapfill2A >kohgpi.map") # neni vubec treba
            os.system("kohgpi i")
            process_GPI_N = subprocess.Popen(['/home/oem/Documents/HTMLproTMaCCaGPIrelatall/kohgpi-1.5/kohgpi','/home/oem/Documents/HTMLproTMaCCaGPIrelatall/FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa'],
                                       stdout = subprocess.PIPE,
                                       stderr = subprocess.PIPE)
            stdout, stderr = process_GPI_N.communicate()
            stdout

            os.chdir(owd) # owd = original working directory

            ###

            from Bio import SeqIO

            if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos')) == False:
                print('Neni *.pos a ...# xcoords.txt soubor, zadane vlakno.')
                GPI_xcoords_P_R = []

            else:

                pos_filename_P = "FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos"

                xcoords_filename_P = "FASTA_output/peptidesP/peptidesP_c00_Standard_Code.xcoords.txt"

                input_handle_P = open(pos_filename_P,"r")

                output_handle_P = open(xcoords_filename_P,"w")

                for seq_record in SeqIO.parse(input_handle_P,"fasta"):
                    output_handle_P.write("%s\n" % (
                        seq_record.id))

                output_handle_P.close()
                input_handle_P.close()

                GPI_xcoords_P = open('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.xcoords.txt')
                GPI_xcoords_P_r = GPI_xcoords_P.read()

                count_GPI_P_n = GPI_xcoords_P_r.count('\n')
                GPI_xcoords_P_rr = GPI_xcoords_P_r.split('\n',count_GPI_P_n)
                GPI_xcoords_P_rrr = GPI_xcoords_P_rr[0:-1]


                iGPiL = []

                for iGP in GPI_xcoords_P_rrr:
                    iGPi = int(iGP)
                    iGPiL += [iGPi]

                GPI_xcoords_P_R = iGPiL

            #print('GPI_xcoords_P_R je: ', GPI_xcoords_P_R)

            ###

            if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos')) == False:
                print('Neni *.pos a ...xcoords.txt soubor, komplementarni vlakno.')
                GPI_xcoords_N_R = []
            else:

                pos_filename_N = "FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos"

                xcoords_filename_N = "FASTA_output/peptidesN/peptidesN_c00_Standard_Code.xcoords.txt"

                input_handle_N = open(pos_filename_N,"r")

                output_handle_N = open(xcoords_filename_N,"w")

                for seq_record in SeqIO.parse(input_handle_N,"fasta"):
                    output_handle_N.write("%s\n" % (
                        seq_record.id))

                output_handle_N.close()
                input_handle_N.close()

                GPI_xcoords_N = open('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.xcoords.txt')
                GPI_xcoords_N_r = GPI_xcoords_N.read()

                count_GPI_N_n = GPI_xcoords_N_r.count('\n')
                GPI_xcoords_N_rr = GPI_xcoords_N_r.split('\n',count_GPI_N_n)
                GPI_xcoords_N_rrr = GPI_xcoords_N_rr[0:-1]


                iGNiL = []

                for iGN in GPI_xcoords_N_rrr:
                    iGNi = int(iGN)
                    iGNiL += [iGNi]

                GPI_xcoords_N_R = iGNiL

            ###

            if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos')) == False:
                print('Neni *.pos a ...ycoords.txt soubor, zadane vlakno.')
            else:

                posy_filename_P = "FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos"

                ycoords_filename_P = "FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt"

                input_handle_Py = open(posy_filename_P,"r")

                output_handle_Py = open(ycoords_filename_P,"w")

                for seq_record in SeqIO.parse(input_handle_Py,"fasta"):
                    output_handle_Py.write("%s\n" % (
                        seq_record.description))

                output_handle_Py.close()
                input_handle_Py.close()

            ###

            if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos')) == False:
                print('Neni *.pos a ...ycoords.txt soubor, komplementarni vlakno.')
            else:

                posy_filename_N = "FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos"

                ycoords_filename_N = "FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt"

                input_handle_Ny = open(posy_filename_N,"r")

                output_handle_Ny = open(ycoords_filename_N,"w")

                for seq_record in SeqIO.parse(input_handle_Ny,"fasta"):
                    output_handle_Ny.write("%s\n" % (
                        seq_record.description))

                output_handle_Ny.close()
                input_handle_Ny.close()

            ###

            import re

            if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt')) == False:
                print('Neni ....ycoords.txt soubor pro ziskani % kvality omega-mista, zadane vlakno.')
                omega_sum_P = []
            else:

                OY_P = open("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt",'r')
                OYR_P = OY_P.readlines()

                omega_sum_P = []

                for line_P in OYR_P:
                    
                    ma_P = re.match(r"(\w+) (\w+) (\W+\w+) (\w+) (\w+\W+\w+) (\w+) (\w+)",line_P)

                    if ma_P == None:
                        omegaprocento_P = 0
                    else:
                        omegaprocento10_P = 10*int(ma_P[7][0])
                        omegaprocento1_P = 1*int(ma_P[7][1])
                        omegaprocento_P = omegaprocento10_P + omegaprocento1_P
                    omega_sum_P += [omegaprocento_P]

            #

            if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt')) == False:
                print('Neni ....ycoords.txt soubor pro ziskani % kvality omega-mista, komplementarni vlakno.')
                omega_sum_N = []
            else:

                OY_N = open("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt",'r')
                OYR_N = OY_N.readlines()

                omega_sum_N = []

                for line_N in OYR_N:
                    
                    ma_N = re.match(r"(\w+) (\w+) (\W+\w+) (\w+) (\w+\W+\w+) (\w+) (\w+)",line_N)

                    if ma_N == None:
                        omegaprocento_N = 0
                    else:
                        omegaprocento10_N = 10*int(ma_N[7][0])
                        omegaprocento1_N = 1*int(ma_N[7][1])
                        omegaprocento_N = omegaprocento10_N + omegaprocento1_N
                    omega_sum_N += [omegaprocento_N]

            # PRIPRAVA pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane i komplementarni vlakno, vsechny ramce - konec.

            # vytvareni mnozin - zacatek:

            # vytvareni mnozin TM, CC, GPI - zadane vlakno:

            set_listTMsP_pozice2_0 = set(listTMsP_pozice2_0)

            set_CC_all_P_coords_0 = set(CC_all_P_coords_0)

            set_GPI_xcoords_P_R = set(GPI_xcoords_P_R)


            # vytvareni mnozin s jen TM, nebo jen s CC, nebo jen s GPI - zadane vlakno:

            jenTMgpi_P = set_listTMsP_pozice2_0.difference(set_CC_all_P_coords_0)
            set_jenTM_P = jenTMgpi_P.difference(set_GPI_xcoords_P_R)
            list_jenTM_P = list(set_jenTM_P) # vysledek - tj. seznam pozic, zacatku ORFu s jen TM strukturou, v zadanem vlaknu

            jenCCgpi_P = set_CC_all_P_coords_0.difference(set_listTMsP_pozice2_0)
            set_jenCC_P = jenCCgpi_P.difference(set_GPI_xcoords_P_R)
            list_jenCC_P = list(set_jenCC_P) # vysledek - tj. seznam pozic, zacatku ORFu s jen CC strukturou, v zadanem vlaknu

            jenGPIcc_P = set_GPI_xcoords_P_R.difference(set_listTMsP_pozice2_0)
            set_jenGPI_P = jenGPIcc_P.difference(set_CC_all_P_coords_0)
            list_jenGPI_P = list(set_jenGPI_P) # vysledek - tj. seznam pozic, zacatku ORFu s jen GPI strukturou, v zadanem vlaknu
            
            # vytvareni mnozin TM, CC, GPI - komplementarni vlakno:

            set_listTMsN_pozice2_0 = set(listTMsN_pozice2_0)

            set_CC_all_N_coords_0 = set(CC_all_N_coords_0)

            set_GPI_xcoords_N_R = set(GPI_xcoords_N_R)


            # vytvareni mnozin s jen TM, nebo jen s CC, nebo jen s GPI - komplementarni vlakno:

            jenTMgpi_N = set_listTMsN_pozice2_0.difference(set_CC_all_N_coords_0)
            set_jenTM_N = jenTMgpi_N.difference(set_GPI_xcoords_N_R)
            list_jenTM_N = list(set_jenTM_N) # vysledek - tj. seznam pozic, zacatku ORFu s jen TM strukturou, v komplementarnim vlaknu


            jenCCgpi_N = set_CC_all_N_coords_0.difference(set_listTMsN_pozice2_0)
            set_jenCC_N = jenCCgpi_N.difference(set_GPI_xcoords_N_R)
            list_jenCC_N = list(set_jenCC_N) # vysledek - tj. seznam pozic, zacatku ORFu s jen CC strukturou, v komplementarnim vlaknu

            jenGPIcc_N = set_GPI_xcoords_N_R.difference(set_listTMsN_pozice2_0)
            set_jenGPI_N = jenGPIcc_N.difference(set_CC_all_N_coords_0)
            list_jenGPI_N = list(set_jenGPI_N) # vysledek - tj. seznam pozic, zacatku ORFu s jen GPI strukturou, v komplementarnim vlaknu


            # vytvareni pruniku vzdy jen dvou mnozin - zadane vlakno:
            # napr. kdyz chci TM a CC, musim udelat jejich prunik, ale odecist od tohoto pruniku, spolecne prvky s GPI

            set_TMaCCsGPI_P = set_listTMsP_pozice2_0.intersection(set_CC_all_P_coords_0)
            set_TMaCC_P = set_TMaCCsGPI_P.difference(set_GPI_xcoords_P_R)
            list_TMaCC_P = list(set_TMaCC_P)

            set_TMaGPIsCC_P = set_listTMsP_pozice2_0.intersection(set_GPI_xcoords_P_R)
            set_TMaGPI_P = set_TMaGPIsCC_P.difference(set_CC_all_P_coords_0)
            list_TMaGPI_P = list(set_TMaGPI_P)

            set_GPIaCCsTM_P = set_GPI_xcoords_P_R.intersection(set_CC_all_P_coords_0)
            set_GPIaCC_P = set_GPIaCCsTM_P.difference(set_listTMsP_pozice2_0)
            list_GPIaCC_P = list(set_GPIaCC_P)


            # vytvareni pruniku vzdy jen dvou mnozin - komplementarni vlakno:
            # napr. kdyz chci TM a CC, musim udelat jejich prunik, ale odecist od tohoto pruniku, spolecne prvky s GPI

            set_TMaCCsGPI_N = set_listTMsN_pozice2_0.intersection(set_CC_all_N_coords_0)
            set_TMaCC_N = set_TMaCCsGPI_N.difference(set_GPI_xcoords_N_R)
            list_TMaCC_N = list(set_TMaCC_N)

            set_TMaGPIsCC_N = set_listTMsN_pozice2_0.intersection(set_GPI_xcoords_N_R)
            set_TMaGPI_N = set_TMaGPIsCC_N.difference(set_CC_all_N_coords_0)
            list_TMaGPI_N = list(set_TMaGPI_N)

            set_GPIaCCsTM_N = set_GPI_xcoords_N_R.intersection(set_CC_all_N_coords_0)
            set_GPIaCC_N = set_GPIaCCsTM_N.difference(set_listTMsN_pozice2_0)
            list_GPIaCC_N = list(set_GPIaCC_N)

            
            # vytvareni pruniku vsech 3 (TM,CC,GPI) mnozin:
            # zadane vlakno:
            set_TMaCCxx_P = set_listTMsP_pozice2_0.intersection(set_CC_all_P_coords_0)
            set_TMaCCaGPI_P = set_GPI_xcoords_P_R.intersection(set_TMaCCxx_P)
            list_TMaCCaGPI_P = list(set_TMaCCaGPI_P)

            # komplementarni vlakno:
            set_TMaCCxx_N = set_listTMsN_pozice2_0.intersection(set_CC_all_N_coords_0)
            set_TMaCCaGPI_N = set_GPI_xcoords_N_R.intersection(set_TMaCCxx_N)
            list_TMaCCaGPI_N = list(set_TMaCCaGPI_N)
                
            # vytvareni mnozin - konec.

            # ziskavani obou souradnic pro GPI strukturu - zadane i komplementarni vlakno - zacatek:

            GPI_XYcoords_P_R = []

            for iLOP in listORFsP_pozice:
                for iGPP in GPI_xcoords_P_R:
                    if iGPP == iLOP[0]:
                        GPI_XYcoords_P_R += [iLOP]
                    else:
                        continue

            GPI_XYcoords_N_R = []

            for iLON in listORFsN_pozice:
                for iGPN in GPI_xcoords_N_R:
                    if iGPN == iLON[0]:
                        GPI_XYcoords_N_R += [iLON]
                    else:
                        continue

            # ziskavani obou souradnic pro GPI strukturu - zadane i komplementarni vlakno - konec.

            # ziskavani souradnic pro ORFy jen s TM, zadane i komplementarni vlakno - zacatek:

            list_jenTM_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLOP in listORFsP_pozice:
                for iTPP in list_jenTM_P:
                    if iTPP == iLOP[0]:
                        list_jenTM_P_01 += [iLOP]
                    else:
                        continue
            # vysledek: list_jenTM_P_01

            list_jenTM_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLON in listORFsN_pozice:
                for iTPN in list_jenTM_N:
                    if iTPN == iLON[0]:
                        list_jenTM_N_01 += [iLON]
                    else:
                        continue
            # vysledek: list_jenTM_N_01

            # ziskavani souradnic pro ORFy jen s TM, zadane i komplementarni vlakno - konec.

            # ziskavani souradnic pro ORFy jen s CC, zadane i komplementarni vlakno - zacatek:

            list_jenCC_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLOP in listORFsP_pozice:
                for iCPP in list_jenCC_P:
                    if iCPP == iLOP[0]:
                        list_jenCC_P_01 += [iLOP]
                    else:
                        continue
            # vysledek: list_jenCC_P_01

            list_jenCC_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLON in listORFsN_pozice:
                for iCPN in list_jenCC_N:
                    if iCPN == iLON[0]:
                        list_jenCC_N_01 += [iLON]
                    else:
                        continue
            # vysledek: list_jenCC_N_01

            # ziskavani souradnic pro ORFy jen s CC, zadane i komplementarni vlakno - konec.

            # ziskavani souradnic pro ORFy jen s GPI, zadane i komplementarni vlakno - zacatek:

            list_jenGPI_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLOP in listORFsP_pozice:
                for iGPP in list_jenGPI_P:
                    if iGPP == iLOP[0]:
                        list_jenGPI_P_01 += [iLOP]
                    else:
                        continue
            # vysledek: list_jenGPI_P_01

            list_jenGPI_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLON in listORFsN_pozice:
                for iGPN in list_jenGPI_N:
                    if iGPN == iLON[0]:
                        list_jenGPI_N_01 += [iLON]
                    else:
                        continue
            # vysledek: list_jenGPI_N_01

            # ziskavani souradnic pro ORFy jen s GPI, zadane i komplementarni vlakno - konec.

            # ziskavani souradnic pro ORFy jen s TM a zaroven s CC, zadane i komplementarni vlakno - zacatek:

            list_TMaCC_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLOP in listORFsP_pozice:
                for iTCPP in list_TMaCC_P:
                    if iTCPP == iLOP[0]:
                        list_TMaCC_P_01 += [iLOP]
                    else:
                        continue
            # vysledek: list_TMaCC_P_01

            list_TMaCC_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLON in listORFsN_pozice:
                for iTCPN in list_TMaCC_N:
                    if iTCPN == iLON[0]:
                        list_TMaCC_N_01 += [iLON]
                    else:
                        continue
            # vysledek: list_TMaCC_N_01

            # ziskavani souradnic pro ORFy jen s TM a zaroven s CC, zadane i komplementarni vlakno - konec.
            
            # ziskavani souradnic pro ORFy jen s TM a zaroven s GPI, zadane i komplementarni vlakno - zacatek:

            list_TMaGPI_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLOP in listORFsP_pozice:
                for iTGPP in list_TMaGPI_P:
                    if iTGPP == iLOP[0]:
                        list_TMaGPI_P_01 += [iLOP]
                    else:
                        continue
            # vysledek: list_TMaGPI_P_01

            list_TMaGPI_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLON in listORFsN_pozice:
                for iTGPN in list_TMaGPI_N:
                    if iTGPN == iLON[0]:
                        list_TMaGPI_N_01 += [iLON]
                    else:
                        continue
            # vysledek: list_TMaGPI_N_01

            # ziskavani souradnic pro ORFy jen s TM a zaroven s GPI, zadane i komplementarni vlakno - konec.

            # ziskavani souradnic pro ORFy jen s CC a zaroven s GPI, zadane i komplementarni vlakno - zacatek:

            list_GPIaCC_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLOP in listORFsP_pozice:
                for iGCPP in list_GPIaCC_P:
                    if iGCPP == iLOP[0]:
                        list_GPIaCC_P_01 += [iLOP]
                    else:
                        continue
            # vysledek: list_GPIaCC_P_01

            list_GPIaCC_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLON in listORFsN_pozice:
                for iGCPN in list_GPIaCC_N:
                    if iGCPN == iLON[0]:
                        list_GPIaCC_N_01 += [iLON]
                    else:
                        continue
            # vysledek: list_GPIaCC_N_01

            # ziskavani souradnic pro ORFy jen s CC a zaroven s GPI, zadane i komplementarni vlakno - konec.

            # ziskavani souradnic pro ORFy s TM a zaroven s CC a zaroven s GPI, zadane i komplementarni vlakno - zacatek:

            list_TMaCCaGPI_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLOP in listORFsP_pozice:
                for iTCGPP in list_TMaCCaGPI_P:
                    if iTCGPP == iLOP[0]:
                        list_TMaCCaGPI_P_01 += [iLOP]
                    else:
                        continue
            # vysledek: list_TMaCCaGPI_P_01

            list_TMaCCaGPI_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

            for iLON in listORFsN_pozice:
                for iTCGPN in list_TMaCCaGPI_N:
                    if iTCGPN == iLON[0]:
                        list_TMaCCaGPI_N_01 += [iLON]
                    else:
                        continue
            # vysledek: list_TMaCCaGPI_N_01

            # ziskavani souradnic pro ORFy s TM a zaroven s CC a zaroven s GPI, zadane i komplementarni vlakno - konec.

            # "linkove" obrazky pro CC (P i N) - zacatek:

            for ilPcc in range(len(CC_all_P_coords)):

                # "linkovy obr.":
                xsP = np.linspace(1,1,lenbP)
                plt.figure(figsize=(lenbP/120,lenbP/4800))
                plt.hlines(y=1,xmin=1,xmax=lenbP,color='blue',linestyle='-',lw=20)
                plt.hlines(y=1,xmin=(CC_all_P_coords[ilPcc][0]),xmax=(CC_all_P_coords[ilPcc][1]),color='r',linestyle='-',lw=30)
                plt.yticks([])
                plt.xticks([])
                '''
                if ((listORFsP_pozice[irpP_c00_4fA][0])%3 == 0):
                    frameP_ORFall_lin = 1
                elif ((listORFsP_pozice[irpP_c00_4fA][0] + 1)%3 == 0):
                    frameP_ORFall_lin = 3
                elif ((listORFsP_pozice[irpP_c00_4fA][0] + 2)%3 == 0):
                    frameP_ORFall_lin = 2
                else:
                    continue
                '''
                #plt.annotate(str(listORFsP_pozice[irpP_c00_4fA][0]),(listORFsP_pozice[irpP_c00_4fA][0],1),fontsize=16)
                #plt.annotate(str(listORFsP_pozice[irpP_c00_4fA][1]),(listORFsP_pozice[irpP_c00_4fA][1],1),fontsize=16)
                #plt.annotate(str(lenbP),(lenbP,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
                #plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
                '''
                plt.annotate('Fr.:'+str(frameP_ORFall_lin),((listORFsP_pozice[irpP_c00_4fA][0]+listORFsP_pozice[irpP_c00_4fA][1])/2,1),fontsize=16,color='green',verticalalignment='top')
                '''
                #plt.savefig('/home/oem/Documents/_html/figures/NA string - fig No. ' + str([irpP_c00_4fA]) + ' provided_string.png')
                plt.savefig('figures/given/NA string - fig No. ' + str(CC_all_P_coords[ilPcc]) + ' provided_string.png')
                #plt.show()
                plt.close()#

            for ilNcc in range(len(CC_all_N_coords)):

                # "linkovy obr.":
                xsN = np.linspace(1,1,lenbN)
                plt.figure(figsize=(lenbN/120,lenbN/4800))
                plt.hlines(y=1,xmin=1,xmax=lenbN,color='blue',linestyle='-',lw=20)
                plt.hlines(y=1,xmin=(CC_all_N_coords[ilNcc][0]),xmax=(CC_all_N_coords[ilNcc][1]),color='r',linestyle='-',lw=30)
                plt.yticks([])
                plt.xticks([])
                '''
                if ((listORFsN_pozice[irpN_c00_4fA][0])%3 == 0):
                    frameN_ORFall_lin = 1
                elif ((listORFsN_pozice[irpN_c00_4fA][0] + 1)%3 == 0):
                    frameN_ORFall_lin = 3
                elif ((listORFsN_pozice[irpN_c00_4fA][0] + 2)%3 == 0):
                    frameN_ORFall_lin = 2
                else:
                    continue
                '''
                #plt.annotate(str(listORFsN_pozice[irpN_c00_4fA][0]),(listORFsN_pozice[irpN_c00_4fA][0],1),fontsize=16)
                #plt.annotate(str(listORFsN_pozice[irpN_c00_4fA][1]),(listORFsN_pozice[irpN_c00_4fA][1],1),fontsize=16)
                #plt.annotate(str(lenbN),(lenbN,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
                #plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
                '''
                plt.annotate('Fr.:'+str(frameN_ORFall_lin),((listORFsN_pozice[irpN_c00_4fA][0]+listORFsN_pozice[irpN_c00_4fA][1])/2,1),fontsize=16,color='green',verticalalignment='top')
                '''
                #plt.savefig('/home/oem/Documents/_html/figures/NA string - fig No. ' + str([irpN_c00_4fA]) + ' complementary_string.png')
                plt.savefig('figures/complementary/NA string - fig No. ' + str(CC_all_N_coords[ilNcc]) + ' complementary_string.png')
                #plt.show()
                plt.close()#

            # "linkove" obrazky pro CC (P i N) - konec.

            # "linkove" obrazky pro GPI (P i N) - zacatek:

            for ilPgpi in range(len(GPI_XYcoords_P_R)):

                # "linkovy obr.":
                xsP = np.linspace(1,1,lenbP)
                plt.figure(figsize=(lenbP/120,lenbP/4800))
                plt.hlines(y=1,xmin=1,xmax=lenbP,color='blue',linestyle='-',lw=20)
                plt.hlines(y=1,xmin=(GPI_XYcoords_P_R[ilPgpi][0]),xmax=(GPI_XYcoords_P_R[ilPgpi][1]),color='r',linestyle='-',lw=30)
                
                plt.yticks([])
                plt.xticks([])

                plt.savefig('figures/given/NA string - fig No. ' + str(GPI_XYcoords_P_R[ilPgpi]) + ' provided_string.png')
                #plt.show()
                plt.close()

            for ilNgpi in range(len(GPI_XYcoords_N_R)):

                # "linkovy obr.":
                xsN = np.linspace(1,1,lenbN)
                plt.figure(figsize=(lenbN/120,lenbN/4800))
                plt.hlines(y=1,xmin=1,xmax=lenbN,color='blue',linestyle='-',lw=20)
                plt.hlines(y=1,xmin=(GPI_XYcoords_N_R[ilNgpi][0]),xmax=(GPI_XYcoords_N_R[ilNgpi][1]),color='r',linestyle='-',lw=30)
                
                plt.yticks([])
                plt.xticks([])

                plt.savefig('figures/complementary/NA string - fig No. ' + str(GPI_XYcoords_N_R[ilNgpi]) + ' complementary_string.png')
                #plt.show()
                plt.close()

            # "linkove" obrazky pro GPI (P i N) - konec.

            # vytvareni slovniku pro vykresleni popisu GPI mista v tabulce v HTML - zacatek: - NEPOUZITO v HTML

            if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos')) == False:
                print('Neni *.pos soubor pro tvorbu slovniku pro GPI strukturu, zadane vlakno. - stejne je to NEpouzito v HTML')
            else:
                record_dict_P = SeqIO.to_dict(SeqIO.parse("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos","fasta"))

            if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos')) == False:
                print('Neni *.pos soubor pro tvorbu slovniku pro GPI strukturu, komplementarni vlakno.- stejne je to NEpouzito v HTML')
            else:
                record_dict_N = SeqIO.to_dict(SeqIO.parse("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos","fasta"))

            # vytvareni slovniku pro vykresleni popisu GPI mista v tabulce v HTML - konec. - NEPOUZITO v HTML

            # vytvareni obrazku pro GPI sestrihove misto - zadane vlakno, vsechny ramce - zacatek:

            import re
            from Bio import SeqIO

            if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt')) == False:
                print('Neni ...ycoords.txt soubor pro tvorbu obrazku pro GPI strukturu, zadane vlakno, vsechny ramce.')
            else:

                oGPIy_P = open("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt",'r')
                GPIy_P = oGPIy_P.readlines()

                distances_AAs_P = []
                cleavages_P = []
                XYcoords_P = []

                for line_gpi_P in GPIy_P:

                    ma_XYcoords_P = re.match(r"(\w+) (\w+)",line_gpi_P).groups()
                    ma_Xcoord_P = int(ma_XYcoords_P[0])
                    ma_Ycoord_P = int(ma_XYcoords_P[1])

                    XYcoords_P += [[ma_Xcoord_P, ma_Ycoord_P]]

                    distance_AAs_P = round((ma_Ycoord_P - ma_Xcoord_P)/3)

                    distances_AAs_P += [distance_AAs_P]

                    ma_Pc = re.match(r"(\w+) (\w+) (\W+\w+) (\w+) (\w+\W+\w+) (\w+) (\w+)",line_gpi_P)

                    if ma_Pc == None:
                        cleavage_int_P = distance_AAs_P
                        cleavages_P += [cleavage_int_P]
                    else:
                        cleavage_P = ma_Pc[5][2:]
                        cleavage_int_P = round(int(cleavage_P))
                        cleavages_P += [cleavage_int_P]
                            
                import matplotlib.pyplot as plt

                for iGPIobrP in range(len(omega_sum_P)):

                    fig = plt.figure()
                    fig, ax = plt.subplots()
                    #plt.figure(figsize=(lenbP/200,lenbP/256))

                    ax.set_xlim([0,distances_AAs_P[iGPIobrP]])
                    ax.set_ylim([0,100])

                    if ( distances_AAs_P[iGPIobrP] == cleavages_P[iGPIobrP] ):

                        plt.bar(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], omega_sum_P[iGPIobrP], width = 1, align='center', color='green')
                        plt.annotate(omega_sum_P[iGPIobrP], xy=(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], omega_sum_P[iGPIobrP]) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'bottom')
                        #plt.annotate('None', xy=(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                        plt.annotate('Neurceno omega-misto', xy=( (distances_AAs_P[iGPIobrP])/2,0 ) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                        #plt.annotate('None', xy=(distances_AAs_P[iGPIobrP],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                        
                    else:
                        
                        plt.bar(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], omega_sum_P[iGPIobrP], width = 1, align='center', color='green')
                        plt.annotate(omega_sum_P[iGPIobrP], xy=(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], omega_sum_P[iGPIobrP]) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], xy=(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')

                    plt.xlabel('Poradi AMK (N-C konec) a omega-misto (-)', fontsize=16)
                    plt.ylabel('Kvalita omega-mista (%)', fontsize=16)
                    plt.title('ORF ' + str(XYcoords_P[iGPIobrP]) +' s GPI strukturou, zadane vlakno.\n', fontsize=14)

                    plt.savefig('figures/GPI/P' + str(XYcoords_P[iGPIobrP]) + '.png')

                    #plt.show()
                    plt.close()

            # vytvareni obrazku pro GPI sestrihove misto - zadane vlakno, vsechny ramce - konec.

            # vytvareni obrazku pro GPI sestrihove misto - komplementarni vlakno, vsechny ramce - zacatek:

            #import re
            #from Bio import SeqIO

            if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt')) == False:
                print('Neni ...ycoords.txt soubor pro tvorbu obrazku pro GPI strukturu, komplementarni vlakno, vsechny ramce.')
            else:

                oGPIy_N = open("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt",'r')
                GPIy_N = oGPIy_N.readlines()

                distances_AAs_N = []
                cleavages_N = []
                XYcoords_N = []

                for line_gpi_N in GPIy_N:

                    ma_XYcoords_N = re.match(r"(\w+) (\w+)",line_gpi_N).groups()
                    ma_Xcoord_N = int(ma_XYcoords_N[0])
                    ma_Ycoord_N = int(ma_XYcoords_N[1])

                    XYcoords_N += [[ma_Xcoord_N, ma_Ycoord_N]]

                    distance_AAs_N = round((ma_Ycoord_N - ma_Xcoord_N)/3)

                    distances_AAs_N += [distance_AAs_N]

                    ma_Nc = re.match(r"(\w+) (\w+) (\W+\w+) (\w+) (\w+\W+\w+) (\w+) (\w+)",line_gpi_N)

                    if ma_Nc == None:
                        cleavage_int_N = distance_AAs_N
                        cleavages_N += [cleavage_int_N]
                    else:
                        cleavage_N = ma_Nc[5][2:]
                        cleavage_int_N = round(int(cleavage_N))
                        cleavages_N += [cleavage_int_N]
                            
                #import matplotlib.pyplot as plt

                for iGPIobrN in range(len(omega_sum_N)):

                    fig = plt.figure()
                    fig, ax = plt.subplots()
                    #plt.figure(figsize=(lenbN/200,lenbN/256))

                    ax.set_xlim([0,distances_AAs_N[iGPIobrN]])
                    ax.set_ylim([0,100])

                    if ( distances_AAs_N[iGPIobrN] == cleavages_N[iGPIobrN] ):

                        plt.bar(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], omega_sum_N[iGPIobrN], width = 1, align='center', color='green')
                        plt.annotate(omega_sum_N[iGPIobrN], xy=(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], omega_sum_N[iGPIobrN]) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'bottom')
                        #plt.annotate('None', xy=(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                        plt.annotate('Neurceno omega-misto', xy=( (distances_AAs_N[iGPIobrN])/2,0 ) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                        #plt.annotate('None', xy=(distances_AAs_N[iGPIobrN],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                        
                    else:
                        
                        plt.bar(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], omega_sum_N[iGPIobrN], width = 1, align='center', color='green')
                        plt.annotate(omega_sum_N[iGPIobrN], xy=(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], omega_sum_N[iGPIobrN]) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], xy=(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')

                    plt.xlabel('Poradi AMK (N-C konec) a omega-misto (-)', fontsize=16)
                    plt.ylabel('Kvalita omega-mista (%)', fontsize=16)
                    plt.title('ORF ' + str(XYcoords_N[iGPIobrN]) +' s GPI strukturou, zadane vlakno.\n', fontsize=14)

                    plt.savefig('figures/GPI/N/' + str(XYcoords_N[iGPIobrN]) + '.png')

                    #plt.show()
                    plt.close()

            # vytvareni obrazku pro GPI sestrihove misto - komplementarni vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - zacatek:
            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist

            if (len(omega_sum_P) == 0):
                print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:
                xPgpi = []
                for ixGPI_P in range(len(GPI_xcoords_P_R)):
                    xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P])]
                    '''
                    if (distances_AAs_P[ixGPI_P] == cleavages_P[ixGPI_P]):
                        xPgpi += [GPI_xcoords_P_R[ixGPI_P]]
                    elif (distances_AAs_P[ixGPI_P] > cleavages_P[ixGPI_P]):
                        xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P]) + 3]
                    else:
                        None
                    '''
                xPgpi
                yPgpi = omega_sum_P
                
                liste_xPgpi = xPgpi
                liste_yPgpi = yPgpi

                liste_xPgpi = [0] + liste_xPgpi
                liste_xPgpi = liste_xPgpi + [lenbP]

                liste_yPgpi = [0] + liste_yPgpi
                liste_yPgpi = liste_yPgpi + [0]

                liste_stringPgpi = liste_yPgpi

                fig = plt.figure()
                
                fig, ax = plt.subplots()

                plt.figure(figsize=(round(lenbP/120),round(lenbP/75)))
                plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, zadane vlakno', fontsize=round(lenbP/65))
                plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbP/65))
                #plt.title('Souhrn ORFu s omega-misty, zadane vlakno\n\n', fontsize=round(lenbP/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxPgpi, dyPgpi in zip(liste_xPgpi, liste_yPgpi):
                 
                    if (sxPgpi == 0):
                        plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    elif ((sxPgpi > 0) and (sxPgpi < lenbP)):
                        for sxPgpiX, dyPgpiX in zip(liste_xPgpi[1:-1], liste_yPgpi[1:-1]):
                            plt.annotate(''+str(dyPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifPgpiX = liste_xPgpi[1:-1].index(sxPgpiX)
                            plt.annotate('     - ORF:'+str(sxPgpiX - 3*(distances_AAs_P[ifPgpiX] - cleavages_P[ifPgpiX]))+' - GPI:'+str(sxPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxPgpi == lenbP):
                        plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                    
                    #plt.scatter(liste_xPgpi, liste_stringPgpi, s = lenbP/1.8, c ='black')
                    plt.bar(liste_xPgpi, liste_stringPgpi, width = 10,  align='center', color='black') # width = lenbP/250, 

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno.svg')
                plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno.png', dpi=32) # , dpi=16
                plt.close()

            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - zacatek:
            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
            # pro BIG verzi - roztazeni obr. pres celou str.
            
            if (len(omega_sum_P) == 0):
                print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:
                xPgpi = []
                for ixGPI_P in range(len(GPI_xcoords_P_R)):
                    xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P])]
                    '''
                    if (distances_AAs_P[ixGPI_P] == cleavages_P[ixGPI_P]):
                        xPgpi += [GPI_xcoords_P_R[ixGPI_P]]
                    elif (distances_AAs_P[ixGPI_P] > cleavages_P[ixGPI_P]):
                        xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P]) + 3]
                    else:
                        None
                    '''
                xPgpi
                yPgpi = omega_sum_P
                
                liste_xPgpi = xPgpi
                liste_yPgpi = yPgpi

                liste_xPgpi = [0] + liste_xPgpi
                liste_xPgpi = liste_xPgpi + [lenbP]

                liste_yPgpi = [0] + liste_yPgpi
                liste_yPgpi = liste_yPgpi + [0]

                liste_stringPgpi = liste_yPgpi

                fig = plt.figure()
                
                fig, ax = plt.subplots()

                plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))
                plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, zadane vlakno', fontsize=round(lenbP/65))
                plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbP/65))
                #plt.title('Souhrn ORFu s omega-misty, zadane vlakno\n\n', fontsize=round(lenbP/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxPgpi, dyPgpi in zip(liste_xPgpi, liste_yPgpi):
                 
                    if (sxPgpi == 0):
                        plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    elif ((sxPgpi > 0) and (sxPgpi < lenbP)):
                        for sxPgpiX, dyPgpiX in zip(liste_xPgpi[1:-1], liste_yPgpi[1:-1]):
                            plt.annotate(''+str(dyPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifPgpiX = liste_xPgpi[1:-1].index(sxPgpiX)
                            plt.annotate('     - ORF:'+str(sxPgpiX - 3*(distances_AAs_P[ifPgpiX] - cleavages_P[ifPgpiX]))+' - GPI:'+str(sxPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxPgpi == lenbP):
                        plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                    
                    #plt.scatter(liste_xPgpi, liste_stringPgpi, s = lenbP/1.8, c ='black')
                    plt.bar(liste_xPgpi, liste_stringPgpi, width = 10,  align='center', color='black') # width = lenbP/250, 

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG.svg')
                plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG.png', dpi=32) # , dpi=16
                plt.close()

            # pro BIG verzi - roztazeni obr. pres celou str.
            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - zacatek:
            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
            # pro BIG verzi - roztazeni obr. pres celou str.
            # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro GPI
            
            if (len(omega_sum_P) == 0):
                print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:

                ## lidsky lokus
                ## popis: >hg38_dna range=chr19:17389648-17406217 5'pad=0 3'pad=0 strand=+
                #GPI = 2463

                ## kureci lokus
                ## popis: >galGal6_dna range=chr28:3531163-3538110 5'pad=0 3'pad=0 strand=-
                #GPI = 2817

                ## mysi lokus
                ## popis: >mm10_dna range=chr8:71519771-71538962 5'pad=0 3'pad=0 strand=-
                #GPI = 4210

                ## luskoun ostrovni
                ## popis:
                ## >ref|NW_023435982.1|:679097-704212 Manis javanica isolate MJ74 unplaced genomic scaffold, YNU_ManJav_2.0 scaffold_88, whole genome shotgun sequence
                #GPI = 7606

                ## kalon vabivy
                ## popis:
                ## >ref|NW_006429864.1|:923553-939337 Pteropus alecto unplaced genomic scaffold, ASM32557v1 scaffold160, whole genome shotgun sequence
                GPI = 5720
                
                xPgpi = []
                for ixGPI_P in range(len(GPI_xcoords_P_R)):
                    xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P])]
                    '''
                    if (distances_AAs_P[ixGPI_P] == cleavages_P[ixGPI_P]):
                        xPgpi += [GPI_xcoords_P_R[ixGPI_P]]
                    elif (distances_AAs_P[ixGPI_P] > cleavages_P[ixGPI_P]):
                        xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P]) + 3]
                    else:
                        None
                    '''
                xPgpi
                yPgpi = omega_sum_P
                
                liste_xPgpi = xPgpi
                liste_yPgpi = yPgpi

                liste_xPgpi = [0] + liste_xPgpi
                liste_xPgpi = liste_xPgpi + [lenbP]

                liste_yPgpi = [0] + liste_yPgpi
                liste_yPgpi = liste_yPgpi + [0]

                liste_stringPgpi = liste_yPgpi

                fig = plt.figure()
                
                fig, ax = plt.subplots()

                plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))
                plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, zadane vlakno', fontsize=round(lenbP/65))
                plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbP/65))
                #plt.title('Souhrn ORFu s omega-misty, zadane vlakno\n\n', fontsize=round(lenbP/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxPgpi, dyPgpi in zip(liste_xPgpi, liste_yPgpi):
                 
                    if (sxPgpi == 0):
                        plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    elif ((sxPgpi > 0) and (sxPgpi < lenbP)):
                        for sxPgpiX, dyPgpiX in zip(liste_xPgpi[1:-1], liste_yPgpi[1:-1]):
                            plt.annotate(''+str(dyPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifPgpiX = liste_xPgpi[1:-1].index(sxPgpiX)
                            plt.annotate('     - ORF:'+str(sxPgpiX - 3*(distances_AAs_P[ifPgpiX] - cleavages_P[ifPgpiX]))+' - GPI:'+str(sxPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxPgpi == lenbP):
                        plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                    
                    #plt.scatter(liste_xPgpi, liste_stringPgpi, s = lenbP/1.8, c ='black')
                    plt.bar(liste_xPgpi, liste_stringPgpi, width = 10,  align='center', color='black') # width = lenbP/250, 

                wth = 400 # tj. width
                hht = 4 # tj. height
                xy = (GPI-0.5*wth,0-0.5*hht)
                GPI_rectangle = plt.Rectangle(xy,wth,hht,angle=0.0,color='yellow')
                fig = plt.gcf()
                ax = fig.gca()
                ax.add_patch(GPI_rectangle)

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbP/80))
                plt.yticks(size=round(lenbP/80))

                #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG-L.svg')
                plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG-L.png', dpi=32) # , dpi=16
                plt.close()

            # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro GPI
            # pro BIG verzi - roztazeni obr. pres celou str.
            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,komplementarni vlakno, vsechny ramce - zacatek:
            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist

            if (len(omega_sum_N) == 0):
                print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:
                xNgpi = []
                for ixGPI_N in range(len(GPI_xcoords_N_R)):
                    xNgpi += [GPI_xcoords_N_R[ixGPI_N] + 3*(distances_AAs_N[ixGPI_N] - cleavages_N[ixGPI_N])]
                    '''
                    if (distances_AAs_N[ixGPI_N] == cleavages_N[ixGPI_N]):
                        xNgpi += [GPI_xcoords_N_R[ixGPI_N]]
                    elif (distances_AAs_N[ixGPI_N] > cleavages_N[ixGPI_N]):
                        xNgpi += [GPI_xcoords_N_R[ixGPI_N] + 3*(distances_AAs_N[ixGPI_N] - cleavages_N[ixGPI_N]) + 3]
                    else:
                        None
                    '''
                xNgpi
                yNgpi = omega_sum_N
                
                liste_xNgpi = xNgpi
                liste_yNgpi = yNgpi

                liste_xNgpi = [0] + liste_xNgpi
                liste_xNgpi = liste_xNgpi + [lenbN]

                liste_yNgpi = [0] + liste_yNgpi
                liste_yNgpi = liste_yNgpi + [0]

                liste_stringNgpi = liste_yNgpi

                fig = plt.figure()
                
                fig, ax = plt.subplots()

                plt.figure(figsize=(round(lenbN/120),round(lenbN/75)))
                plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, komplementarni vlakno', fontsize=round(lenbN/65))
                plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbN/65))
                #plt.title('Souhrn ORFu s omega-misty, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxNgpi, dyNgpi in zip(liste_xNgpi, liste_yNgpi):
                 
                    if (sxNgpi == 0):
                        plt.annotate(''+str(dyNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    elif ((sxNgpi > 0) and (sxNgpi < lenbN)):
                        for sxNgpiX, dyNgpiX in zip(liste_xNgpi[1:-1], liste_yNgpi[1:-1]):
                            plt.annotate(''+str(dyNgpiX), xy = (sxNgpiX, dyNgpiX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifNgpiX = liste_xNgpi[1:-1].index(sxNgpiX)
                            plt.annotate('     - ORF:'+str(sxNgpiX - 3*(distances_AAs_N[ifNgpiX] - cleavages_N[ifNgpiX]))+' - GPI:'+str(sxNgpiX), xy = (sxNgpiX, dyNgpiX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxNgpi == lenbN):
                        plt.annotate(''+str(dyNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                    
                    #plt.scatter(liste_xNgpi, liste_stringNgpi, s = lenbN/1.8, c ='black')
                    plt.bar(liste_xNgpi, liste_stringNgpi, width = 10,  align='center', color='black') # width = lenbN/250, 

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbN/80))
                plt.yticks(size=round(lenbN/80))

                #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno.svg')
                plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno.png', dpi=32) # , dpi=16
                plt.close()

            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,komplementarni vlakno, vsechny ramce - konec.

            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,komplementarni vlakno, vsechny ramce - zacatek:
            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
            # pro BIG verzi - roztazeni obr. pres celou str.
            
            if (len(omega_sum_N) == 0):
                print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
                print('\n')
            else:
                xNgpi = []
                for ixGPI_N in range(len(GPI_xcoords_N_R)):
                    xNgpi += [GPI_xcoords_N_R[ixGPI_N] + 3*(distances_AAs_N[ixGPI_N] - cleavages_N[ixGPI_N])]
                    '''
                    if (distances_AAs_N[ixGPI_N] == cleavages_N[ixGPI_N]):
                        xNgpi += [GPI_xcoords_N_R[ixGPI_N]]
                    elif (distances_AAs_N[ixGPI_N] > cleavages_N[ixGPI_N]):
                        xNgpi += [GPI_xcoords_N_R[ixGPI_N] + 3*(distances_AAs_N[ixGPI_N] - cleavages_N[ixGPI_N]) + 3]
                    else:
                        None
                    '''
                xNgpi
                yNgpi = omega_sum_N
                
                liste_xNgpi = xNgpi
                liste_yNgpi = yNgpi

                liste_xNgpi = [0] + liste_xNgpi
                liste_xNgpi = liste_xNgpi + [lenbN]

                liste_yNgpi = [0] + liste_yNgpi
                liste_yNgpi = liste_yNgpi + [0]

                liste_stringNgpi = liste_yNgpi

                fig = plt.figure()
                
                fig, ax = plt.subplots()

                plt.figure(figsize=(round(lenbN/40),round(lenbN/75)))
                plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, komplementarni vlakno', fontsize=round(lenbN/65))
                plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbN/65))
                #plt.title('Souhrn ORFu s omega-misty, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

                #ax = plt.gca()
                #fig = plt.gcf()

                for sxNgpi, dyNgpi in zip(liste_xNgpi, liste_yNgpi):
                 
                    if (sxNgpi == 0):
                        plt.annotate(''+str(dyNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    elif ((sxNgpi > 0) and (sxNgpi < lenbN)):
                        for sxNgpiX, dyNgpiX in zip(liste_xNgpi[1:-1], liste_yNgpi[1:-1]):
                            plt.annotate(''+str(dyNgpiX), xy = (sxNgpiX, dyNgpiX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                            ifNgpiX = liste_xNgpi[1:-1].index(sxNgpiX)
                            plt.annotate('     - ORF:'+str(sxNgpiX - 3*(distances_AAs_N[ifNgpiX] - cleavages_N[ifNgpiX]))+' - GPI:'+str(sxNgpiX), xy = (sxNgpiX, dyNgpiX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                        break
                    elif (sxNgpi == lenbN):
                        plt.annotate(''+str(dyNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        plt.annotate('     -  '+str(sxNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    else:
                        None
                    
                    #plt.scatter(liste_xNgpi, liste_stringNgpi, s = lenbN/1.8, c ='black')
                    plt.bar(liste_xNgpi, liste_stringNgpi, width = 10,  align='center', color='black') # width = lenbN/250, 

                plt.xticks(rotation=90)
                plt.xticks(size=round(lenbN/80))
                plt.yticks(size=round(lenbN/80))

                #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno-BIG.svg')
                plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno-BIG.png', dpi=32) # , dpi=16
                plt.close()

            # pro BIG verzi - roztazeni obr. pres celou str.
            # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
            # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,komplementarni vlakno, vsechny ramce - konec.

            index1 = open("index1.html","w")            

            from airium import Airium

            a = Airium()

            a('<!DOCTYPE html>')
            with a.html(lang='en'):
              with a.head():
                  a.meta(charset="utf-8")
                  a.title(_t="TM, CC a GPI predikce")
                  a.link(href='styly.css', rel='stylesheet')

              with a.body():
                  with a.h3(id="id17", klass='main_header'):
                    a("TM, CC a GPI predikce:")
                  with a.p():
                    a('Skript použitelný pro: standardní kód, kód archeí a rostlinných plastidů a kvasinkový alternativní jaderný kód.')
                  with a.p():
                    a('dolní limit délky vstupní sekvence je (bp):')
                    with a.b():a(dolni_limit)
                  with a.p():
                      a('horní limit délky vstupní sekvence je (bp):')
                      with a.b():a(horni_limit)
                  with a.p():
                    a('nukleotidů (délka zadané NK) je (bp):')
                    with a.b():a(lenbP)
                  with a.p():
                    a('minimální délka ORFu je (bp):')
                    with a.b():a(limitORF)
                  with a.p():
                    a('počet ORFů v zadaném vláknu je:')
                    with a.b():a(len(listORFsP))
                  with a.p():
                    a('počet ORFů v komplementárním vláknu je:')
                    with a.b():a(len(listORFsN))
                  with a.p():
                    a('celkový počet ORFů v obou vláknech je:')
                    with a.b():a(len(listORFsP) + len(listORFsN))
                  with a.p():
                    a('počet ORFů s TM sekundární strukturou (zadané vlákno): ')
                    with a.b():a(count_TM_dvoj_obrazku_P)
                  with a.p():
                    a('počet ORFů s TM sekundární strukturou (komplementární vlákno): ')
                    with a.b():a(count_TM_dvoj_obrazku_N)
                  with a.p():
                    a('počet ORFů s CC sekundární strukturou (zadané vlákno): ')
                    with a.b():a(CC_count_P)
                  with a.p():
                    a('počet ORFů s CC sekundární strukturou (komplementární vlákno): ')
                    with a.b():a(CC_count_N)
                  with a.p():
                    a('práh pro detekci CC sekundární struktury - zadané vlákno (%): ')
                    with a.b():a(CC_threshold_P*100)
                  with a.p():
                    a('práh pro detekci CC sekundární struktury - komplementární vlákno (%): ')
                    with a.b():a(CC_threshold_N*100)
                  with a.p():
                    a('počet ORFů s GPI sekundární strukturou (zadané vlákno): ')
                    with a.b():a(len(GPI_xcoords_P_R))
                  with a.p():
                    a('počet ORFů s GPI sekundární strukturou (komplementární vlákno): ')
                    with a.b():a(len(GPI_xcoords_N_R))

                  with a.p():
                    a('počet ORFů jen s TM sekundární strukturou (zadané vlákno): ')
                    with a.b():a(len(list_jenTM_P))
                  with a.p():
                    a('počet ORFů jen s TM sekundární strukturou (komplementární vlákno): ')
                    with a.b():a(len(list_jenTM_N))

                  with a.p():
                    a('počet ORFů jen s CC sekundární strukturou (zadané vlákno): ')
                    with a.b():a(len(list_jenCC_P))
                  with a.p():
                    a('počet ORFů jen s CC sekundární strukturou (komplementární vlákno): ')
                    with a.b():a(len(list_jenCC_N))

                  with a.p():
                    a('počet ORFů jen s GPI sekundární strukturou (zadané vlákno): ')
                    with a.b():a(len(list_jenGPI_P))
                  with a.p():
                    a('počet ORFů jen s GPI sekundární strukturou (komplementární vlákno): ')
                    with a.b():a(len(list_jenGPI_N))
                    
                  with a.p():
                    a('počet ORFů jen s TM a CC sekundární strukturou (zadané vlákno): ')
                    with a.b():a(len(list_TMaCC_P))
                  with a.p():
                    a('počet ORFů jen s TM a CC sekundární strukturou (komplementární vlákno): ')
                    with a.b():a(len(list_TMaCC_N))

                  with a.p():
                    a('počet ORFů jen s TM a GPI sekundární strukturou (zadané vlákno): ')
                    with a.b():a(len(list_TMaGPI_P))
                  with a.p():
                    a('počet ORFů jen s TM a GPI sekundární strukturou (komplementární vlákno): ')
                    with a.b():a(len(list_TMaGPI_N))

                  with a.p():
                    a('počet ORFů jen s CC a GPI sekundární strukturou (zadané vlákno): ')
                    with a.b():a(len(list_GPIaCC_P))
                  with a.p():
                    a('počet ORFů jen s CC a GPI sekundární strukturou (komplementární vlákno): ')
                    with a.b():a(len(list_GPIaCC_N))

                  with a.p():
                    a('počet ORFů s TM, CC a GPI sekundární strukturou současně (zadané vlákno): ')
                    with a.b():a(len(list_TMaCCaGPI_P))
                  with a.p():
                    a('počet ORFů s TM, CC a GPI sekundární strukturou současně (komplementární vlákno): ')
                    with a.b():a(len(list_TMaCCaGPI_N))

                  #with a.p(align='center'):a.a(href='#(konec stránky)', _t='(přejít na konec stránky)')

                  a.hr(width="100%", size="5", color="yellow", align="center")
                  a.hr(width="100%", size="5", color="yellow", align="center")
                  with a.h2(id="id8", klass='main_header'):a('ORFy - zadané vlákno DNA:')
                  with a.table(align="center"):
                      with a.tr():
                          a.td(a.img(src="figures/NA-ORFs-Pframe1-given_string-2D.png", alt="", width="33%"))
                          a.td(a.img(src="figures/NA-ORFs-Pframe2-given_string-2D.png", alt="", width="33%"))
                          a.td(a.img(src="figures/NA-ORFs-Pframe3-given_string-2D.png", alt="", width="33%"))
                      with a.tr():
                          a.td(a.img(src="figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno.png", alt="", width="33%"))
                          a.td(a.img(src="figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno.png", alt="", width="33%"))
                          a.td(a.img(src="figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno.png", alt="", width="33%"))

                  with a.table(align="center"):
                      with a.tr():
                          a.td(a.img(src="figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG.png", alt="", width="100%"))
                      with a.tr():
                          a.td(a.img(src="figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG.png", alt="", width="100%"))
                      with a.tr():
                          a.td(a.img(src="figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG.png", alt="", width="100%"))

                          
                  with a.table(align="center", cellpadding='5', border='2', bordercolor='black'): #
                       with a.tr():
                          with a.th(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s TM strukturou - zadané vlákno')
                          with a.th(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s CC strukturou - zadané vlákno')
                          with a.th(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s GPI strukturou - zadané vlákno')
                       with a.tr(): #
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('1. čtecí rámec:')
                                  for itmP1all in listTMsP_pozice:
                                          if ((itmP1all[0])%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(itmP1all), _t='ORF-pozice: '+str(itmP1all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('1. čtecí rámec:')
                                  for iccP1all in CC_all_P_coords:
                                          if ((iccP1all[0])%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(iccP1all), _t='ORF-pozice: '+str(iccP1all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('1. čtecí rámec:')
                                  for igpiP1all in GPI_XYcoords_P_R:
                                          if ((igpiP1all[0])%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(igpiP1all), _t='ORF-pozice: '+str(igpiP1all)+'\n')
                                          else:
                                              None
                       with a.tr(): ##
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('2. čtecí rámec:')
                                  for itmP2all in listTMsP_pozice:
                                          if ((itmP2all[0]+2)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(itmP2all), _t='ORF-pozice: '+str(itmP2all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('2. čtecí rámec:')
                                  for iccP2all in CC_all_P_coords:
                                          if ((iccP2all[0]+2)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(iccP2all), _t='ORF-pozice: '+str(iccP2all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('2. čtecí rámec:')
                                  for igpiP2all in GPI_XYcoords_P_R:
                                          if ((igpiP2all[0]+2)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(igpiP2all), _t='ORF-pozice: '+str(igpiP2all)+'\n')
                                          else:
                                              None
                       with a.tr(): ###
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('3. čtecí rámec:')
                                  for itmP3all in listTMsP_pozice:
                                          if ((itmP3all[0]+1)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(itmP3all), _t='ORF-pozice: '+str(itmP3all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('3. čtecí rámec:')
                                  for iccP3all in CC_all_P_coords:
                                          if ((iccP3all[0]+1)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(iccP3all), _t='ORF-pozice: '+str(iccP3all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('3. čtecí rámec:')
                                  for igpiP3all in GPI_XYcoords_P_R:
                                          if ((igpiP3all[0]+1)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(igpiP3all), _t='ORF-pozice: '+str(igpiP3all)+'\n')
                                          else:
                                              None
                  with a.table(align="center", cellpadding='5', border='2', bordercolor='black'): # width='99%'
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM strukturou - zadané vlákno - 1. čtecí rámec:')
                                  for itmP1 in list_jenTM_P_01:
                                      if ((itmP1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmP1), _t='ORF-pozice: '+str(itmP1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM strukturou - zadané vlákno - 2. čtecí rámec:')
                                  for itmP2 in list_jenTM_P_01:
                                      if ((itmP2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmP2), _t='ORF-pozice: '+str(itmP2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM strukturou - zadané vlákno - 3. čtecí rámec:')
                                  for itmP3 in list_jenTM_P_01:
                                      if ((itmP3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmP3), _t='ORF-pozice: '+str(itmP3)+'\n')
                                      else:
                                          None # # #
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC strukturou - zadané vlákno - 1. čtecí rámec:')
                                  for iccP1 in list_jenCC_P_01:
                                      if ((iccP1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccP1), _t='ORF-pozice: '+str(iccP1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC strukturou - zadané vlákno - 2. čtecí rámec:')
                                  for iccP2 in list_jenCC_P_01:
                                      if ((iccP2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccP2), _t='ORF-pozice: '+str(iccP2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC strukturou - zadané vlákno - 3. čtecí rámec:')
                                  for iccP3 in list_jenCC_P_01:
                                      if ((iccP3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccP3), _t='ORF-pozice: '+str(iccP3)+'\n')
                                      else:
                                          None # # # 
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s GPI strukturou - zadané vlákno - 1. čtecí rámec:')
                                  for igpiP1 in list_jenGPI_P_01:
                                      if ((igpiP1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiP1), _t='ORF-pozice: '+str(igpiP1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s GPI strukturou - zadané vlákno - 2. čtecí rámec:')
                                  for igpiP2 in list_jenGPI_P_01:
                                      if ((igpiP2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiP2), _t='ORF-pozice: '+str(igpiP2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s GPI strukturou - zadané vlákno - 3. čtecí rámec:')
                                  for igpiP3 in list_jenGPI_P_01:
                                      if ((igpiP3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiP3), _t='ORF-pozice: '+str(igpiP3)+'\n')
                                      else:
                                          None
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a CC strukturou - zadané vlákno - 1. čtecí rámec:')
                                  for itmccP1 in list_TMaCC_P_01:
                                      if ((itmccP1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccP1), _t='ORF-pozice: '+str(itmccP1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a CC strukturou - zadané vlákno - 2. čtecí rámec:')
                                  for itmccP2 in list_TMaCC_P_01:
                                      if ((itmccP2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccP2), _t='ORF-pozice: '+str(itmccP2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a CC strukturou - zadané vlákno - 3. čtecí rámec:')
                                  for itmccP3 in list_TMaCC_P_01:
                                      if ((itmccP3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccP3), _t='ORF-pozice: '+str(itmccP3)+'\n')
                                      else:
                                          None # # #
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a GPI strukturou - zadané vlákno - 1. čtecí rámec:')
                                  for itmgpiP1 in list_TMaGPI_P_01:
                                      if ((itmgpiP1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmgpiP1), _t='ORF-pozice: '+str(itmgpiP1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a GPI strukturou - zadané vlákno - 2. čtecí rámec:')
                                  for itmgpiP2 in list_TMaGPI_P_01:
                                      if ((itmgpiP2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmgpiP2), _t='ORF-pozice: '+str(itmgpiP2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a GPI strukturou - zadané vlákno - 3. čtecí rámec:')
                                  for itmgpiP3 in list_TMaGPI_P_01:
                                      if ((itmgpiP3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmgpiP3), _t='ORF-pozice: '+str(itmgpiP3)+'\n')
                                      else:
                                          None # # #
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC a GPI strukturou - zadané vlákno - 1. čtecí rámec:')
                                  for iccgpiP1 in list_GPIaCC_P_01:
                                      if ((iccgpiP1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccgpiP1), _t='ORF-pozice: '+str(iccgpiP1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC a GPI strukturou - zadané vlákno - 2. čtecí rámec:')
                                  for iccgpiP2 in list_GPIaCC_P_01:
                                      if ((iccgpiP2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccgpiP2), _t='ORF-pozice: '+str(iccgpiP2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC a GPI strukturou - zadané vlákno - 3. čtecí rámec:')
                                  for iccgpiP3 in list_GPIaCC_P_01:
                                      if ((iccgpiP3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccgpiP3), _t='ORF-pozice: '+str(iccgpiP3)+'\n')
                                      else:
                                          None # # #
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s TM, CC a GPI strukturou - zadané vlákno - 1. čtecí rámec:')
                                  for itmccgpiP1 in list_TMaCCaGPI_P_01:
                                      if ((itmccgpiP1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccgpiP1), _t='ORF-pozice: '+str(itmccgpiP1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s TM, CC a GPI strukturou - zadané vlákno - 2. čtecí rámec:')
                                  for itmccgpiP2 in list_TMaCCaGPI_P_01:
                                      if ((itmccgpiP2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccgpiP2), _t='ORF-pozice: '+str(itmccgpiP2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s TM, CC a GPI strukturou - zadané vlákno - 3. čtecí rámec:')
                                  for itmccgpiP3 in list_TMaCCaGPI_P_01:
                                      if ((itmccgpiP3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccgpiP3), _t='ORF-pozice: '+str(itmccgpiP3)+'\n')
                                      else:
                                          None # # #
                  a.hr(width="100%", size="5", color="yellow", align="center")
                  with a.h2(id="id9", klass='main_header'):a('ORFy - komplementární vlákno DNA:')
                  with a.table(align="center"):
                      with a.tr():
                          a.td(a.img(src="figures/NA-ORFs-Nframe1-complementary_string-2D.png", alt="", width="33%"))
                          a.td(a.img(src="figures/NA-ORFs-Nframe2-complementary_string-2D.png", alt="", width="33%"))
                          a.td(a.img(src="figures/NA-ORFs-Nframe3-complementary_string-2D.png", alt="", width="33%"))
                      with a.tr():
                          a.td(a.img(src="figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno.png", alt="", width="33%"))
                          a.td(a.img(src="figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno.png", alt="", width="33%"))
                          a.td(a.img(src="figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno.png", alt="", width="33%"))

                  with a.table(align="center"):
                      with a.tr():
                          a.td(a.img(src="figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno-BIG.png", alt="", width="100%"))
                      with a.tr():
                          a.td(a.img(src="figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno-BIG.png", alt="", width="100%"))
                      with a.tr():
                          a.td(a.img(src="figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno-BIG.png", alt="", width="100%"))
                          
                  with a.table(align="center", cellpadding='5', border='2', bordercolor='black'): #
                       with a.tr():
                          with a.th(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s TM strukturou - komplementární vlákno')
                          with a.th(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s CC strukturou - komplementární vlákno')
                          with a.th(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s GPI strukturou - komplementární vláno')
                       with a.tr(): #
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('1. čtecí rámec:')
                                  for itmN1all in listTMsN_pozice:
                                          if ((itmN1all[0])%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(itmN1all), _t='ORF-pozice: '+str(itmN1all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('1. čtecí rámec:')
                                  for iccN1all in CC_all_N_coords:
                                          if ((iccN1all[0])%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(iccN1all), _t='ORF-pozice: '+str(iccN1all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('1. čtecí rámec:')
                                  for igpiN1all in GPI_XYcoords_N_R:
                                          if ((igpiN1all[0])%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(igpiN1all), _t='ORF-pozice: '+str(igpiN1all)+'\n')
                                          else:
                                              None
                       with a.tr(): ##
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('2. čtecí rámec:')
                                  for itmN2all in listTMsN_pozice:
                                          if ((itmN2all[0]+2)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(itmN2all), _t='ORF-pozice: '+str(itmN2all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('2. čtecí rámec:')
                                  for iccN2all in CC_all_N_coords:
                                          if ((iccN2all[0]+2)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(iccN2all), _t='ORF-pozice: '+str(iccN2all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('2. čtecí rámec:')
                                  for igpiN2all in GPI_XYcoords_N_R:
                                          if ((igpiN2all[0]+2)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(igpiN2all), _t='ORF-pozice: '+str(igpiN2all)+'\n')
                                          else:
                                              None
                       with a.tr(): ###
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('3. čtecí rámec:')
                                  for itmN3all in listTMsN_pozice:
                                          if ((itmN3all[0]+1)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(itmN3all), _t='ORF-pozice: '+str(itmN3all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('3. čtecí rámec:')
                                  for iccN3all in CC_all_N_coords:
                                          if ((iccN3all[0]+1)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(iccN3all), _t='ORF-pozice: '+str(iccN3all)+'\n')
                                          else:
                                              None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('3. čtecí rámec:')
                                  for igpiN3all in GPI_XYcoords_N_R:
                                          if ((igpiN3all[0]+1)%3 == 0)==True:
                                              with a.p(align='center'):a.a(href='#'+str(igpiN3all), _t='ORF-pozice: '+str(igpiN3all)+'\n')
                                          else:
                                              None

                  with a.table(align="center", cellpadding='5', border='2', bordercolor='black'): # width='99%'
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM strukturou - komplementární vlákno - 1. čtecí rámec:')
                                  for itmN1 in list_jenTM_N_01:
                                      if ((itmN1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmN1), _t='ORF-pozice: '+str(itmN1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM strukturou - komplementární vlákno - 2. čtecí rámec:')
                                  for itmN2 in list_jenTM_N_01:
                                      if ((itmN2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmN2), _t='ORF-pozice: '+str(itmN2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM strukturou - komplementární vlákno - 3. čtecí rámec:')
                                  for itmN3 in list_jenTM_N_01:
                                      if ((itmN3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmN3), _t='ORF-pozice: '+str(itmN3)+'\n')
                                      else:
                                          None # # #
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC strukturou - komplementární vlákno - 1. čtecí rámec:')
                                  for iccN1 in list_jenCC_N_01:
                                      if ((iccN1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccN1), _t='ORF-pozice: '+str(iccN1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC strukturou - komplementární vlákno - 2. čtecí rámec:')
                                  for iccN2 in list_jenCC_N_01:
                                      if ((iccN2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccN2), _t='ORF-pozice: '+str(iccN2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC strukturou - komplementární vlákno - 3. čtecí rámec:')
                                  for iccN3 in list_jenCC_N_01:
                                      if ((iccN3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccN3), _t='ORF-pozice: '+str(iccN3)+'\n')
                                      else:
                                          None # # # 
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s GPI strukturou - komplementární vlákno - 1. čtecí rámec:')
                                  for igpiN1 in list_jenGPI_N_01:
                                      if ((igpiN1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiN1), _t='ORF-pozice: '+str(igpiN1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s GPI strukturou - komplementární vlákno - 2. čtecí rámec:')
                                  for igpiN2 in list_jenGPI_N_01:
                                      if ((igpiN2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiN2), _t='ORF-pozice: '+str(igpiN2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s GPI strukturou - komplementární vlákno - 3. čtecí rámec:')
                                  for igpiN3 in list_jenGPI_N_01:
                                      if ((igpiN3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiN3), _t='ORF-pozice: '+str(igpiN3)+'\n')
                                      else:
                                          None
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a CC strukturou - komplementární vlákno - 1. čtecí rámec:')
                                  for itmccN1 in list_TMaCC_N_01:
                                      if ((itmccN1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccN1), _t='ORF-pozice: '+str(itmccN1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a CC strukturou - komplementární vlákno - 2. čtecí rámec:')
                                  for itmccN2 in list_TMaCC_N_01:
                                      if ((itmccN2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccN2), _t='ORF-pozice: '+str(itmccN2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a CC strukturou - komplementární vlákno - 3. čtecí rámec:')
                                  for itmccN3 in list_TMaCC_N_01:
                                      if ((itmccN3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccN3), _t='ORF-pozice: '+str(itmccN3)+'\n')
                                      else:
                                          None # # #
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a GPI strukturou - komplementární vlákno - 1. čtecí rámec:')
                                  for itmgpiN1 in list_TMaGPI_N_01:
                                      if ((itmgpiN1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmgpiN1), _t='ORF-pozice: '+str(itmgpiN1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a GPI strukturou - komplementární vlákno - 2. čtecí rámec:')
                                  for itmgpiN2 in list_TMaGPI_N_01:
                                      if ((itmgpiN2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmgpiN2), _t='ORF-pozice: '+str(itmgpiN2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s TM a GPI strukturou - komplementární vlákno - 3. čtecí rámec:')
                                  for itmgpiN3 in list_TMaGPI_N_01:
                                      if ((itmgpiN3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmgpiN3), _t='ORF-pozice: '+str(itmgpiN3)+'\n')
                                      else:
                                          None # # #
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC a GPI strukturou - komplementární vlákno - 1. čtecí rámec:')
                                  for iccgpiN1 in list_GPIaCC_N_01:
                                      if ((iccgpiN1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccgpiN1), _t='ORF-pozice: '+str(iccgpiN1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC a GPI strukturou - komplementární vlákno - 2. čtecí rámec:')
                                  for iccgpiN2 in list_GPIaCC_N_01:
                                      if ((iccgpiN2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccgpiN2), _t='ORF-pozice: '+str(iccgpiN2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy jen s CC a GPI strukturou - komplementární vlákno - 3. čtecí rámec:')
                                  for iccgpiN3 in list_GPIaCC_N_01:
                                      if ((iccgpiN3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccgpiN3), _t='ORF-pozice: '+str(iccgpiN3)+'\n')
                                      else:
                                          None # # #
                       with a.tr():
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s TM, CC a GPI strukturou - komplementární vlákno - 1. čtecí rámec:')
                                  for itmccgpiN1 in list_TMaCCaGPI_N_01:
                                      if ((itmccgpiN1[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccgpiN1), _t='ORF-pozice: '+str(itmccgpiN1)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s TM, CC a GPI strukturou - komplementární vlákno - 2. čtecí rámec:')
                                  for itmccgpiN2 in list_TMaCCaGPI_N_01:
                                      if ((itmccgpiN2[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccgpiN2), _t='ORF-pozice: '+str(itmccgpiN2)+'\n')
                                      else:
                                          None
                          with a.td(width='33%'):
                              with a.p(align='center'):
                                  with a.b():a('ORFy s TM, CC a GPI strukturou - komplementární vlákno - 3. čtecí rámec:')
                                  for itmccgpiN3 in list_TMaCCaGPI_N_01:
                                      if ((itmccgpiN3[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmccgpiN3), _t='ORF-pozice: '+str(itmccgpiN3)+'\n')
                                      else:
                                          None # # #

                  a.hr(width="100%", size="5", color="yellow", align="center")
                  with a.h2(id="id10", klass='main_header'):
                    a("jednotlivé ORFy (zadané vlákno):")
                    if (len(listORFsP) == 0):
                        with a.p():a('Žádný ORF a peptid v zadaném vláknu.')
                    else:
                        for iximgP1 in listORFsP_pozice:
                              if ((iximgP1[0])%3 == 0)==True:
                                if (os.path.isfile('figures/given/NA string - fig No. ' + str(iximgP1) + ' provided_string.png'))==True:
                                      with a.p():a('ORF, zadané vlákno, 1. čtecí rámec - pozice: ')
                                      with a.b(id=str(iximgP1)):a(str(iximgP1) + ' z: ' + str(lenbP) + ':')
                                      a.img(src="figures/given/NA string - fig No. " + str(iximgP1) + " provided_string.png", width = "100%", alt="")
                                      with a.table(align="center", cellpadding='2', border='2', bordercolor='black'):
                                          with a.tr():
                                              with a.td():
                                                  if (os.path.isfile('figures/given/Fig No. ' + str(iximgP1) + ' provided_string.png'))==True:
                                                      a.img(src="figures/given/Fig No. " + str(iximgP1) + " provided_string.png", alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td():
                                                  if (os.path.isfile('CC_plots/plotsP-realCC/' + str(iximgP1[0]) + '.png'))==True:
                                                      a.img(src='CC_plots/plotsP-realCC/' + str(iximgP1[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td(align='center'): # bgcolor='green', ='white'
                                                  if (os.path.isfile('figures/GPI/P' + str(iximgP1) + '.png'))==True:
                                                      a.img(src='figures/GPI/P' + str(iximgP1) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None

                                      with a.p(align='center'):
                                        a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                        #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>ORF_'+str(iximgP1[0])+'-'+str(iximgP1[1]))
                                              a.br()
                                              a(str(listORFsP[listORFsP_pozice.index(iximgP1)]))
                                          #with a.td(width='100%'):a([listORFsP[listORFsP_pozice.index(iximgP1)][iSekP1:iSekP1+90] for iSekP1 in range(0,len(listORFsP[listORFsP_pozice.index(iximgP1)]),90)])

                                      with a.p(align='center'):
                                        a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>trORF_'+str(iximgP1[0])+'-'+str(iximgP1[1]))
                                              a.br()
                                              a(str(PepsP_c00[listORFsP_pozice.index(iximgP1)]))

                                      a.hr(width="100%", size="5", color="yellow", align="center")
                                    
                                #else:
                                    #with a.p():a('Neni ORF a peptid v 1. ctecim ramci v zadanem vlaknu.')
                        with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                        a.hr(width="100%", size="5", color="yellow", align="center")

                    if (len(listORFsP) == 0):
                        with a.p():a('Zadny ORF a peptid v zadanem vlaknu.')
                    else:
                        for iximgP2 in listORFsP_pozice:
                              if ((iximgP2[0]+2)%3 == 0)==True:
                                if (os.path.isfile('figures/given/NA string - fig No. ' + str(iximgP2) + ' provided_string.png'))==True:
                                      with a.p():a('ORF, zadané vlákno, 2. čtecí rámec - pozice: ')
                                      with a.b(id=str(iximgP2)):a(str(iximgP2) + ' z: ' + str(lenbP) + ':')
                                      a.img(src="figures/given/NA string - fig No. " + str(iximgP2) + " provided_string.png", width = "100%", alt="")
                                      with a.table(align="center", cellpadding='2', border='2', bordercolor='black'):
                                          with a.tr():
                                              with a.td():
                                                  if (os.path.isfile('figures/given/Fig No. ' + str(iximgP2) + ' provided_string.png'))==True:
                                                      a.img(src="figures/given/Fig No. " + str(iximgP2) + " provided_string.png", alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td():
                                                  if (os.path.isfile('CC_plots/plotsP-realCC/' + str(iximgP2[0]) + '.png'))==True:
                                                      a.img(src='CC_plots/plotsP-realCC/' + str(iximgP2[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td(align='center'): # bgcolor='green'
                                                  if (os.path.isfile('figures/GPI/P' + str(iximgP2) + '.png'))==True:
                                                      a.img(src='figures/GPI/P' + str(iximgP2) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                                      
                                      with a.p(align='center'):
                                        a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                        #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>ORF_'+str(iximgP2[0])+'-'+str(iximgP2[1]))
                                              a.br()
                                              a(str(listORFsP[listORFsP_pozice.index(iximgP2)]))
                                          #with a.td(width='100%'):a([listORFsP[listORFsP_pozice.index(iximgP2)][iSekP2:iSekP2+90] for iSekP2 in range(0,len(listORFsP[listORFsP_pozice.index(iximgP2)]),90)])

                                      with a.p(align='center'):
                                        a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>trORF_'+str(iximgP2[0])+'-'+str(iximgP2[1]))
                                              a.br()
                                              a(str(PepsP_c00[listORFsP_pozice.index(iximgP2)]))

                                      a.hr(width="100%", size="5", color="yellow", align="center") 
                                #else:
                                    #with a.p():a('Neni ORF a peptid ve 2. ctecim ramci v zadanem vlaknu.')
                        with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                        a.hr(width="100%", size="5", color="yellow", align="center")

                    if (len(listORFsP) == 0):
                        with a.p():a('Zadny ORF a peptid v zadanem vlaknu.')
                    else:
                        for iximgP3 in listORFsP_pozice:
                              if ((iximgP3[0]+1)%3 == 0)==True:
                                if (os.path.isfile('figures/given/NA string - fig No. ' + str(iximgP3) + ' provided_string.png'))==True:
                                      with a.p():a('ORF, zadané vlákno, 3. čtecí rámec - pozice: ')
                                      with a.b(id=str(iximgP3)):a(str(iximgP3) + ' z: ' + str(lenbP) + ':')
                                      a.img(src="figures/given/NA string - fig No. " + str(iximgP3) + " provided_string.png", width = "100%", alt="")
                                      with a.table(align="center", cellpadding='2', border='2', bordercolor='black'):
                                          with a.tr():
                                              with a.td():
                                                  if (os.path.isfile('figures/given/Fig No. ' + str(iximgP3) + ' provided_string.png'))==True:
                                                      a.img(src="figures/given/Fig No. " + str(iximgP3) + " provided_string.png", alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td():
                                                  if (os.path.isfile('CC_plots/plotsP-realCC/' + str(iximgP3[0]) + '.png'))==True:
                                                      a.img(src='CC_plots/plotsP-realCC/' + str(iximgP3[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td(align='center'): # bgcolor='green'
                                                  if (os.path.isfile('figures/GPI/P' + str(iximgP3) + '.png'))==True:
                                                      a.img(src='figures/GPI/P' + str(iximgP3) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None

                                      with a.p(align='center'):
                                        a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                        #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>ORF_'+str(iximgP3[0])+'-'+str(iximgP3[1]))
                                              a.br()
                                              a(str(listORFsP[listORFsP_pozice.index(iximgP3)]))
                                          #with a.td(width='100%'):a([listORFsP[listORFsP_pozice.index(iximgP3)][iSekP3:iSekP3+90] for iSekP3 in range(0,len(listORFsP[listORFsP_pozice.index(iximgP3)]),90)])


                                      with a.p(align='center'):
                                        a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>trORF_'+str(iximgP3[0])+'-'+str(iximgP3[1]))
                                              a.br()
                                              a(str(PepsP_c00[listORFsP_pozice.index(iximgP3)]))

                                                      
                                      a.hr(width="100%", size="5", color="yellow", align="center")
                                #else:
                                    #with a.p():a('Neni ORF a peptid ve 3. ctecim ramci v zadanem vlaknu.')
                        with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                        a.hr(width="100%", size="5", color="yellow", align="center")

                  a.hr(width="100%", size="5", color="red", align="center")
                  a.hr(width="100%", size="5", color="yellow", align="center")
                  with a.h2(id="id11", klass='main_header', align='center'):
                    a("jednotlivé ORFy (komplementární vlákno):")
                    if (len(listORFsN) == 0):
                        with a.p():a('Žádný ORF a peptid v komplementárním vláknu.')
                    else:
                        for iximgN1 in listORFsN_pozice:
                              if ((iximgN1[0])%3 == 0)==True:
                                if (os.path.isfile('figures/complementary/NA string - fig No. ' + str(iximgN1) + ' complementary_string.png'))==True:
                                      with a.p():a('ORF, komplementární vlákno, 1. čtecí rámec - pozice: ')
                                      with a.b(id=str(iximgN1)):a(str(iximgN1) + ' z: ' + str(lenbN) + ':')
                                      a.img(src="figures/complementary/NA string - fig No. " + str(iximgN1) + " complementary_string.png", width = "100%", alt="")
                                      with a.table(align="center", cellpadding='2', border='2', bordercolor='black'):
                                          with a.tr():
                                              with a.td():
                                                  if (os.path.isfile('figures/complementary/Fig No. ' + str(iximgN1) + ' complementary_string.png'))==True:
                                                      a.img(src="figures/complementary/Fig No. " + str(iximgN1) + " complementary_string.png", alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td():
                                                  if (os.path.isfile('CC_plots/plotsN-realCC/' + str(iximgN1[0]) + '.png'))==True:
                                                      a.img(src='CC_plots/plotsN-realCC/' + str(iximgN1[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td(align='center'): # bgcolor='green' , bgcolor='white'
                                                  if (os.path.isfile('figures/GPI/N/' + str(iximgN1) + '.png'))==True:
                                                      a.img(src='figures/GPI/N/' + str(iximgN1) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None

                                      with a.p(align='center'):
                                        a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                        #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>ORF_'+str(iximgN1[0])+'-'+str(iximgN1[1]))
                                              a.br()
                                              a(str(listORFsN[listORFsN_pozice.index(iximgN1)]))
                                          #with a.td(width='100%'):a([listORFsN[listORFsN_pozice.index(iximgN1)][iSekN1:iSekN1+90] for iSekN1 in range(0,len(listORFsN[listORFsN_pozice.index(iximgN1)]),90)])

                                      with a.p(align='center'):
                                        a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>trORF_'+str(iximgN1[0])+'-'+str(iximgN1[1]))
                                              a.br()
                                              a(str(PepsN_c00[listORFsN_pozice.index(iximgN1)]))
                                                      
                                      a.hr(width="100%", size="5", color="yellow", align="center")
                                #else:
                                    #with a.p():a('Neni ORF a peptid v 1. ctecim ramci v komplementarnim vlaknu.')
                        with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                        a.hr(width="100%", size="5", color="yellow", align="center")

                    if (len(listORFsN) == 0):
                        with a.p():a('Žádný ORF a peptid v komplementárním vláknu.')
                    else:
                        for iximgN2 in listORFsN_pozice:
                              if ((iximgN2[0]+2)%3 == 0)==True:
                                if (os.path.isfile('figures/complementary/NA string - fig No. ' + str(iximgN2) + ' complementary_string.png'))==True:
                                      with a.p():a('ORF, komplementární vlákno, 2. čtecí rámec - pozice: ')
                                      with a.b(id=str(iximgN2)):a(str(iximgN2) + ' z: ' + str(lenbN) + ':')
                                      a.img(src="figures/complementary/NA string - fig No. " + str(iximgN2) + " complementary_string.png", width = "100%", alt="")
                                      with a.table(align="center", cellpadding='2', border='2', bordercolor='black'):
                                          with a.tr():
                                              with a.td():
                                                  if (os.path.isfile('figures/complementary/Fig No. ' + str(iximgN2) + ' complementary_string.png'))==True:
                                                      a.img(src="figures/complementary/Fig No. " + str(iximgN2) + " complementary_string.png", alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td():
                                                  if (os.path.isfile('CC_plots/plotsN-realCC/' + str(iximgN2[0]) + '.png'))==True:
                                                      a.img(src='CC_plots/plotsN-realCC/' + str(iximgN2[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td(align='center'): # bgcolor='green' , bgcolor='white'
                                                  if (os.path.isfile('figures/GPI/N/' + str(iximgN2) + '.png'))==True:
                                                      a.img(src='figures/GPI/N/' + str(iximgN2) + '.png', alt="", width="550", height="450", align="absmiddle")
                                                  else:
                                                      None

                                      with a.p(align='center'):
                                        a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                        #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>ORF_'+str(iximgN2[0])+'-'+str(iximgN2[1]))
                                              a.br()
                                              a(str(listORFsN[listORFsN_pozice.index(iximgN2)]))
                                          #with a.td(width='100%'):a([listORFsN[listORFsN_pozice.index(iximgN2)][iSekN2:iSekN2+90] for iSekN2 in range(0,len(listORFsN[listORFsN_pozice.index(iximgN2)]),90)])

                                      with a.p(align='center'):
                                        a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>trORF_'+str(iximgN2[0])+'-'+str(iximgN2[1]))
                                              a.br()
                                              a(str(PepsN_c00[listORFsN_pozice.index(iximgN2)]))

                                      a.hr(width="100%", size="5", color="yellow", align="center")
                                #else:
                                    #with a.p():a('Neni ORF a peptid ve 2. ctecim ramci v komplementarnim vlaknu.')
                        with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                        a.hr(width="100%", size="5", color="yellow", align="center")

                    if (len(listORFsN) == 0):
                        with a.p():a('Žádný ORF a peptid v komplementárním vláknu.')
                    else:
                        for iximgN3 in listORFsN_pozice:
                              if ((iximgN3[0]+1)%3 == 0)==True:
                                if (os.path.isfile('figures/complementary/NA string - fig No. ' + str(iximgN3) + ' complementary_string.png'))==True:
                                      with a.p():a('ORF, komplementární vlákno, 3. čtecí rámec - pozice: ')
                                      with a.b(id=str(iximgN3)):a(str(iximgN3) + ' z: ' + str(lenbN) + ':')
                                      a.img(src="figures/complementary/NA string - fig No. " + str(iximgN3) + " complementary_string.png", width = "100%", alt="")
                                      with a.table(align="center", cellpadding='2', border='2', bordercolor='black'):
                                          with a.tr(width = "100%"):
                                              with a.td():
                                                  if (os.path.isfile('figures/complementary/Fig No. ' + str(iximgN3) + ' complementary_string.png'))==True:
                                                      a.img(src="figures/complementary/Fig No. " + str(iximgN3) + " complementary_string.png", alt="", width = "550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td():
                                                  if (os.path.isfile('CC_plots/plotsN-realCC/' + str(iximgN3[0]) + '.png'))==True:
                                                      a.img(src='CC_plots/plotsN-realCC/' + str(iximgN3[0]) + '.png', alt="", width = "550", height="450", align="absmiddle")
                                                  else:
                                                      None
                                              with a.td(align='center'): # bgcolor='white'
                                                  if (os.path.isfile('figures/GPI/N/' + str(iximgN3) + '.png'))==True:
                                                      a.img(src='figures/GPI/N/' + str(iximgN3) + '.png', alt="", width = "550", height="450", align="absmiddle")
                                                  else:
                                                      None

                                      with a.p(align='center'):
                                        a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                        #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>ORF_'+str(iximgN3[0])+'-'+str(iximgN3[1]))
                                              a.br()
                                              a(str(listORFsN[listORFsN_pozice.index(iximgN3)]))
                                              #with a.pre():a(str(listORFsN[listORFsN_pozice.index(iximgN3)]))
                                          #with a.td(width='100%'):a([listORFsN[listORFsN_pozice.index(iximgN3)][iSekN3:iSekN3+90] for iSekN3 in range(0,len(listORFsN[listORFsN_pozice.index(iximgN3)]),90)])

                                      with a.p(align='center'):
                                        a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                      with a.table(align="center", cellpadding='5', border='2', bordercolor='black', width='100%'):
                                        with a.tr(width='100%'):
                                          with a.td(width='100%'):
                                              a('>trORF_'+str(iximgN3[0])+'-'+str(iximgN3[1]))
                                              a.br()
                                              a(str(PepsN_c00[listORFsN_pozice.index(iximgN3)]))
                                              #with a.pre():a(str(PepsN_c00[listORFsN_pozice.index(iximgN3)]))
                    
                                      a.hr(width="100%", size="5", color="yellow", align="center")
                                #else:
                                    #with a.p():a('Neni ORF a peptid ve 3. ctecim ramci v komplementarnim vlaknu.')
                        with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                        a.hr(width="100%", size="5", color="yellow", align="center")
                        
                  #with a.p(align="center"):
                  #a('(konec stránky)')

                  a.hr(width="100%", size="5", color="red", align="center")
                  a.hr(width="100%", size="5", color="red", align="center")

            index1_str_a = str(a)

            index1.write(index1_str_a)
            index1.close()
            
            #break

            print('\n')
            
            print('Vysledek byl ulozen do HTML souboru: "index1.html".\n')
            print('V dialogovem okne je doporuceno pro setrvani v Python 3.x.x Shellu zmacknout tlacitko: "CANCEL".\n')

            elapsed = time.time() - t
            print('Elapsed time (seconds): ',elapsed)

            quit()

        elif (proces123 == 2):

            basesRP_2 = basesRP

            for iBP in basesRP:
                if ((iBP !=str("A"))&(iBP !=str("G"))&(iBP !=str("T"))&(iBP !=str("C"))):
                    basesRP_2 = basesRP.replace(iBP,'')
                else:
                    continue
                    #basesRP_2 = basesRP.replace(iBP,iBP)

            basesP = basesRP_2

            lenbP = len(basesP)

            #print("Nukleotidu (delka NA) je (po redukci o pripadny znak noveho radku) = ",lenbP)

            if (lenbP < dolni_limit):
                print('Delka sekvence je pod dolnim povolenym limitem ',dolni_limit,' paru bazi.')
            elif (lenbP > horni_limit):
                print('Delka sekvence je nad hornim povolenym limitem, ',horni_limit,' paru bazi.')
            else:
                print('\n')

                #print('Vypis zadane sekvence:')
                #print(bases)
                #print('\n')

                #lenb = len(bases)
                #print("Nukleotidu (delka NA) je = ",lenb)
                #print('\n')

                import re

                # vyhledani ORFu (v zadanem vlakne)
                # def.2 podle :
                '''
                Def. 2: ORF je usek DNA mezi dvemi STOP kodony a pocet jeho nukleotidu je deliteny 3 beze zbytku
                '''
                stop_codon_findP = [mPstop.start() for mPstop in re.finditer('(TAA|TAG|TGA)', basesP)]
                #print('Stop kodon v zadanem vlakne je na techto pozicich: \n',stop_codon_findP)
                #print('Stop kodon je tolikrat:',len(stop_codon_findP))
                ## najde opravdu, kde zacina dana sekvence, tedy kde je "T"
                ## vyhodi cislo, kde je prave to pismeno, napr. je-li "T" 6. zleva pak: "5"

                limitORFP = limitORF # DULEZITY nastavitelny parametr

                stop_kodon_1_ramec_P = []
                stop_kodon_2_ramec_P = []
                stop_kodon_3_ramec_P = []

                listORFsP = []
                
                listORFsP_1 = []
                listORFsP_2 = []
                listORFsP_3 = []

                listORFsP_1_pozice = []
                listORFsP_2_pozice = []
                listORFsP_3_pozice = []

                for ixyzP in stop_codon_findP:
                    if ((ixyzP)%3 == 0):
                        stop_kodon_1_ramec_P += [ixyzP]
                    elif ((ixyzP+2)%3 == 0):
                        stop_kodon_2_ramec_P += [ixyzP]
                    elif ((ixyzP+1)%3 == 0):
                        stop_kodon_3_ramec_P += [ixyzP]

                #print('Zadane vlakno - stop_kodon_1_ramec: \n',stop_kodon_1_ramec_P)
                #print('Zadane vlakno - stop_kodon_2_ramec: \n',stop_kodon_2_ramec_P)
                #print('Zadane vlakno - stop_kodon_3_ramec: \n',stop_kodon_3_ramec_P)


                if (len(stop_kodon_1_ramec_P) < 2):
                    print('Nedostecny pocet STOP kodonu v 1. ramci v zadanem vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
                else:
                    for ix1P in range(0,len(stop_kodon_1_ramec_P)-1,1):
                        if ((stop_kodon_1_ramec_P[ix1P+1] - stop_kodon_1_ramec_P[ix1P]) < limitORFP):
                            continue
                        elif ((stop_kodon_1_ramec_P[ix1P+1] - stop_kodon_1_ramec_P[ix1P]) >= limitORFP):
                            listORFsP_1 += [basesP[(stop_kodon_1_ramec_P[ix1P]+3):(stop_kodon_1_ramec_P[ix1P+1]+3)]]
                            listORFsP_1_pozice += [[stop_kodon_1_ramec_P[ix1P], stop_kodon_1_ramec_P[ix1P+1]]]
                    #print(listORFsP_1)
                    #print('Pozice ORFu v zadanem vlakne v 1. ramci jsou (v kodu: "listORFsP_1_pozice"): \n',listORFsP_1_pozice)

                if (len(stop_kodon_2_ramec_P) < 2):
                    print('Nedostecny pocet STOP kodonu v 2. ramci v zadanem vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
                else:
                    for ix2P in range(0,len(stop_kodon_2_ramec_P)-1,1):
                        if ((stop_kodon_2_ramec_P[ix2P+1] - stop_kodon_2_ramec_P[ix2P]) < limitORFP):
                            continue
                        elif ((stop_kodon_2_ramec_P[ix2P+1] - stop_kodon_2_ramec_P[ix2P]) >= limitORFP):
                            listORFsP_2 += [basesP[(stop_kodon_2_ramec_P[ix2P]+3):(stop_kodon_2_ramec_P[ix2P+1]+3)]]
                            listORFsP_2_pozice += [[stop_kodon_2_ramec_P[ix2P], stop_kodon_2_ramec_P[ix2P+1]]]
                    #print(listORFsP_2)
                    #print('Pozice ORFu v zadanem vlakne ve 2. ramci jsou (v kodu: "listORFsP_2_pozice"): \n',listORFsP_2_pozice)            

                if (len(stop_kodon_3_ramec_P) < 2):
                    print('Nedostecny pocet STOP kodonu v 3. ramci v zadanem vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
                else:
                    for ix3P in range(0,len(stop_kodon_3_ramec_P)-1,1):
                        if ((stop_kodon_3_ramec_P[ix3P+1] - stop_kodon_3_ramec_P[ix3P]) < limitORFP):
                            continue
                        elif ((stop_kodon_3_ramec_P[ix3P+1] - stop_kodon_3_ramec_P[ix3P]) >= limitORFP):
                            listORFsP_3 += [basesP[(stop_kodon_3_ramec_P[ix3P]+3):(stop_kodon_3_ramec_P[ix3P+1]+3)]]
                            listORFsP_3_pozice += [[stop_kodon_3_ramec_P[ix3P], stop_kodon_3_ramec_P[ix3P+1]]]
                    #print(listORFsP_3)
                    #print('Pozice ORFu v zadanem vlakne ve 3. ramci jsou (v kodu: "listORFsP_3_pozice"): \n',listORFsP_3_pozice)     

                listORFsP = listORFsP_1 + listORFsP_2 + listORFsP_3
                #print('Pocet ORFu pro min. delku ORFu ',limitORFP,'je: ',len(listORFsP))
                #print('Seznam ORFu (pro "P"), tj. zadane vlakno: ',listORFsP)

                listORFsP_pozice = listORFsP_1_pozice + listORFsP_2_pozice + listORFsP_3_pozice ## !!! ve skutecnosti se hodi <-- toto
                #sorted_listORFsP_pozice = sorted(listORFsP_pozice) # bude se hodit pro vyhledavani pozice pouziteho ORFu pro vysetrovany peptid/protein

                basesRN_2 = basesRN
                
                for iBN in basesRN:
                  if ((iBN !=str("A"))&(iBN !=str("G"))&(iBN !=str("T"))&(iBN !=str("C"))):
                    basesRN_2 = basesRN.replace(iBN,'')
                  else:
                    continue
                    #basesRN_2 = basesRN.replace(iBN,iBN)

                basesN = basesRN_2

                lenbN = len(basesN)
                ####print("Nukleotidu (delka NA) je po redukci o znaky pro novy radek: ",lenbN)

                #####
                #print('Vypis komplementarni sekvence k zadane sekvenci:')
                #print(basesN)
                #####

                # vyhledani ORFu (v komplementarnim vlakne)
                # def.2 podle:
                '''
                Def. 2: ORF je usek DNA mezi dvemi STOP kodony a pocet jeho nukleotidu je deliteny 3 beze zbytku
                '''
                stop_codon_findN = [mNstop.start() for mNstop in re.finditer('(TAA|TAG|TGA)', basesN)]
                #print('Stop kodon v komplementarnim vlakne je na techto pozicich: \n',stop_codon_findN)
                #print('Stop kodon je tolikrat:',len(stop_codon_findN))
                ## najde opravdu, kde zacina dana sekvence, tedy kde je "T"
                ## vyhodi cislo, kde je prave to pismeno, napr. je-li "T" 6. zleva pak: "5"

                limitORFN = limitORF # DULEZITY nastavitelny parametr

                stop_kodon_1_ramec_N = []
                stop_kodon_2_ramec_N = []
                stop_kodon_3_ramec_N = []

                listORFsN = []

                listORFsN_1 = []
                listORFsN_2 = []
                listORFsN_3 = []

                listORFsN_1_pozice = []
                listORFsN_2_pozice = []
                listORFsN_3_pozice = []

                for ixyzN in stop_codon_findN:
                    if ((ixyzN)%3 == 0):
                        stop_kodon_1_ramec_N += [ixyzN]
                    elif ((ixyzN+2)%3 == 0):
                        stop_kodon_2_ramec_N += [ixyzN]
                    elif ((ixyzN+1)%3 == 0):
                        stop_kodon_3_ramec_N += [ixyzN]

                #print('stop_kodon_1_ramec: \n',stop_kodon_1_ramec_N)
                #print('stop_kodon_2_ramec: \n',stop_kodon_2_ramec_N)
                #print('stop_kodon_3_ramec: \n',stop_kodon_3_ramec_N)


                if (len(stop_kodon_1_ramec_N) < 2):
                    print('Nedostecny pocet STOP kodonu v 1. ramci v komplemtarnim vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
                else:
                    for ix1N in range(0,len(stop_kodon_1_ramec_N)-1,1):
                        if ((stop_kodon_1_ramec_N[ix1N+1] - stop_kodon_1_ramec_N[ix1N]) < limitORFN):
                            continue
                        elif ((stop_kodon_1_ramec_N[ix1N+1] - stop_kodon_1_ramec_N[ix1N]) >= limitORFN):
                            listORFsN_1 += [basesN[(stop_kodon_1_ramec_N[ix1N]+3):(stop_kodon_1_ramec_N[ix1N+1]+3)]]
                            listORFsN_1_pozice += [[stop_kodon_1_ramec_N[ix1N], stop_kodon_1_ramec_N[ix1N+1]]]
                    #print(listORFsN_1)
                    #print('Pozice ORFu v komplementarnim vlakne v 1. ramci jsou (v kodu: "listORFsN_1_pozice"): \n',listORFsN_1_pozice)

                if (len(stop_kodon_2_ramec_N) < 2):
                    print('Nedostecny pocet STOP kodonu v 2. ramci v komplemtarnim vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
                else:
                    for ix2N in range(0,len(stop_kodon_2_ramec_N)-1,1):
                        if ((stop_kodon_2_ramec_N[ix2N+1] - stop_kodon_2_ramec_N[ix2N]) < limitORFN):
                            continue
                        elif ((stop_kodon_2_ramec_N[ix2N+1] - stop_kodon_2_ramec_N[ix2N]) >= limitORFN):
                            listORFsN_2 += [basesN[(stop_kodon_2_ramec_N[ix2N]+3):(stop_kodon_2_ramec_N[ix2N+1]+3)]]
                            listORFsN_2_pozice += [[stop_kodon_2_ramec_N[ix2N], stop_kodon_2_ramec_N[ix2N+1]]]
                    #print(listORFsN_2)
                    #print('Pozice ORFu v komplementarnim vlakne ve 2. ramci jsou (v kodu: "listORFsN_2_pozice"): \n',listORFsN_2_pozice)            

                if (len(stop_kodon_3_ramec_N) < 2):
                    print('Nedostecny pocet STOP kodonu v 3. ramci v komplemtarnim vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
                else:
                    for ix3N in range(0,len(stop_kodon_3_ramec_N)-1,1):
                        if ((stop_kodon_3_ramec_N[ix3N+1] - stop_kodon_3_ramec_N[ix3N]) < limitORFN):
                            continue
                        elif ((stop_kodon_3_ramec_N[ix3N+1] - stop_kodon_3_ramec_N[ix3N]) >= limitORFN):
                            listORFsN_3 += [basesN[(stop_kodon_3_ramec_N[ix3N]+3):(stop_kodon_3_ramec_N[ix3N+1]+3)]]
                            listORFsN_3_pozice += [[stop_kodon_3_ramec_N[ix3N], stop_kodon_3_ramec_N[ix3N+1]]]
                    #print(listORFsN_3)
                    #print('Pozice ORFu v komplementarnim vlakne ve 3. ramci jsou (v kodu: "listORFsN_3_pozice"): \n',listORFsN_3_pozice)    

                listORFsN = listORFsN_1 + listORFsN_2 + listORFsN_3
                #print('Pocet ORFu pro min. delku ORFu ',limitORFN,'je: ',len(listORFsN))
                #print('Seznam ORFu (pro "N"), tj. komplemtarni vlakno: ',listORFsN)

                listORFsN_pozice = listORFsN_1_pozice + listORFsN_2_pozice + listORFsN_3_pozice ## !!! ve skutecnosti se hodi <-- toto
                #sorted_listORFsP_pozice = sorted(listORFsP_pozice) # bude se hodit pro vyhledavani pozice pouziteho ORFu pro vysetrovany peptid/protein

                #print('\n')

        elif (proces123 == 3):

            stop_kodon_1_ramec_P = []
            stop_kodon_2_ramec_P = []
            stop_kodon_3_ramec_P = []

            listORFsP = []
            #listORFsN = [] ##### - je nize ... 

            listORFsP_1 = []
            listORFsP_2 = []
            listORFsP_3 = []

            listORFsP_1_pozice = []
            listORFsP_2_pozice = []
            listORFsP_3_pozice = []

            odpad_ORF_P_1 = []
            odpad_ORF_P_2 = []
            odpad_ORF_P_3 = []

            odpad_ORF_P_1_pozice = []
            odpad_ORF_P_2_pozice = []
            odpad_ORF_P_3_pozice = []


            for ixyzP in stop_codon_findP:
                if ((ixyzP)%3 == 0):
                    stop_kodon_1_ramec_P += [ixyzP]
                elif ((ixyzP+2)%3 == 0):
                    stop_kodon_2_ramec_P += [ixyzP]
                elif ((ixyzP+1)%3 == 0):
                    stop_kodon_3_ramec_P += [ixyzP]

            #print('Zadane vlakno - stop_kodon_1_ramec: \n',stop_kodon_1_ramec_P)
            #print('Zadane vlakno - stop_kodon_2_ramec: \n',stop_kodon_2_ramec_P)
            #print('Zadane vlakno - stop_kodon_3_ramec: \n',stop_kodon_3_ramec_P)

            basesPW = basesP

            for iReplaceP in basesP:
                if ((iReplaceP != str("A"))&(iReplaceP != str("G"))&(iReplaceP != str("T"))&(iReplaceP != str("C")))==True:
                    basesPW = basesP.replace(iReplaceP,str("W"))
                else:
                    continue
                    #basesPW = basesP.replace(iReplaceP,iReplaceP)
            ##print(basesPW)

            #patternP = 'W+'

            if (len(stop_kodon_1_ramec_P) < 2):
                print('Nedostecny pocet STOP kodonu v 1. ramci v zadanem vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix1P in range(0,len(stop_kodon_1_ramec_P)-1,1):
                    if ((stop_kodon_1_ramec_P[ix1P+1] - stop_kodon_1_ramec_P[ix1P]) < limitORFP):
                        continue
                    elif ((stop_kodon_1_ramec_P[ix1P+1] - stop_kodon_1_ramec_P[ix1P]) >= limitORFP):
                        for iP1_nns in [basesPW[(stop_kodon_1_ramec_P[ix1P]+3):(stop_kodon_1_ramec_P[ix1P+1]+3)]]:
                            if re.search('W+', iP1_nns) == None:
                                listORFsP_1 += [basesPW[(stop_kodon_1_ramec_P[ix1P]+3):(stop_kodon_1_ramec_P[ix1P+1]+3)]]
                                listORFsP_1_pozice += [[stop_kodon_1_ramec_P[ix1P], stop_kodon_1_ramec_P[ix1P+1]]]
                            else:
                                odpad_ORF_P_1 += [basesPW[(stop_kodon_1_ramec_P[ix1P]+3):(stop_kodon_1_ramec_P[ix1P+1]+3)]]
                                odpad_ORF_P_1_pozice += [[stop_kodon_1_ramec_P[ix1P], stop_kodon_1_ramec_P[ix1P+1]]]
                            
                    else:
                        break
                #print(listORFsP_1)
                #print('Pozice ORFu v zadanem vlakne v 1. ramci jsou (v kodu: "listORFsP_1_pozice"): \n',listORFsP_1_pozice)

            if (len(stop_kodon_2_ramec_P) < 2):
                print('Nedostecny pocet STOP kodonu v 2. ramci v zadanem vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix2P in range(0,len(stop_kodon_2_ramec_P)-1,1):
                    if ((stop_kodon_2_ramec_P[ix2P+1] - stop_kodon_2_ramec_P[ix2P]) < limitORFP):
                        continue
                    elif ((stop_kodon_2_ramec_P[ix2P+1] - stop_kodon_2_ramec_P[ix2P]) >= limitORFP):
                        for iP2_nns in [basesPW[(stop_kodon_2_ramec_P[ix2P]+3):(stop_kodon_2_ramec_P[ix2P+1]+3)]]:
                            if re.search('W+', iP2_nns) == None:
                                listORFsP_2 += [basesPW[(stop_kodon_2_ramec_P[ix2P]+3):(stop_kodon_2_ramec_P[ix2P+1]+3)]]
                                listORFsP_2_pozice += [[stop_kodon_2_ramec_P[ix2P], stop_kodon_2_ramec_P[ix2P+1]]]
                            else:
                                odpad_ORF_P_2 += [basesPW[(stop_kodon_2_ramec_P[ix2P]+3):(stop_kodon_2_ramec_P[ix2P+1]+3)]]
                                odpad_ORF_P_2_pozice += [[stop_kodon_2_ramec_P[ix2P], stop_kodon_2_ramec_P[ix2P+1]]]
                    else:
                        break
                #print(listORFsP_2)
                #print('Pozice ORFu v zadanem vlakne ve 2. ramci jsou (v kodu: "listORFsP_2_pozice"): \n',listORFsP_2_pozice)            


            if (len(stop_kodon_3_ramec_P) < 2):
                print('Nedostecny pocet STOP kodonu v 3. ramci v zadanem vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix3P in range(0,len(stop_kodon_3_ramec_P)-1,1):
                    if ((stop_kodon_3_ramec_P[ix3P+1] - stop_kodon_3_ramec_P[ix3P]) < limitORFP):
                        continue
                    elif ((stop_kodon_3_ramec_P[ix3P+1] - stop_kodon_3_ramec_P[ix3P]) >= limitORFP):
                        for iP3_nns in [basesPW[(stop_kodon_3_ramec_P[ix3P]+3):(stop_kodon_3_ramec_P[ix3P+1]+3)]]: # mozna by melo byt: ":(stop_kodon_3_ramec_P[ix3P+1]+2)", tj. "..+2..--> ale pak neni len(ORF)mod3 = 0
                            if re.search('W+', iP3_nns) == None:
                                listORFsP_3 += [basesPW[(stop_kodon_3_ramec_P[ix3P]+3):(stop_kodon_3_ramec_P[ix3P+1]+3)]]
                                listORFsP_3_pozice += [[stop_kodon_3_ramec_P[ix3P], stop_kodon_3_ramec_P[ix3P+1]]]
                            else:
                                odpad_ORF_P_3 += [basesPW[(stop_kodon_3_ramec_P[ix3P]+3):(stop_kodon_3_ramec_P[ix3P+1]+3)]]
                                odpad_ORF_P_3_pozice += [[stop_kodon_3_ramec_P[ix3P], stop_kodon_3_ramec_P[ix3P+1]]]
                    else:
                        break
                #print(listORFsP_3)
                #print('Pozice ORFu v zadanem vlakne ve 3. ramci jsou (v kodu: "listORFsP_3_pozice"): \n',listORFsP_3_pozice)     

            listORFsP = listORFsP_1 + listORFsP_2 + listORFsP_3
            #print('Pocet ORFu pro min. delku ORFu ',limitORFP,'je: ',len(listORFsP))
            #print('Seznam ORFu (pro "P"), tj. zadane vlakno: ',listORFsP)

            listORFsP_pozice = listORFsP_1_pozice + listORFsP_2_pozice + listORFsP_3_pozice ## !!! ve skutecnosti se hodi <-- toto
            #sorted_listORFsP_pozice = sorted(listORFsP_pozice) # bude se hodit pro vyhledavani pozice pouziteho ORFu pro vysetrovany peptid/protein
            
            ####print("Nukleotidu (delka NA) je po redukci o znaky pro novy radek: ",lenbN)

            #####
            #print('Vypis komplementarni sekvence k zadane sekvenci:')
            #print(basesN)
            #####

            # vyhledani ORFu (v komplementarnim vlakne)
            # def.2 podle:
            '''
            Def. 2: ORF je usek DNA mezi dvemi STOP kodony a pocet jeho nukleotidu je deliteny 3 beze zbytku
            '''
            stop_codon_findN = [mNstop.start() for mNstop in re.finditer('(TAA|TAG|TGA)', basesN)]
            #print('Stop kodon v komplementarnim vlakne je na techto pozicich: \n',stop_codon_findN)
            #print('Stop kodon je tolikrat:',len(stop_codon_findN))
            ## najde opravdu, kde zacina dana sekvence, tedy kde je "T"
            ## vyhodi cislo, kde je prave to pismeno, napr. je-li "T" 6. zleva pak: "5"

            #limitORF = 30 # DULEZITY nastavitelny parametr # - uz tam je ...
            limitORFN = limitORF # DULEZITY nastavitelny parametr

            stop_kodon_1_ramec_N = []
            stop_kodon_2_ramec_N = []
            stop_kodon_3_ramec_N = []

            listORFsN = []

            listORFsN_1 = []
            listORFsN_2 = []
            listORFsN_3 = []

            listORFsN_1_pozice = []
            listORFsN_2_pozice = []
            listORFsN_3_pozice = []

            odpad_ORF_N_1 = []
            odpad_ORF_N_2 = []
            odpad_ORF_N_3 = []

            odpad_ORF_N_1_pozice = []
            odpad_ORF_N_2_pozice = []
            odpad_ORF_N_3_pozice = []


            for ixyzN in stop_codon_findN:
                if ((ixyzN)%3 == 0):
                    stop_kodon_1_ramec_N += [ixyzN]
                elif ((ixyzN+2)%3 == 0):
                    stop_kodon_2_ramec_N += [ixyzN]
                elif ((ixyzN+1)%3 == 0):
                    stop_kodon_3_ramec_N += [ixyzN]

            #print('stop_kodon_1_ramec: \n',stop_kodon_1_ramec_N)
            #print('stop_kodon_2_ramec: \n',stop_kodon_2_ramec_N)
            #print('stop_kodon_3_ramec: \n',stop_kodon_3_ramec_N)

            basesNW = basesN

            for iReplaceN in basesN:
                if ((iReplaceN != str("A"))&(iReplaceN != str("G"))&(iReplaceN != str("T"))&(iReplaceN != str("C")))==True:
                    basesNW = basesN.replace(iReplaceN,str("W"))
                else:
                    continue
                    #basesNW = basesN.replace(iReplaceN,iReplaceN)
            ##print(basesNW)

            #patternN = 'W+'

            if (len(stop_kodon_1_ramec_N) < 2):
                print('Nedostecny pocet STOP kodonu v 1. ramci v komplemtarnim vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix1N in range(0,len(stop_kodon_1_ramec_N)-1,1):
                    if ((stop_kodon_1_ramec_N[ix1N+1] - stop_kodon_1_ramec_N[ix1N]) < limitORFN):
                        continue
                    elif ((stop_kodon_1_ramec_N[ix1N+1] - stop_kodon_1_ramec_N[ix1N]) >= limitORFN):
                        for iN1_nns in [basesNW[(stop_kodon_1_ramec_N[ix1N]+3):(stop_kodon_1_ramec_N[ix1N+1]+3)]]:
                             if re.search('W+', iN1_nns) == None:
                                listORFsN_1 += [basesNW[(stop_kodon_1_ramec_N[ix1N]+3):(stop_kodon_1_ramec_N[ix1N+1]+3)]]
                                listORFsN_1_pozice += [[stop_kodon_1_ramec_N[ix1N], stop_kodon_1_ramec_N[ix1N+1]]]
                             else:
                                odpad_ORF_N_1 += [basesNW[(stop_kodon_1_ramec_N[ix1N]+3):(stop_kodon_1_ramec_N[ix1N+1]+3)]]
                                odpad_ORF_N_1_pozice += [[stop_kodon_1_ramec_N[ix1N], stop_kodon_1_ramec_N[ix1N+1]]]
                    else:
                        break
                #print(listORFsN_1)
                #print('Pozice ORFu v komplementarnim vlakne v 1. ramci jsou (v kodu: "listORFsN_1_pozice"): \n',listORFsN_1_pozice)

            if (len(stop_kodon_2_ramec_N) < 2):
                print('Nedostecny pocet STOP kodonu v 2. ramci v komplemtarnim vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix2N in range(0,len(stop_kodon_2_ramec_N)-1,1):
                    if ((stop_kodon_2_ramec_N[ix2N+1] - stop_kodon_2_ramec_N[ix2N]) < limitORFN):
                        continue
                    elif ((stop_kodon_2_ramec_N[ix2N+1] - stop_kodon_2_ramec_N[ix2N]) >= limitORFN):
                        for iN2_nns in [basesNW[(stop_kodon_2_ramec_N[ix2N]+3):(stop_kodon_2_ramec_N[ix2N+1]+3)]]:
                            if re.search('W+', iN2_nns) == None:
                                listORFsN_2 += [basesNW[(stop_kodon_2_ramec_N[ix2N]+3):(stop_kodon_2_ramec_N[ix2N+1]+3)]]
                                listORFsN_2_pozice += [[stop_kodon_2_ramec_N[ix2N], stop_kodon_2_ramec_N[ix2N+1]]]
                            else:
                                odpad_ORF_N_2 += [basesNW[(stop_kodon_2_ramec_N[ix2N]+3):(stop_kodon_2_ramec_N[ix2N+1]+3)]]
                                odpad_ORF_N_2_pozice += [[stop_kodon_2_ramec_N[ix2N], stop_kodon_2_ramec_N[ix2N+1]]]
                    else:
                        break
                #print(listORFsN_2)
                #print('Pozice ORFu v komplementarnim vlakne ve 2. ramci jsou (v kodu: "listORFsN_2_pozice"): \n',listORFsN_2_pozice)            

            if (len(stop_kodon_3_ramec_N) < 2):
                print('Nedostecny pocet STOP kodonu v 3. ramci v komplemtarnim vlakne: je 1 nebo zadny! Jsou treba alespon 2.')
            else:
                for ix3N in range(0,len(stop_kodon_3_ramec_N)-1,1):
                    if ((stop_kodon_3_ramec_N[ix3N+1] - stop_kodon_3_ramec_N[ix3N]) < limitORFN):
                        continue
                    elif ((stop_kodon_3_ramec_N[ix3N+1] - stop_kodon_3_ramec_N[ix3N]) >= limitORFN):
                        for iN3_nns in [basesNW[(stop_kodon_3_ramec_N[ix3N]+3):(stop_kodon_3_ramec_N[ix3N+1]+3)]]:
                            if re.search('W+', iN3_nns) == None:
                                listORFsN_3 += [basesNW[(stop_kodon_3_ramec_N[ix3N]+3):(stop_kodon_3_ramec_N[ix3N+1]+3)]]
                                listORFsN_3_pozice += [[stop_kodon_3_ramec_N[ix3N], stop_kodon_3_ramec_N[ix3N+1]]]
                            else:
                                odpad_ORF_N_3 += [basesNW[(stop_kodon_3_ramec_N[ix3N]+3):(stop_kodon_3_ramec_N[ix3N+1]+3)]]
                                odpad_ORF_N_3_pozice += [[stop_kodon_3_ramec_N[ix3N], stop_kodon_3_ramec_N[ix3N+1]]]
                    else:
                        break
                #print(listORFsN_3)
                #print('Pozice ORFu v komplementarnim vlakne ve 3. ramci jsou (v kodu: "listORFsN_3_pozice"): \n',listORFsN_3_pozice)    

            listORFsN = listORFsN_1 + listORFsN_2 + listORFsN_3
            #print('Pocet ORFu pro min. delku ORFu ',limitORFN,'je: ',len(listORFsN))
            #print('Seznam ORFu (pro "N"), tj. komplemtarni vlakno: ',listORFsN)

            listORFsN_pozice = listORFsN_1_pozice + listORFsN_2_pozice + listORFsN_3_pozice ## !!! ve skutecnosti se hodi <-- toto
            #sorted_listORFsP_pozice = sorted(listORFsP_pozice) # bude se hodit pro vyhledavani pozice pouziteho ORFu pro vysetrovany peptid/protein

        else:
            print('Pro pokracovani je treba zadat: 1 nebo 2 nebo 3.')
            quit()


        #print('\n')

        #print('Minimalni delka ORFu je: ',limitORF,' nukleotidu.') # Vcetne start a stop kodonu.
        #if countP:
        #print('Pocet ORFu v zadanem vlaknu je: ', len(listORFsP))
        #if countN:
        #print('Pocet ORFu v komplementarnim vlaknu je: ', len(listORFsN))

        #print('Celkovy pocet ORFu v obou vlaknech je: ', len(listORFsP) + len(listORFsN))
        #print('\n')
        '''
        print('ORFy v zadanem vlaknu:')
        print('\n')
        print(listORFsP)
        print('\n')
        print('ORFy v komplementarnim vlaknu:')
        print('\n')
        print(listORFsN)
        print('\n')
        '''

        # Nyni preklad do peptidu (udane kody, pak zadane vlakno, zn. "P", pak komplementarni vlakno, zn. "N"):

        # Standardni geneticky kod:
        Base1 = ['TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG']
        Base2 = ['TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG']
        Base3 = ['TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG']
        Standard_Code = ['FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']

        # Alternativni geneticke kody:
        Vertebrate_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG']
        Yeast_Mitochondrial_Code = ['FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Mold_Protozoan_Coelenterate_Mitochondrial_and_Mycoplasma_Spiroplasma_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Invertebrate_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG']
        Ciliate_Dasycladacean_and_Hexamita_Nuclear_Code = ['FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Echinoderm_and_Flatworm_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG']
        Euplotid_Nuclear_Code = ['FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Bacterial_Archaeal_and_Plant_Plastid_Code = ['FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Alternative_Yeast_Nuclear_Code = ['FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Ascidian_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG']
        Alternative_Flatworm_Mitochondrial_Code = ['FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG']
        Blepharisma_Nuclear_Code = ['FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Chlorophycean_Mitochondrial_Code = ['FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Trematode_Mitochondrial_Code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG']
        Scenedesmus_obliquus_mitochondrial_Code = ['FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Thraustochytrium_Mitochondrial_Code = ['FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']
        Pterobranchia_mitochondrial_code = ['FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG']
        Candidate_Division_SR1_and_Gracilibacteria_Code = ['FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG']

        listBasesAll = Base1 + Base2 + Base3

        ##### Nejprve pro zadane vlakno ("P"):

        peptidP_c00 = []
        peptidP_c01 = []
        peptidP_c02 = []
        peptidP_c03 = []
        peptidP_c04 = []
        peptidP_c05 = []
        peptidP_c06 = []
        peptidP_c07 = []
        peptidP_c08 = []
        peptidP_c09 = []
        peptidP_c10 = []
        peptidP_c11 = []
        peptidP_c12 = []
        peptidP_c13 = []
        peptidP_c14 = []
        peptidP_c15 = []
        peptidP_c16 = []
        peptidP_c17 = []
        peptidP_c18 = []

        #peptidP = []
        #PeptidyP = []


        if (len(listORFsP) == 0):
            print('Zadny ORF a peptid v zadanem vlaknu.')
            print('\n')
        else:

            for iorfxP in range(len(listORFsP)):
             for kodonP in range(int(len(listORFsP[iorfxP][:])/3)):

              # print('kodon je: ',kodon+1)
              res1P = [i1P for i1P in range(len(listBasesAll[0])) if listBasesAll[0].startswith(listORFsP[iorfxP][kodonP+0+2*kodonP], i1P)]
              res2P = [i2P for i2P in range(len(listBasesAll[1])) if listBasesAll[1].startswith(listORFsP[iorfxP][kodonP+1+2*kodonP], i2P)]
              res3P = [i3P for i3P in range(len(listBasesAll[2])) if listBasesAll[2].startswith(listORFsP[iorfxP][kodonP+2+2*kodonP], i3P)]

              set_res1P = set(res1P)
              set_res2P = set(res2P)
              set_res3P = set(res3P)

              #print('set_res1P je: ',set_res1P)
              #print('set_res2P je: ',set_res2P)
              #print('set_res3P je: ',set_res3P)

              prunik12P = set_res1P.intersection(set_res2P)
              prunik123P = prunik12P.intersection(set_res3P)

              # print('prunik je: ',prunik123P)

              indexAAsP = int(list(prunik123P)[0])

              AMKP_c00 = Standard_Code[0][indexAAsP]

              #print('AMK_c00 je: ',AMKP_c00)

              peptidP_c00 += [AMKP_c00]

              sPepP_c00 = ''
              seqPepP_c00 = peptidP_c00
              JseqPepP_c00 = sPepP_c00.join(seqPepP_c00)
              strPepP_c00 = JseqPepP_c00

              countP_STOP_c00 = strPepP_c00.count('*')
              PepsP_c00 = strPepP_c00.split('*',countP_STOP_c00)

              #PepsP = strPepP.split('*',len(listORFsP))
              #PepsP = PepsP[0:-1]
              
              #if AMKP:
              #print('peptidy jsou: ', PepsP)
              #PeptidyP += [PepsP]

              #preramec = []

              ###

              OpeptidesP_c00 = open("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa","w+")
              for iOpP_c00 in range(len(PepsP_c00)-1):
                  preramecP_c00 = listORFsP_pozice[iorfxP][0]
                  if ((preramecP_c00)%3 == 0):
                      frameP_c00 = 1
                  elif ((preramecP_c00 + 1)%3 == 0):
                      frameP_c00 = 3
                  elif ((preramecP_c00 + 2)%3 == 0):
                      frameP_c00 = 2
                  else:
                      print('Nenasel jsem cteci ramec...')
                  #strORF = str(ORF)
                  #OpeptidesP_c00.write(">" + "frame: " + str(frameP_c00) + ". ," + "ORF-position: " + str(listORFsP_pozice[iOpP_c00]) + "\n" + PepsP_c00[iOpP_c00] + "\n")
                  OpeptidesP_c00.write(">"+str(listORFsP_pozice[iOpP_c00][0])+" "+str(listORFsP_pozice[iOpP_c00][1]) + "\n" + PepsP_c00[iOpP_c00] + "\n")
              OpeptidesP_c00.close()

              ### XXX      

              '''
              with open('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.fasta','w+') as fPepPfasta_c00:
                for Pf_c00 in range(len(PepsP_c00)-1):
                    fPepPfasta_c00.write("> ORF: " + str(listORFsP_pozice[Pf_c00]) + "\n" + PepsP_c00[Pf_c00] + "\n")
              fPepPfasta_c00.close()
              '''
            
              ### znacka: "alternativa P (zacatek)"

              #
              AMKP_c08 = Bacterial_Archaeal_and_Plant_Plastid_Code[0][indexAAsP]

              peptidP_c08 += [AMKP_c08]
              sPepP_c08 = ''
              seqPepP_c08 = peptidP_c08
              JseqPepP_c08 = sPepP_c08.join(seqPepP_c08)
              strPepP_c08 = JseqPepP_c08

              countP_STOP_c08 = strPepP_c08.count('*')
              PepsP_c08 = strPepP_c08.split('*',countP_STOP_c08)
              
              OpeptidesP_c08 = open("FASTA_output/peptidesP/peptidesP_c08_Bacterial_Archaeal_and_Plant_Plastid_Code.faa","w+")
              for iOpP_c08 in range(len(PepsP_c08)-1):
                    preramecP_c08 = listORFsP_pozice[iorfxP][0]
                    if ((preramecP_c08)%3 == 0):
                        frameP_c08 = 1
                    elif ((preramecP_c08 + 1)%3 == 0):
                        frameP_c08 = 3
                    elif ((preramecP_c08 + 2)%3 == 0):
                        frameP_c08 = 2
                    else:
                        print('Nenasel jsem cteci ramec...')
                        
                    #OpeptidesP_c08.write(">" + "frame: " + str(frameP_c08) + ". ," + "ORF-position: " + str(listORFsP_pozice[iOpP_c08]) + "\n" + PepsP_c08[iOpP_c08] + "\n")
                    OpeptidesP_c08.write(">"+str(listORFsP_pozice[iOpP_c08][0])+" "+str(listORFsP_pozice[iOpP_c08][1]) + "\n" + PepsP_c08[iOpP_c08] + "\n")
              OpeptidesP_c08.close()

              #
              AMKP_c09 = Alternative_Yeast_Nuclear_Code[0][indexAAsP]

              peptidP_c09 += [AMKP_c09]
              sPepP_c09 = ''
              seqPepP_c09 = peptidP_c09
              JseqPepP_c09 = sPepP_c09.join(seqPepP_c09)
              strPepP_c09 = JseqPepP_c09

              countP_STOP_c09 = strPepP_c09.count('*')
              PepsP_c09 = strPepP_c09.split('*',countP_STOP_c09)
              
              OpeptidesP_c09 = open("FASTA_output/peptidesP/peptidesP_c09_Alternative_Yeast_Nuclear_Code.faa","w+")
              for iOpP_c09 in range(len(PepsP_c09)-1):
                    preramecP_c09 = listORFsP_pozice[iorfxP][0]
                    if ((preramecP_c09)%3 == 0):
                        frameP_c09 = 1
                    elif ((preramecP_c09 + 1)%3 == 0):
                        frameP_c09 = 3
                    elif ((preramecP_c09 + 2)%3 == 0):
                        frameP_c09 = 2
                    else:
                        print('Nenasel jsem cteci ramec...')
                        
                    #OpeptidesP_c09.write(">" + "frame: " + str(frameP_c09) + ". ," + "ORF-position: " + str(listORFsP_pozice[iOpP_c09]) + "\n" + PepsP_c09[iOpP_c09] + "\n")
                    OpeptidesP_c09.write(">"+str(listORFsP_pozice[iOpP_c09][0])+" "+str(listORFsP_pozice[iOpP_c09][1]) + "\n" + PepsP_c09[iOpP_c09] + "\n")
              OpeptidesP_c09.close()


              # ukladam sekvence peptiduu pro zadane vlakno ("P"):
              # pri kazdem behu programu se soubor cely prepise...
              # ukladam jako txt a "fasta"
              
              with open('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.txt','w+') as fPepP_c00:
                  for P_c00 in PepsP_c00:
                      fPepP_c00.write(str(P_c00) + '\n')
                  fPepP_c00.close()
              '''
              #zaloha prechoziho postupu pro "fasta" format:

              with open('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.fasta','w+') as fPepPfasta:
                  for Pf in PepsP:
                      fPepPfasta.write(str(Pf) + '\n')
                  fPepPfasta.close()
              '''
              OpeptidesP_c00 = open("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.fasta", "w")
              for iOpP_c00 in range(len(PepsP_c00)):
                  OpeptidesP_c00.write(">" + "\n" +PepsP_c00[iOpP_c00] + "\n")
              #napoveda pro hlavicku:
              #mame 2 seznamy:
              # list_seq = [sequence1, sequence2, sequence3, sequence4]
              # list_name = [name1, name2, name3, name4]
              #udelame dle:
              # OpeptidesP.write(">" + list_name[i] + "\n" +list_seq[i] + "\n")
              OpeptidesP_c00.close()

              #print('Seznam peptidu pro zadane vlakno je v TXT souboru /home/oem/Documents/PacesDr/peptidesP.txt')
              #print('Seznam peptidu pro zadane vlakno je ve FASTA souboru /home/oem/Documents/PacesDr/peptidesP.fasta')

        ### Nyni preklad do peptidu pro komplementarni vlakno ("N"):

        print('\n')

        peptidN_c00 = []
        peptidN_c01 = []
        peptidN_c02 = []
        peptidN_c03 = []
        peptidN_c04 = []
        peptidN_c05 = []
        peptidN_c06 = []
        peptidN_c07 = []
        peptidN_c08 = []
        peptidN_c09 = []
        peptidN_c10 = []
        peptidN_c11 = []
        peptidN_c12 = []
        peptidN_c13 = []
        peptidN_c14 = []
        peptidN_c15 = []
        peptidN_c16 = []
        peptidN_c17 = []
        peptidN_c18 = []

        #peptidN = []
        #PeptidyN = []

        if (len(listORFsN) == 0):
            print('Zadny ORF a peptid v komplementarnim vlaknu.')
            print('\n')
        else:

            for iorfxN in range(len(listORFsN)):
             for kodonN in range(int(len(listORFsN[iorfxN][:])/3)):
                 
              res1N = [i1N for i1N in range(len(listBasesAll[0])) if listBasesAll[0].startswith(listORFsN[iorfxN][kodonN+0+2*kodonN], i1N)]
              res2N = [i2N for i2N in range(len(listBasesAll[1])) if listBasesAll[1].startswith(listORFsN[iorfxN][kodonN+1+2*kodonN], i2N)]
              res3N = [i3N for i3N in range(len(listBasesAll[2])) if listBasesAll[2].startswith(listORFsN[iorfxN][kodonN+2+2*kodonN], i3N)]

              set_res1N = set(res1N)
              set_res2N = set(res2N)
              set_res3N = set(res3N)

              #print('set_res1N je: ',set_res1N)
              #print('set_res2N je: ',set_res2N)
              #print('set_res3N je: ',set_res3N)

              prunik12N = set_res1N.intersection(set_res2N)
              prunik123N = prunik12N.intersection(set_res3N)

              # print('prunik je: ',prunik123N)

              indexAAsN = int(list(prunik123N)[0])

              # pro standardni kod (komplementarni, "N" vlakno):

              #
              AMKN_c00 = Standard_Code[0][indexAAsN]
              # print('AMK je: ',AMKN_c00)

              peptidN_c00 += [AMKN_c00]

              sPepN_c00 = ''
              seqPepN_c00 = peptidN_c00
              JseqPepN_c00 = sPepN_c00.join(seqPepN_c00)
              strPepN_c00 = JseqPepN_c00

              countN_STOP_c00 = strPepN_c00.count('*')
              PepsN_c00 = strPepN_c00.split('*',countN_STOP_c00)
              
              OpeptidesN_c00 = open("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa","w+")
              for iOpN_c00 in range(len(PepsN_c00)-1):
                    preramecN_c00 = listORFsN_pozice[iorfxN][0]
                    if ((preramecN_c00)%3 == 0):
                        frameN_c00 = 1
                    elif ((preramecN_c00 + 1)%3 == 0):
                        frameN_c00 = 3
                    elif ((preramecN_c00 + 2)%3 == 0):
                        frameN_c00 = 2
                    else:
                        print('Nenasel jsem cteci ramec...')

                    #OpeptidesN_c00.write(">" + "frame: " + str(frameN_c00) + ". ," + "ORF-position: " + str(listORFsN_pozice[iOpN_c00]) + "\n" + PepsN_c00[iOpN_c00] + "\n")
                    OpeptidesN_c00.write(">"+str(listORFsN_pozice[iOpN_c00][0])+" "+str(listORFsN_pozice[iOpN_c00][1]) + "\n" + PepsN_c00[iOpN_c00] + "\n")
              OpeptidesN_c00.close()

              #
              AMKN_c08 = Bacterial_Archaeal_and_Plant_Plastid_Code[0][indexAAsN]

              peptidN_c08 += [AMKN_c08]
              sPepN_c08 = ''
              seqPepN_c08 = peptidN_c08
              JseqPepN_c08 = sPepN_c08.join(seqPepN_c08)
              strPepN_c08 = JseqPepN_c08

              countN_STOP_c08 = strPepN_c08.count('*')
              PepsN_c08 = strPepN_c08.split('*',countN_STOP_c08)
              
              OpeptidesN_c08 = open("FASTA_output/peptidesN/peptidesN_c08_Bacterial_Archaeal_and_Plant_Plastid_Code.faa","w+")
              for iOpN_c08 in range(len(PepsN_c08)-1):
                    preramecN_c08 = listORFsN_pozice[iorfxN][0]
                    if ((preramecN_c08)%3 == 0):
                        frameN_c08 = 1
                    elif ((preramecN_c08 + 1)%3 == 0):
                        frameN_c08 = 3
                    elif ((preramecN_c08 + 2)%3 == 0):
                        frameN_c08 = 2
                    else:
                        print('Nenasel jsem cteci ramec...')
                        
                    #OpeptidesN_c08.write(">" + "frame: " + str(frameN_c08) + ". ," + "ORF-position: " + str(listORFsN_pozice[iOpN_c08]) + "\n" + PepsN_c08[iOpN_c08] + "\n")
                    OpeptidesN_c08.write(">"+str(listORFsN_pozice[iOpN_c08][0])+" "+str(listORFsN_pozice[iOpN_c08][1]) + "\n" + PepsN_c08[iOpN_c08] + "\n")
              OpeptidesN_c08.close()

              #
              AMKN_c09 = Alternative_Yeast_Nuclear_Code[0][indexAAsN]

              peptidN_c09 += [AMKN_c09]
              sPepN_c09 = ''
              seqPepN_c09 = peptidN_c09
              JseqPepN_c09 = sPepN_c09.join(seqPepN_c09)
              strPepN_c09 = JseqPepN_c09

              countN_STOP_c09 = strPepN_c09.count('*')
              PepsN_c09 = strPepN_c09.split('*',countN_STOP_c09)
              
              OpeptidesN_c09 = open("FASTA_output/peptidesN/peptidesN_c09_Alternative_Yeast_Nuclear_Code.faa","w+")
              for iOpN_c09 in range(len(PepsN_c09)-1):
                    preramecN_c09 = listORFsN_pozice[iorfxN][0]
                    if ((preramecN_c09)%3 == 0):
                        frameN_c09 = 1
                    elif ((preramecN_c09 + 1)%3 == 0):
                        frameN_c09 = 3
                    elif ((preramecN_c09 + 2)%3 == 0):
                        frameN_c09 = 2
                    else:
                        print('Nenasel jsem cteci ramec...')
                        
                    #OpeptidesN_c09.write(">" + "frame: " + str(frameN_c09) + ". ," + "ORF-position: " + str(listORFsN_pozice[iOpN_c09]) + "\n" + PepsN_c09[iOpN_c09] + "\n")
                    OpeptidesN_c09.write(">"+str(listORFsN_pozice[iOpN_c09][0])+" "+str(listORFsN_pozice[iOpN_c09][1]) + "\n" + PepsN_c09[iOpN_c09] + "\n")
              OpeptidesN_c09.close()

              # ukladam sekvence peptiduu pro komplementarni vlakno ("N"):
              # pri kazdem behu programu se soubor cely prepise...
              #  ukladam jako txt a "fasta"

              with open('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.txt','w+') as fPepN_c00:
                  for N_c00 in PepsN_c00:
                      fPepN_c00.write(str(N_c00) + '\n')
                  fPepN_c00.close()
              '''
              #zaloha prechoziho postupu pro "fasta" format:

              with open('/home/oem/Documents/PacesDr/peptidesN.fasta','w+') as fPepNfasta:
                  for Nf in PepsN:
                      fPepNfasta.write(str(Nf) + '\n')
                  fPepNfasta.close()
              '''
              OpeptidesN_c00 = open("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.fasta", "w")
              for iOpN_c00 in range(len(PepsN_c00)):
                  OpeptidesN_c00.write(">" + "\n" +PepsN_c00[iOpN_c00] + "\n")
              #napoveda pro hlavicku:
              #mame 2 seznamy:
              # list_seq = [sequence1, sequence2, sequence3, sequence4]
              # list_name = [name1, name2, name3, name4]
              #udelame dle:# OpeptidesP.write(">" + list_name[i] + "\n" +list_seq[i] + "\n")
              OpeptidesN_c00.close()

              #print('\n')

              #print('Seznam peptidu pro komplementarni vlakno je v TXT souboru /home/oem/Documents/PacesDr/peptidesN.txt')
              #print('Seznam peptidu pro komplementarni vlakno je ve FASTA souboru /home/oem/Documents/PacesDr/peptidesN.fasta')

        count_TM_dvoj_obrazku_P = 0
        count_TM_dvoj_obrazku_N = 0

        PepsP_c00_4finalA = []
        
        if (len(listORFsP) == 0):
            print('Zadny ORF a peptid v zadanem vlaknu.')
            print('\n')
        else:
            for irpP_c00 in range(len(PepsP_c00)):
              if (len(PepsP_c00[irpP_c00]) <= 0):
                  #print(irpP_c00,'je pod limitem')
                  continue
              elif (len(PepsP_c00[irpP_c00]) > 0):
                  PepsP_c00_4finalA += [PepsP_c00[irpP_c00]]
                  #print(irpP_c00,'je nad limitem')
              else:
                  break
        
        PepsN_c00_4finalA = []

        
        if (len(listORFsN) == 0):
            print('Zadny ORF a peptid v komplementarnim vlaknu.')
            print('\n')
        else:
            for irpN_c00 in range(len(PepsN_c00)):
              if (len(PepsN_c00[irpN_c00]) <= 0):
                  #print(irpN_c00,'je pod limitem')
                  continue
              elif (len(PepsN_c00[irpN_c00]) > 0):
                  PepsN_c00_4finalA += [PepsN_c00[irpN_c00]]
                  #print(irpN_c00,'je nad limitem')
              else:
                  break
        
        import os.path

        import tmhmm

        import matplotlib.pyplot as plt
        import matplotlib.text
        import matplotlib.cm as cm
        import numpy as np
        import matplotlib.mlab as mlab
        from matplotlib.artist import Artist

        from matplotlib.lines import Line2D
        from matplotlib.patches import Rectangle
        from matplotlib.text import Text
        from matplotlib.image import AxesImage

        import tkinter

        if (len(listORFsP) == 0):
            print('Zadny ORF a peptid v zadanem vlaknu.')
            listTMsP_pozice = []
            print('\n')
        else:

            listTMsP_pozice = []

            count_M_P = []

            apM_P_coords = []

            for irpP_c00_4fA in range(len(PepsP_c00_4finalA)):

                annotationP, posteriorP = tmhmm.predict(PepsP_c00_4finalA[irpP_c00_4fA],'3','TM_model/TMHMM20.model')

                print(annotationP, posteriorP)

                #print(listORFsP_pozice[irpP_c00_4fA])

                apP = annotationP, posteriorP

                #import matplotlib.pyplot as plt
                #import numpy as np

                data1P = []
                data2P = []
                data3P = []

                for iapP in range(len((annotationP, posteriorP)[0])):
                    data1P += [apP[1][iapP][0]]
                    data2P += [apP[1][iapP][1]]
                    data3P += [apP[1][iapP][2]]

                    if (((len(data1P) & len(data2P) & len(data3P)) == len((annotationP, posteriorP)[0])) & (str('M') in ((annotationP, posteriorP)[0]))):

                      apP0 = (annotationP, posteriorP)[0]
                      apM_P = [apMP for apMP in range(len(apP0)) if apP0.startswith('M', apMP)]
                      min_apM_P = min(apM_P)
                      max_apM_P = max(apM_P)
                      apM_P_coords += [[min_apM_P,max_apM_P]]

                      count_TM_dvoj_obrazku_P += 1 # pocet ORFu s TM sek. strukturou, pocita oba obr. (tento a "linkovy") pro jeden OFR s TM jako jeden
                      #print('pocet ORFu s TM sekundarni strukturou (zadane vlakno): ', count_TM_dvoj_obrazku_P)

                      count_M_P_i = ((annotationP, posteriorP)[0]).count('M')

                      XaxP = np.arange(0,len(PepsP_c00_4finalA[irpP_c00_4fA]),1.0)

                      #plt.plot(XaxP,data1P,'g--',label='intracellular',data2P,'b--',label='in membrane',data3P,'y--',label='extracellular')
                      plt.plot(XaxP,data1P,'g--', label='intracellular')
                      plt.plot(XaxP,data2P,'b--', label='in membrane')
                      plt.plot(XaxP,data3P,'y--', label='extracellular')
                      plt.xlabel('poradi aminokyseliny (N-C konec)')
                      plt.ylabel('pravdepodobnost umisteni')
                      plt.title('TM predikce' + ' - pozice ORFu: ' + str(listORFsP_pozice[irpP_c00_4fA]) + ' bp.')
                      plt.legend()
                      #plt.show()
                      plt.savefig('figures/given/Fig No. ' + str(listORFsP_pozice[irpP_c00_4fA]) + ' provided_string.png')
                      #plt.savefig('/home/oem/Documents/_html/figures/Fig No. ' + str([irpP_c00_4fA]) + ' provided_string.png')
                      plt.close()#

                      # "linkovy obr.":
                      xsP = np.linspace(1,1,lenbP)
                      plt.figure(figsize=(lenbP/120,lenbP/4800))
                      plt.hlines(y=1,xmin=1,xmax=lenbP,color='blue',linestyle='-',lw=20)
                      plt.hlines(y=1,xmin=(listORFsP_pozice[irpP_c00_4fA][0]),xmax=(listORFsP_pozice[irpP_c00_4fA][1]),color='r',linestyle='-',lw=30)
                      plt.yticks([])
                      plt.xticks([])
                      '''
                      if ((listORFsP_pozice[irpP_c00_4fA][0])%3 == 0):
                          frameP_ORFall_lin = 1
                      elif ((listORFsP_pozice[irpP_c00_4fA][0] + 1)%3 == 0):
                          frameP_ORFall_lin = 3
                      elif ((listORFsP_pozice[irpP_c00_4fA][0] + 2)%3 == 0):
                          frameP_ORFall_lin = 2
                      else:
                          continue
                      '''
                      #plt.annotate(str(listORFsP_pozice[irpP_c00_4fA][0]),(listORFsP_pozice[irpP_c00_4fA][0],1),fontsize=16)
                      #plt.annotate(str(listORFsP_pozice[irpP_c00_4fA][1]),(listORFsP_pozice[irpP_c00_4fA][1],1),fontsize=16)
                      #plt.annotate(str(lenbP),(lenbP,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
                      #plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
                      '''
                      plt.annotate('Fr.:'+str(frameP_ORFall_lin),((listORFsP_pozice[irpP_c00_4fA][0]+listORFsP_pozice[irpP_c00_4fA][1])/2,1),fontsize=16,color='green',verticalalignment='top')
                      '''
                      #plt.savefig('/home/oem/Documents/_html/figures/NA string - fig No. ' + str([irpP_c00_4fA]) + ' provided_string.png')
                      plt.savefig('figures/given/NA string - fig No. ' + str(listORFsP_pozice[irpP_c00_4fA]) + ' provided_string.png')
                      #plt.show()
                      plt.close()#

                      listTMsP_pozice += [listORFsP_pozice[irpP_c00_4fA]]

                      count_M_P += [[listORFsP_pozice[irpP_c00_4fA][0], count_M_P_i]]
                      
                      ## listTMsP_pozice += listORFsP_pozice[irpP_c00_4fA]
                      ## pro souhrn typu [2310, 2460, 3150, 3306, 8037, 8202, 8304, 8460, 2368, 2581, 5236, 5443, 6997, 7456, 10210, 10441, 1847, 2000, 3035, 3260, 4175, 4394, 7958, 8330]
                      ## tj. bez vzeti do dvojic
            #print('ORFy s TM strukturou, zadane vlakno:\n',listTMsP_pozice)

        if (len(listORFsN) == 0):
            print('Zadny ORF a peptid v komplementarnim vlaknu.')
            listTMsN_pozice = []
            print('\n')
        else:

            listTMsN_pozice = []

            count_M_N = []

            apM_N_coords = []

            for irpN_c00_4fA in range(len(PepsN_c00_4finalA)):

                annotationN, posteriorN = tmhmm.predict(PepsN_c00_4finalA[irpN_c00_4fA],'4','TM_model/TMHMM20.model')

                print(annotationN, posteriorN)

                apN = annotationN, posteriorN

                #import matplotlib.pyplot as plt
                #import numpy as np

                data1N = []
                data2N = []
                data3N = []

                for iapN in range(len((annotationN, posteriorN)[0])):
                    data1N += [apN[1][iapN][0]]
                    data2N += [apN[1][iapN][1]]
                    data3N += [apN[1][iapN][2]]

                    if (((len(data1N) & len(data2N) & len(data3N)) == len((annotationN, posteriorN)[0])) & (str('M') in ((annotationN, posteriorN)[0]))):

                      apN0 = (annotationN, posteriorN)[0]
                      apM_N = [apMN for apMN in range(len(apN0)) if apN0.startswith('M', apMN)]
                      min_apM_N = min(apM_N)
                      max_apM_N = max(apM_N)
                      apM_N_coords += [[min_apM_N,max_apM_N]]

                      count_TM_dvoj_obrazku_N += 1 # pocet ORFu s TM sek. strukturou, pocita oba obr. pro jeden OFR s TM jako jeden
                      #print('pocet ORFu s TM sekundarni strukturou (komplementarni vlakno): ', count_TM_dvoj_obrazku_N)

                      count_M_N_i = ((annotationN, posteriorN)[0]).count('M')

                      XaxN = np.arange(0,len(PepsN_c00_4finalA[irpN_c00_4fA]),1.0)

                      #plt.plot(XaxN,data1N,'g--',label='intracellular',data2N,'b--',label='in membrane',data3N,'y--',label='extracellular')
                      plt.plot(XaxN,data1N,'g--', label='intracellular')
                      plt.plot(XaxN,data2N,'b--', label='in membrane')
                      plt.plot(XaxN,data3N,'y--', label='extracellular')
                      plt.xlabel('poradi aminokyseliny (N-C konec)')
                      plt.ylabel('pravdepodobnost umisteni')
                      plt.title('TM predikce' + ' - pozice ORFu: ' + str(listORFsN_pozice[irpN_c00_4fA]) + ' bp.')
                      plt.legend()
                      #plt.show()
                      #plt.savefig('/home/oem/Documents/_html/figures/Fig No. ' + str([irpN_c00_4fA]) + ' complementary_string.png')
                      plt.savefig('figures/complementary/Fig No. ' + str(listORFsN_pozice[irpN_c00_4fA]) + ' complementary_string.png')
                      plt.close()#

                      # "linkovy obr.":
                      xsN = np.linspace(1,1,lenbN)
                      plt.figure(figsize=(lenbN/120,lenbN/4800))
                      plt.hlines(y=1,xmin=1,xmax=lenbN,color='blue',linestyle='-',lw=20)
                      plt.hlines(y=1,xmin=(listORFsN_pozice[irpN_c00_4fA][0]),xmax=(listORFsN_pozice[irpN_c00_4fA][1]),color='r',linestyle='-',lw=30)
                      plt.yticks([])
                      plt.xticks([])
                      '''
                      if ((listORFsN_pozice[irpN_c00_4fA][0])%3 == 0):
                          frameN_ORFall_lin = 1
                      elif ((listORFsN_pozice[irpN_c00_4fA][0] + 1)%3 == 0):
                          frameN_ORFall_lin = 3
                      elif ((listORFsN_pozice[irpN_c00_4fA][0] + 2)%3 == 0):
                          frameN_ORFall_lin = 2
                      else:
                          continue
                      '''
                      #plt.annotate(str(listORFsN_pozice[irpN_c00_4fA][0]),(listORFsN_pozice[irpN_c00_4fA][0],1),fontsize=16)
                      #plt.annotate(str(listORFsN_pozice[irpN_c00_4fA][1]),(listORFsN_pozice[irpN_c00_4fA][1],1),fontsize=16)
                      #plt.annotate(str(lenbN),(lenbN,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
                      #plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
                      '''
                      plt.annotate('Fr.:'+str(frameN_ORFall_lin),((listORFsN_pozice[irpN_c00_4fA][0]+listORFsN_pozice[irpN_c00_4fA][1])/2,1),fontsize=16,color='green',verticalalignment='top')
                      '''
                      plt.savefig('figures/complementary/NA string - fig No. ' + str(listORFsN_pozice[irpN_c00_4fA]) + ' complementary_string.png')
                      #plt.savefig('/home/oem/Documents/_html/figures/NA string - fig No. ' + str([irpN_c00_4fA]) + ' complementary_string.png')
                      #plt.show()
                      plt.close()#

                      listTMsN_pozice += [listORFsN_pozice[irpN_c00_4fA]]

                      count_M_N += [[listORFsN_pozice[irpN_c00_4fA][0], count_M_N_i]]
                      
                      ## listTMsN_pozice += listORFsN_pozice[irpN_c00_4fA]
                      ## pro souhrn typu [2310, 2460, 3150, 3306, 8037, 8202, 8304, 8460, 2368, 2581, 5236, 5443, 6997, 7456, 10210, 10441, 1847, 2000, 3035, 3260, 4175, 4394, 7958, 8330]
                      ## tj. bez vzeti do dvojic
            #print('ORFy s TM strukturou, komplementarni vlakno:\n',listTMsN_pozice)

        # urcovani CC struktury (zadane i komlementarni vlakno) - zacatek:

        #CC_threshold_P = 0.10
        #CC_threshold_N = 0.10

        CC_count_P = 0
        CC_count_N = 0

        CC_ORFstart_P = []
        CC_ORFstart_N = []

        max_results_CC_P = []
        max_results_CC_N = []

        ind_MrkcP_coords = []
        ind_MrkcN_coords = []
        
        if (len(listORFsP) == 0):
            print('Zadny ORF a peptid v zadanem vlaknu. Nelze hledat CC strukturu.')
            print('\n')
        else:
            from deepcoil import DeepCoil
            from deepcoil.utils import plot_preds
            from Bio import SeqIO

            dc = DeepCoil(use_gpu=False)
            inp = {str(entry.id): str(entry.seq) for entry in SeqIO.parse('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa', 'fasta')}
            results = dc.predict(inp)

            for key in results.keys():
                plot_preds(results[key], out_file='CC_plots/plotsP/{}.png'.format(key))

            for key in results.keys():
                #print('Max z ['+key+'] CC (zadane vlakno) je: ',max(results[key]['cc']))
                if (max(results[key]['cc']) >= CC_threshold_P):
                    plot_preds(results[key], out_file='CC_plots/plotsP-realCC/{}.png'.format(key))
                    CC_count_P += 1
                    CC_ORFstart_P += [int(key)]
                    max_results_CC_P += [[int(key), max(results[key]['cc'])]]

                    MrkcP = max(results[key]['cc'])
                    list_rkcP = (results[key]['cc']).tolist()
                    ind_MrkcP = list_rkcP.index(MrkcP) # udava jen 1 cislo - prvni vyskyt z leva, cislovano od 0 (nuly)
                    ###min_ind_MrkcP = min(ind_MrkcP) # nelze - min jde ze seznamu, ne z jednoho cisla 'int' - neni treba - viz koment. k ind_MrkcP o 1 radek vyse
                    #print(ind_MrkcP)
                    ind_MrkcP_coords += [ind_MrkcP]

                else:
                    None
        
        if (len(listORFsN) == 0):
            print('Zadny ORF a peptid v komplementarnim vlaknu. Nelze hledat CC strukturu.')
            print('\n')
        else:
            from deepcoil import DeepCoil
            from deepcoil.utils import plot_preds
            from Bio import SeqIO

            dc = DeepCoil(use_gpu=False)
            inp = {str(entry.id): str(entry.seq) for entry in SeqIO.parse('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa', 'fasta')}
            results = dc.predict(inp)

            for key in results.keys():
                plot_preds(results[key], out_file='CC_plots/plotsN/{}.png'.format(key))

            for key in results.keys():
                #print('Max z ['+key+'] CC (komplementarni vlakno) je: ',max(results[key]['cc']))
                if (max(results[key]['cc']) >= CC_threshold_N):
                    plot_preds(results[key], out_file='CC_plots/plotsN-realCC/{}.png'.format(key))
                    CC_count_N += 1
                    CC_ORFstart_N += [int(key)]
                    max_results_CC_N += [[int(key), max(results[key]['cc'])]]

                    MrkcN = max(results[key]['cc'])
                    list_rkcN = (results[key]['cc']).tolist()
                    ind_MrkcN = list_rkcN.index(MrkcN) # udava jen 1 cislo - prvni vyskyt z leva, cislovano od 0 (nuly)
                    ###min_ind_MrkcN = min(ind_MrkcN) # nelze - min jde ze seznamu, ne z jednoho cisla 'int' - neni treba - viz koment. k ind_MrkcN o 1 radek vyse
                    #print(ind_MrkcN)
                    ind_MrkcN_coords += [ind_MrkcN]
                    
                else:
                    None
        
        # urcovani CC struktury (zadane i komlementarni vlakno) - konec.

        # urcovani obrazku/ORFu s jen CC strukturou (ale muze tam byt i GPI struktura), jen zadane, "P" vlakno - zacatek:

        listTMsP_pozice_0 = []

        if len(listTMsP_pozice) == 0 :
            print('Zadny ORF s TM strukturou v zadanem vlaknu.')

        for itmPp in range(len(listTMsP_pozice)):
            listTMsP_pozice_0 += [listTMsP_pozice[itmPp][0]]

        set_listTMsP_pozice_0 = set(listTMsP_pozice_0)
        set_CC_ORFstart_P = set(CC_ORFstart_P)

        CC_notTM_P_set = set_CC_ORFstart_P.difference(set_listTMsP_pozice_0)
        CC_notTM_P_list = list(CC_notTM_P_set)

        CC_notTM_P_coords = []

        for iOPp in listORFsP_pozice:
            for iCP in CC_notTM_P_list:
                if (iCP == iOPp[0]):
                    CC_notTM_P_coords += [iOPp]
                else:
                    None

        # cilem zde je hlavne: CC_notTM_P_coords
        
        TM_notCC_P_set = set_listTMsP_pozice_0.difference(set_CC_ORFstart_P)    

        TM_a_CC_P_set = set_listTMsP_pozice_0.intersection(set_CC_ORFstart_P)    

        # urcovani obrazku/ORFu s jen CC strukturou (ale muze tam byt i GPI struktura), jen zadane, "P" vlakno - konec.

        # urcovani ORFu s CC strukturou (neresi se zda je pritomna dalsi sek. struktura), zadane vlakno - zacatek:

        CC_all_P_coords = []

        for iOPcc in listORFsP_pozice:
            for iPcc in CC_ORFstart_P:
                if (iPcc == iOPcc[0]):
                    CC_all_P_coords += [iOPcc]
                else:
                    None

        # urcovani ORFu s CC strukturou (neresi se zda je pritomna dalsi sek. struktura), zadane vlakno - konec.

        # urcovani obrazku/ORFu s jen CC strukturou (ale muze tam byt i GPI struktura), jen komplementarni, "N" vlakno - zacatek:

        listTMsN_pozice_0 = []

        for itmNp in range(len(listTMsN_pozice)):
            listTMsN_pozice_0 += [listTMsN_pozice[itmNp][0]]

        set_listTMsN_pozice_0 = set(listTMsN_pozice_0)
        set_CC_ORFstart_N = set(CC_ORFstart_N)

        CC_notTM_N_set = set_CC_ORFstart_N.difference(set_listTMsN_pozice_0)
        CC_notTM_N_list = list(CC_notTM_N_set)

        CC_notTM_N_coords = []

        for iONp in listORFsN_pozice:
            for iCN in CC_notTM_N_list:
                if (iCN == iONp[0]):
                    CC_notTM_N_coords += [iONp]
                else:
                    None

        # cilem zde je hlavne: CC_notTM_N_coords
        
        TM_notCC_N_set = set_listTMsN_pozice_0.difference(set_CC_ORFstart_N)    

        TM_a_CC_N_set = set_listTMsN_pozice_0.intersection(set_CC_ORFstart_N)    

        # urcovani obrazku/ORFu s jen CC strukturou (ale muze tam byt i GPI struktura), jen komplementarni, "N" vlakno - konec.

        # urcovani ORFu s CC strukturou (neresi se zda je pritomna dalsi sek. struktura), komplementarni vlakno - zacatek:

        CC_all_N_coords = []

        for iONcc in listORFsN_pozice:
            for iNcc in CC_ORFstart_N:
                if (iNcc == iONcc[0]):
                    CC_all_N_coords += [iONcc]
                else:
                    None

        # urcovani ORFu s CC strukturou (neresi se zda je pritomna dalsi sek. struktura), komplementarni vlakno - konec.    

        # tvorba "linkoveho" obr. pro ORF (mozna) jen s CC (je mozna GPI) strukturou, zadane, "P" vlakno - zacatek:

        for iLccP in range(len(CC_notTM_P_coords)):

            xsP = np.linspace(1,1,lenbP)
            plt.figure(figsize=(lenbP/120,lenbP/4800))
            plt.hlines(y=1,xmin=1,xmax=lenbP,color='blue',linestyle='-',lw=20)
            plt.hlines(y=1,xmin=(CC_notTM_P_coords[iLccP][0]),xmax=(CC_notTM_P_coords[iLccP][1]),color='r',linestyle='-',lw=30)
            plt.yticks([])
            plt.xticks([])
            if ((CC_notTM_P_coords[iLccP][0])%3 == 0):
              frameP_iLccP = 1
            elif ((CC_notTM_P_coords[iLccP][0] + 1)%3 == 0):
              frameP_iLccP = 3
            elif ((CC_notTM_P_coords[iLccP][0] + 2)%3 == 0):
              frameP_iLccP = 2
            else:
              continue
            plt.annotate(str(CC_notTM_P_coords[iLccP][0]),(CC_notTM_P_coords[iLccP][0],1),fontsize=16)
            plt.annotate(str(CC_notTM_P_coords[iLccP][1]),(CC_notTM_P_coords[iLccP][1],1),fontsize=16)
            plt.annotate(str(lenbP),(lenbP,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
            plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
            #
            plt.annotate('Fr.:'+str(frameP_iLccP),((CC_notTM_P_coords[iLccP][0]+CC_notTM_P_coords[iLccP][1])/2,1),fontsize=16,color='green',verticalalignment='top')
            plt.savefig('CC_plots/plotsP-realCC/line-picture-CC-P/NA string - fig No. ' + str([iLccP]) + ' provided_string.png')
            #plt.show()
            plt.close()#

        # tvorba "linkoveho" obr. pro ORF (mozna) jen s CC (je mozna GPI) strukturou, zadane, "P" vlakno - konec.

        # tvorba "linkoveho" obr. pro ORF (mozna) jen s CC (je mozna GPI) strukturou, komplementarni, "N" vlakno - zacatek:

        for iLccN in range(len(CC_notTM_N_coords)):

            xsN = np.linspace(1,1,lenbN)
            plt.figure(figsize=(lenbN/120,lenbN/4800))
            plt.hlines(y=1,xmin=1,xmax=lenbN,color='blue',linestyle='-',lw=20)
            plt.hlines(y=1,xmin=(CC_notTM_N_coords[iLccN][0]),xmax=(CC_notTM_N_coords[iLccN][1]),color='r',linestyle='-',lw=30)
            plt.yticks([])
            plt.xticks([])
            if ((CC_notTM_N_coords[iLccN][0])%3 == 0):
              frameN_iLccN = 1
            elif ((CC_notTM_N_coords[iLccN][0] + 1)%3 == 0):
              frameN_iLccN = 3
            elif ((CC_notTM_N_coords[iLccN][0] + 2)%3 == 0):
              frameN_iLccN = 2
            else:
              continue
            plt.annotate(str(CC_notTM_N_coords[iLccN][0]),(CC_notTM_N_coords[iLccN][0],1),fontsize=16)
            plt.annotate(str(CC_notTM_N_coords[iLccN][1]),(CC_notTM_N_coords[iLccN][1],1),fontsize=16)
            plt.annotate(str(lenbN),(lenbN,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
            plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
            #
            plt.annotate('Fr.:'+str(frameN_iLccN),((CC_notTM_N_coords[iLccN][0]+CC_notTM_N_coords[iLccN][1])/2,1),fontsize=16,color='green',verticalalignment='top')
            plt.savefig('CC_plots/plotsN-realCC/line-picture-CC-N/NA string - fig No. ' + str([iLccN]) + ' complementary_string.png')
            #plt.show()
            plt.close()#

        # tvorba "linkoveho" obr. pro ORF (mozna) jen s CC (je mozna GPI) strukturou, komplementarni, "N" vlakno - konec.

        # 2D obrazky pro souhrn ORFu, pro zadane vlakno, po jednotlivych ctecich ramcich - zacatek :

        from urllib import *
        import json
        from matplotlib import *
        from tkinter import *
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import Rectangle
        from matplotlib.text import Text
        import matplotlib.image
        import numpy as np 

        import webbrowser

        ###

        if (len(listORFsP_1) == 0):
            print('Zadny ORF a peptid v zadanem vlaknu, v 1. ramci.')
            print('\n')
        else:

            xP1 = []
            yP1 = []

            for iAOP1 in range(len(listORFsP_1_pozice)):

                xP1 += [listORFsP_1_pozice[iAOP1][1]]
                yP1 += [listORFsP_1_pozice[iAOP1][0]]

            liste_xP1 = xP1
            liste_yP1 = yP1

            liste_stringP1 = liste_yP1

            fig = plt.figure()

            fig, ax = plt.subplots() #

            ax.set_xlim([0,lenbP])
            ax.set_ylim([0,lenbP])


            plt.figure(figsize=(lenbP/120,lenbP/120))
            plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
            plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
            plt.title('Souhrn ORFu - zadane vlakno \n (1. cteci ramec)', fontsize=round(lenbP/65))

            #ax = plt.gca()
            #fig = plt.gcf()
                
            for sxP1, dyP1 in zip(liste_xP1, liste_yP1):
                plt.annotate(' '+str(dyP1), xy = (sxP1, dyP1), fontsize=round(lenbP/120))
                plt.annotate('            -  '+str(sxP1), xy = (sxP1, dyP1), fontsize=round(lenbP/120))
                plt.scatter(liste_xP1, liste_stringP1, s = round(lenbP/5), c ='r') # picker=(True & 10000)

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/NA-ORFs-Pframe1-given_string-2D.svg')
            plt.savefig('figures/NA-ORFs-Pframe1-given_string-2D.png', dpi=32) # , dpi=16
            plt.close()

        ###

        if (len(listORFsP_2) == 0):
            print('Zadny ORF a peptid v zadanem vlaknu, v 2. ramci.')
            print('\n')
        else:

            xP2 = []
            yP2 = []

            for iAOP2 in range(len(listORFsP_2_pozice)):

                xP2 += [listORFsP_2_pozice[iAOP2][1]]
                yP2 += [listORFsP_2_pozice[iAOP2][0]]

            liste_xP2 = xP2
            liste_yP2 = yP2

            liste_stringP2 = liste_yP2

            fig = plt.figure()

            fig, ax = plt.subplots() ###

            plt.figure(figsize=(lenbP/120,lenbP/120))
            plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
            plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
            plt.title('Souhrn ORFu - zadane vlakno \n (2. cteci ramec)', fontsize=round(lenbP/65))

            ax = plt.gca()
            fig = plt.gcf()

            for sxP2, dyP2 in zip(liste_xP2, liste_yP2):
                plt.annotate(' '+str(dyP2), xy = (sxP2, dyP2), fontsize=round(lenbP/120))
                plt.annotate('            -  '+str(sxP2), xy = (sxP2, dyP2), fontsize=round(lenbP/120))
                plt.scatter(liste_xP2, liste_stringP2, s = round(lenbP/5), c ='b') # picker=(True & 10000)

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/NA-ORFs-Pframe2-given_string-2D.svg')
            plt.savefig('figures/NA-ORFs-Pframe2-given_string-2D.png', dpi=32) # , dpi=16
            plt.close()

        ###

        if (len(listORFsP_3) == 0):
            print('Zadny ORF a peptid v zadanem vlaknu, v 3. ramci.')
            print('\n')
        else:

            xP3 = []
            yP3 = []

            for iAOP3 in range(len(listORFsP_3_pozice)):

                xP3 += [listORFsP_3_pozice[iAOP3][1]]
                yP3 += [listORFsP_3_pozice[iAOP3][0]]

            liste_xP3 = xP3
            liste_yP3 = yP3

            liste_stringP3 = liste_yP3

            fig = plt.figure()

            fig, ax = plt.subplots() ###

            plt.figure(figsize=(lenbP/120,lenbP/120))
            plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
            plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
            plt.title('Souhrn ORFu - zadane vlakno \n (3. cteci ramec)', fontsize=round(lenbP/65))

            ax = plt.gca()
            fig = plt.gcf()

            for sxP3, dyP3 in zip(liste_xP3, liste_yP3):
                plt.annotate(' '+str(dyP3), xy = (sxP3, dyP3), fontsize=round(lenbP/120))
                plt.annotate('            -  '+str(sxP3), xy = (sxP3, dyP3), fontsize=round(lenbP/120))
                plt.scatter(liste_xP3, liste_stringP3, s = round(lenbP/5), c ='g') # picker=(True & 10000)


            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            ##plt.scatter(liste_x, liste_string, s = 10000, c ='r', alpha = 0.8).set_urls(['/home/oem/Documents/_html/index_123.html/#:~:text=['+str(s)+', '+str(d)+']'])
            #plt.plot(liste_x, liste_string,'x', markersize=128, color='black')

            #plt.show()
            #plt.savefig('figures/NA-ORFs-Pframe3-given_string-2D.svg')
            plt.savefig('figures/NA-ORFs-Pframe3-given_string-2D.png', dpi=32) # , dpi=16
            plt.close()

        # 2D obrazky pro souhrn ORFu, pro zadane vlakno, po jednotlivych ctecich ramcich - konec.

        # 2D obrazky pro souhrn ORFu, pro komplementarni vlakno, po jednotlivych ctecich ramcich - zacatek:
        # Tj. nyni obrazky s ORFy pro komplementarni vlakno - zacatek:

        if (len(listORFsN_1) == 0):
            print('Zadny ORF a peptid v komplementarnim vlaknu, v 1. ramci.')
            print('\n')
        else:

            xN1 = []
            yN1 = []

            for iAON1 in range(len(listORFsN_1_pozice)):

                xN1 += [listORFsN_1_pozice[iAON1][1]]
                yN1 += [listORFsN_1_pozice[iAON1][0]]

            liste_xN1 = xN1
            liste_yN1 = yN1

            liste_stringN1 = liste_yN1

            fig = plt.figure()

            fig, ax = plt.subplots() ###

            plt.figure(figsize=(lenbN/120,lenbN/120))
            plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
            plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
            plt.title('Souhrn ORFu - komplementarni vlakno \n (1. cteci ramec)', fontsize=round(lenbP/65))

            ax = plt.gca()
            fig = plt.gcf()
                
            for sxN1, dyN1 in zip(liste_xN1, liste_yN1):
                plt.annotate(' '+str(dyN1), xy = (sxN1, dyN1), fontsize=round(lenbP/120))
                plt.annotate('            -  '+str(sxN1), xy = (sxN1, dyN1), fontsize=round(lenbP/120))
                plt.scatter(liste_xN1, liste_stringN1, s = round(lenbP/5), c ='r') # picker=(True & 10000)

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/NA-ORFs-Nframe1-complementary_string-2D.svg')
            plt.savefig('figures/NA-ORFs-Nframe1-complementary_string-2D.png', dpi=32) # , dpi=16
            plt.close()

        ###

        if (len(listORFsN_2) == 0):
            print('Zadny ORF a peptid v komplementarnim vlaknu, v 2. ramci.')
            print('\n')
        else:

            xN2 = []
            yN2 = []

            for iAON2 in range(len(listORFsN_2_pozice)):

                xN2 += [listORFsN_2_pozice[iAON2][1]]
                yN2 += [listORFsN_2_pozice[iAON2][0]]

            liste_xN2 = xN2
            liste_yN2 = yN2

            liste_stringN2 = liste_yN2

            fig = plt.figure()

            fig, ax = plt.subplots() ###

            plt.figure(figsize=(lenbN/120,lenbN/120))
            plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
            plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
            plt.title('Souhrn ORFu - komplementarni vlakno \n (2. cteci ramec)', fontsize=round(lenbP/65))

            ax = plt.gca()
            fig = plt.gcf()

            for sxN2, dyN2 in zip(liste_xN2, liste_yN2):
                plt.annotate(' '+str(dyN2), xy = (sxN2, dyN2), fontsize=round(lenbP/120))
                plt.annotate('            -  '+str(sxN2), xy = (sxN2, dyN2), fontsize=round(lenbP/120))
                plt.scatter(liste_xN2, liste_stringN2, s = round(lenbP/5), c ='b') # picker=(True & 10000)

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/NA-ORFs-Nframe2-complementary_string-2D.svg')
            plt.savefig('figures/NA-ORFs-Nframe2-complementary_string-2D.png', dpi=32) # , dpi=16
            plt.close()

        ###

        if (len(listORFsN_3) == 0):
            print('Zadny ORF a peptid v komplementarnim vlaknu, v 3. ramci.')
            print('\n')
        else:

            xN3 = []
            yN3 = []

            for iAON3 in range(len(listORFsN_3_pozice)):

                xN3 += [listORFsN_3_pozice[iAON3][1]]
                yN3 += [listORFsN_3_pozice[iAON3][0]]

            liste_xN3 = xN3
            liste_yN3 = yN3

            liste_stringN3 = liste_yN3

            fig = plt.figure()

            fig, ax = plt.subplots() ###

            plt.figure(figsize=(lenbN/120,lenbN/120))
            plt.xlabel('Konec ORFu (v bp)', fontsize=round(lenbP/65))
            plt.ylabel('Pocatek ORFu (v bp)', fontsize=round(lenbP/65))
            plt.title('Souhrn ORFu - komplementarni vlakno \n (3. cteci ramec)', fontsize=round(lenbP/65))

            ax = plt.gca()
            fig = plt.gcf()

            for sxN3, dyN3 in zip(liste_xN3, liste_yN3):
                plt.annotate(' '+str(dyN3), xy = (sxN3, dyN3), fontsize=round(lenbP/120))
                plt.annotate('            -  '+str(sxN3), xy = (sxN3, dyN3), fontsize=round(lenbP/120))
                plt.scatter(liste_xN3, liste_stringN3, s = round(lenbP/5), c ='g') # picker=(True & 10000)


            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            ##plt.scatter(liste_x, liste_string, s = 10000, c ='r', alpha = 0.8).set_urls(['/home/oem/Documents/_html/index_123.html/#:~:text=['+str(s)+', '+str(d)+']'])
            #plt.plot(liste_x, liste_string,'x', markersize=128, color='black')

            #plt.show()
            #plt.savefig('figures/NA-ORFs-Nframe3-complementary_string-2D.svg')
            plt.savefig('figures/NA-ORFs-Nframe3-complementary_string-2D.png', dpi=32) # , dpi=16
            plt.close()

        # 2D obrazky pro souhrn ORFu, pro komplementarni vlakno, po jednotlivych ctecich ramcich - konec.
        # Tj. nyni pro obrazky s OFRy pro komplementarni vlako - konec.
        
        # puvodni 2D obr. pro souhrn vsech ORFu, nebo jednoho/1. ramce, zadane vlakno - konec.

        # urocvani  ORFu (pozice) s jen TM, s TM a CC, a jen s CC, sek. strukturou (mnozinove operace), zadane vlakno - zacatek:
        # dulezite vstupni promenne:
        # listTMsP_pozice
        # CC_all_P_coords
        # listORFsP_pozice

        listTMsP_pozice2_0 = []
        CC_all_P_coords_0 = []

        for itmPp2 in range(len(listTMsP_pozice)):
            listTMsP_pozice2_0 += [listTMsP_pozice[itmPp2][0]]

        set_listTMsP_pozice_0 = set(listTMsP_pozice2_0)

        for iccPp2 in range(len(CC_all_P_coords)):
            CC_all_P_coords_0 += [CC_all_P_coords[iccPp2][0]]    
        
        set_CC_all_P_coords_0 = set(CC_all_P_coords_0)
            
        set_only_TMsP_pozice_0 = set_listTMsP_pozice_0.difference(set_CC_all_P_coords_0)

        set_CC_only_P_coords_0 = set_CC_all_P_coords_0.difference(set_listTMsP_pozice_0)

        set_intersection_TM_CC_P_0 = set_listTMsP_pozice_0.intersection(set_CC_all_P_coords_0)

        intersection_TM_CC_P_0 = list(set_intersection_TM_CC_P_0)

        only_TMsP_pozice_0 = list(set_only_TMsP_pozice_0)

        CC_only_P_coords_0 = list(set_CC_only_P_coords_0)

        intersection_TM_CC_P = []

        for iOPx1 in listORFsP_pozice:
            for itmccP0 in intersection_TM_CC_P_0:
                if (itmccP0 == iOPx1[0]):
                    intersection_TM_CC_P += [iOPx1]
                else:
                    None

        only_TMsP_pozice = []

        for iOPx2 in listORFsP_pozice:
            for itmP0 in only_TMsP_pozice_0:
                if (itmP0 == iOPx2[0]):
                    only_TMsP_pozice += [iOPx2]
                else:
                    None

        CC_only_P_coords = []

        for iOPx3 in listORFsP_pozice:
            for iccP0 in CC_only_P_coords_0:
                if (iccP0 == iOPx3[0]):
                    CC_only_P_coords += [iOPx3]
                else:
                    None

        # ziskane promenne , hlavne:
        # intersection_TM_CC_P
        # only_TMsP_pozice
        # CC_only_P_coords

        # urocvani ORFu (pozice) s jen TM, s TM a CC, a jen s CC, sek. strukturou (mnozinove operace), zadane vlakno - konec.

        # urocvani ORFu (pozice) s jen TM, s TM a CC, a jen s CC, sek. strukturou (mnozinove operace), komplementarni vlakno - zacatek:
        # dulezite vstupni promenne:
        # listTMsN_pozice
        # CC_all_N_coords
        # listORFsN_pozice

        listTMsN_pozice2_0 = []
        CC_all_N_coords_0 = []
        
        for itmNp2 in range(len(listTMsN_pozice)):
            listTMsN_pozice2_0 += [listTMsN_pozice[itmNp2][0]]

        set_listTMsN_pozice_0 = set(listTMsN_pozice2_0)

        for iccNp2 in range(len(CC_all_N_coords)):
            CC_all_N_coords_0 += [CC_all_N_coords[iccNp2][0]]    
        
        set_CC_all_N_coords_0 = set(CC_all_N_coords_0)
            
        set_only_TMsN_pozice_0 = set_listTMsN_pozice_0.difference(set_CC_all_N_coords_0)

        set_CC_only_N_coords_0 = set_CC_all_N_coords_0.difference(set_listTMsN_pozice_0)

        set_intersection_TM_CC_N_0 = set_listTMsN_pozice_0.intersection(set_CC_all_N_coords_0)

        intersection_TM_CC_N_0 = list(set_intersection_TM_CC_N_0)    

        only_TMsN_pozice_0 = list(set_only_TMsN_pozice_0)

        CC_only_N_coords_0 = list(set_CC_only_N_coords_0)

        intersection_TM_CC_N = []

        for iONx1 in listORFsN_pozice:
            for itmccN0 in intersection_TM_CC_N_0:
                if (itmccN0 == iONx1[0]):
                    intersection_TM_CC_N += [iONx1]
                else:
                    None

        only_TMsN_pozice = []

        for iONx2 in listORFsN_pozice:
            for itmN0 in only_TMsN_pozice_0:
                if (itmN0 == iONx2[0]):
                    only_TMsN_pozice += [iONx2]
                else:
                    None

        CC_only_N_coords = []

        for iONx3 in listORFsN_pozice:
            for iccN0 in CC_only_N_coords_0:
                if (iccN0 == iONx3[0]):
                    CC_only_N_coords += [iONx3]
                else:
                    None

        # ziskane promenne , hlavne:
        # intersection_TM_CC_N
        # only_TMsN_pozice
        # CC_only_N_coords
        
        # urcovani ORFu (pozice) s jen TM, s TM a CC, a jen s CC, sek. strukturou (mnozinove operace), komplementarni vlakno - konec.

        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - zacatek:
        # dava souradnice ORFu a CC struktury
        
        if (len(CC_all_P_coords) == 0):
            print('Zadny ORF s CC sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xPccp = []
            yPccp = []

            for iPccp in range(len(CC_all_P_coords)):

                xPccp += [max_results_CC_P[iPccp][0] + 3*ind_MrkcP_coords[iPccp] + 3]
                yPccp += [round(max_results_CC_P[iPccp][1]*1000)/10]

            liste_xPccp = xPccp
            liste_yPccp = yPccp

            liste_xPccp = [0] + liste_xPccp
            liste_xPccp = liste_xPccp + [lenbP]

            liste_yPccp = [0] + liste_yPccp
            liste_yPccp = liste_yPccp + [0]

            liste_stringPccp = liste_yPccp

            fig = plt.figure()
            
            fig, ax = plt.subplots() ###

            plt.figure(figsize=(round(lenbP/120),round(lenbP/75)))
            plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, zadane vlakno', fontsize=round(lenbP/65))
            plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbP/65))
            #plt.title('Souhrn ORFu s CC strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxPccp, dyPccp in zip(liste_xPccp, liste_yPccp):
                if (sxPccp == 0):
                    plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                    plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                elif (sxPccp > 0) and (sxPccp < lenbP):
                    for sxPccpX, dyPccpX in zip(liste_xPccp[1:-1], liste_yPccp[1:-1]):
                        plt.annotate(''+str(dyPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifPccX = liste_xPccp[1:-1].index(sxPccpX)
                        plt.annotate('     - ORF:'+str(sxPccpX - 3*ind_MrkcP_coords[ifPccX] - 3)+' - CC:'+str(sxPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxPccp == lenbP):
                    plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                    plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                else:
                    None
                #plt.scatter(liste_xPccp, liste_stringPccp, s = lenbP/1.5, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xPccp, liste_stringPccp, width = 10, align='center', color='black') # width = lenbP/120
                    
            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno.svg')
            plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno.png', dpi=32) # , dpi=16
            plt.close()

        # dava souradnice ORFu a CC struktury
        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - zacatek:
        # dava souradnice ORFu a CC struktury
        # pro BIG verzi - roztazeni obr. pres celou str.
        
        if (len(CC_all_P_coords) == 0):
            print('Zadny ORF s CC sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xPccp = []
            yPccp = []

            for iPccp in range(len(CC_all_P_coords)):

                xPccp += [max_results_CC_P[iPccp][0] + 3*ind_MrkcP_coords[iPccp] + 3]
                yPccp += [round(max_results_CC_P[iPccp][1]*1000)/10]

            liste_xPccp = xPccp
            liste_yPccp = yPccp

            liste_xPccp = [0] + liste_xPccp
            liste_xPccp = liste_xPccp + [lenbP]

            liste_yPccp = [0] + liste_yPccp
            liste_yPccp = liste_yPccp + [0]

            liste_stringPccp = liste_yPccp

            fig = plt.figure()
            
            fig, ax = plt.subplots() ###

            plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))
            plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, zadane vlakno', fontsize=round(lenbP/65))
            plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbP/65))
            #plt.title('Souhrn ORFu s CC strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxPccp, dyPccp in zip(liste_xPccp, liste_yPccp):
                if (sxPccp == 0):
                    plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                    plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                elif (sxPccp > 0) and (sxPccp < lenbP):
                    for sxPccpX, dyPccpX in zip(liste_xPccp[1:-1], liste_yPccp[1:-1]):
                        plt.annotate(''+str(dyPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifPccX = liste_xPccp[1:-1].index(sxPccpX)
                        plt.annotate('     - ORF:'+str(sxPccpX - 3*ind_MrkcP_coords[ifPccX] - 3)+' - CC:'+str(sxPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxPccp == lenbP):
                    plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                    plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                else:
                    None
                #plt.scatter(liste_xPccp, liste_stringPccp, s = lenbP/1.5, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xPccp, liste_stringPccp, width = 10, align='center', color='black') # width = lenbP/120
                    
            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG.svg')
            plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG.png', dpi=32) # , dpi=16
            plt.close()

        # pro BIG verzi - roztazeni obr. pres celou str.
        # dava souradnice ORFu a CC struktury
        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - zacatek:
        # dava souradnice ORFu a CC struktury
        # pro BIG verzi - roztazeni obr. pres celou str.
        # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro CC
        
        if (len(CC_all_P_coords) == 0):
            print('Zadny ORF s CC sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:

            ## lidsky lokus
            ## popis: >hg38_dna range=chr19:17389648-17406217 5'pad=0 3'pad=0 strand=+
            #CCs = [[816,924],[1780,1847],[2028,2089],[2393,2445],[12236,12383],[14371,14468],[14874,15007],[15558,15574]]

            ## kureci lokus
            ## popis: >galGal6_dna range=chr28:3531163-3538110 5'pad=0 3'pad=0 strand=-
            #CCs = [[2122,2314],[2387,2484],[2601,2671],[2741,2802],[3384,3597],[4366,4463],[4545,4678],[4760,4806]]

            ## mysi lokus
            ## popis: >mm10_dna range=chr8:71519771-71538962 5'pad=0 3'pad=0 strand=-
            #CCs = [[1712,1835],[2528,2604],[3137,3177],[4170,4204],[13428,13578],[16406,16503],[16991,17124],[17580,17626]]

            ## luskoun ostrovni
            ## popis:
            ## >ref|NW_023435982.1|:679097-704212 Manis javanica isolate MJ74 unplaced genomic scaffold, YNU_ManJav_2.0 scaffold_88, whole genome shotgun sequence
            #CCs = [[5987,6080],[6935,7011],[7258,7319],[7581,7603],[15298,15439],[16489,16586],[16958,17091],[17551,17591]]

            ## kalon vabivy
            ## popis:
            ## >ref|NW_006429864.1|:923553-939337 Pteropus alecto unplaced genomic scaffold, ASM32557v1 scaffold160, whole genome shotgun sequence
            CCs = [[4140,4263],[5090,5166],[5375,5436],[5686,5717],[10693,10843],[11931,12021]]

            xPccp = []
            yPccp = []

            for iPccp in range(len(CC_all_P_coords)):

                xPccp += [max_results_CC_P[iPccp][0] + 3*ind_MrkcP_coords[iPccp] + 3]
                yPccp += [round(max_results_CC_P[iPccp][1]*1000)/10]

            liste_xPccp = xPccp
            liste_yPccp = yPccp

            liste_xPccp = [0] + liste_xPccp
            liste_xPccp = liste_xPccp + [lenbP]

            liste_yPccp = [0] + liste_yPccp
            liste_yPccp = liste_yPccp + [0]

            liste_stringPccp = liste_yPccp

            fig = plt.figure()
            
            fig, ax = plt.subplots() ###

            plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))
            plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, zadane vlakno', fontsize=round(lenbP/65))
            plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbP/65))
            #plt.title('Souhrn ORFu s CC strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxPccp, dyPccp in zip(liste_xPccp, liste_yPccp):
                if (sxPccp == 0):
                    plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                    plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                elif (sxPccp > 0) and (sxPccp < lenbP):
                    for sxPccpX, dyPccpX in zip(liste_xPccp[1:-1], liste_yPccp[1:-1]):
                        plt.annotate(''+str(dyPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifPccX = liste_xPccp[1:-1].index(sxPccpX)
                        plt.annotate('     - ORF:'+str(sxPccpX - 3*ind_MrkcP_coords[ifPccX] - 3)+' - CC:'+str(sxPccpX), xy = (sxPccpX, dyPccpX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxPccp == lenbP):
                    plt.annotate(' '+str(dyPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95))
                    plt.annotate('     -  '+str(sxPccp), xy = (sxPccp, dyPccp), fontsize=round(lenbP/95), rotation=90)
                else:
                    None
                #plt.scatter(liste_xPccp, liste_stringPccp, s = lenbP/1.5, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xPccp, liste_stringPccp, width = 10, align='center', color='black') # width = lenbP/120

            for iCCs in range(len(CCs)):
                plt.hlines(y=0,xmin=(CCs[iCCs][0]),xmax=(CCs[iCCs][1]),color='blue',linestyle='-',lw=1200)
                    
            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG-L.svg')
            plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG-L.png', dpi=32) # , dpi=16
            plt.close()

        # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro CC
        # pro BIG verzi - roztazeni obr. pres celou str.
        # dava souradnice ORFu a CC struktury
        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,zadane vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,komplementarni vlakno, vsechny ramce - zacatek:
        # dava souradnice ORFu a CC struktury
        
        if (len(CC_all_N_coords) == 0):
            print('Zadny ORF s CC sekundarni strukturou v komplementarnim vlaknu.')
            print('\n')
        else:
            xNccp = []
            yNccp = []

            for iNccp in range(len(CC_all_N_coords)):

                xNccp += [max_results_CC_N[iNccp][0] + 3*ind_MrkcN_coords[iNccp] + 3]
                yNccp += [round(max_results_CC_N[iNccp][1]*1000)/10]

            liste_xNccp = xNccp
            liste_yNccp = yNccp

            liste_xNccp = [0] + liste_xNccp
            liste_xNccp = liste_xNccp + [lenbN]

            liste_yNccp = [0] + liste_yNccp
            liste_yNccp = liste_yNccp + [0]

            liste_stringNccp = liste_yNccp

            fig = plt.figure()
            
            fig, ax = plt.subplots() ###

            plt.figure(figsize=(round(lenbN/120),round(lenbN/75)))
            plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, komplementarni vlakno', fontsize=round(lenbN/65))
            plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbN/65))
            #plt.title('Souhrn ORFu s CC strukturou, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxNccp, dyNccp in zip(liste_xNccp, liste_yNccp):
                if (sxNccp == 0):
                    plt.annotate(' '+str(dyNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95))
                    plt.annotate('     -  '+str(sxNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95), rotation=90)
                elif (sxNccp > 0) and (sxNccp < lenbN):
                    for sxNccpX, dyNccpX in zip(liste_xNccp[1:-1], liste_yNccp[1:-1]):
                        plt.annotate(''+str(dyNccpX), xy = (sxNccpX, dyNccpX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifNccX = liste_xNccp[1:-1].index(sxNccpX)
                        plt.annotate('     - ORF:'+str(sxNccpX - 3*ind_MrkcN_coords[ifNccX] - 3)+' - CC:'+str(sxNccpX), xy = (sxNccpX, dyNccpX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxNccp == lenbN):
                    plt.annotate(' '+str(dyNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95))
                    plt.annotate('     -  '+str(sxNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95), rotation=90)
                else:
                    None
                #plt.scatter(liste_xNccp, liste_stringNccp, s = lenbN/1.5, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xNccp, liste_stringNccp, width = 10, align='center', color='black') # width = lenbN/120
                    
            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbN/80))
            plt.yticks(size=round(lenbN/80))

            #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno.svg')
            plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno.png', dpi=32) # , dpi=16
            plt.close()

        # dava souradnice ORFu a CC struktury
        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,komplementarni vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,komplementarni vlakno, vsechny ramce - zacatek:
        # dava souradnice ORFu a CC struktury
        # pro BIG verzi - roztazeni obr. pres celou str.
        
        if (len(CC_all_N_coords) == 0):
            print('Zadny ORF s CC sekundarni strukturou v komplementarnim vlaknu.')
            print('\n')
        else:
            xNccp = []
            yNccp = []

            for iNccp in range(len(CC_all_N_coords)):

                xNccp += [max_results_CC_N[iNccp][0] + 3*ind_MrkcN_coords[iNccp] + 3]
                yNccp += [round(max_results_CC_N[iNccp][1]*1000)/10]

            liste_xNccp = xNccp
            liste_yNccp = yNccp

            liste_xNccp = [0] + liste_xNccp
            liste_xNccp = liste_xNccp + [lenbN]

            liste_yNccp = [0] + liste_yNccp
            liste_yNccp = liste_yNccp + [0]

            liste_stringNccp = liste_yNccp

            fig = plt.figure()
            
            fig, ax = plt.subplots() ###

            plt.figure(figsize=(round(lenbN/40),round(lenbN/75)))
            plt.xlabel('Pozice maxima CC struktury (v bp) \n\n Souhrn ORFu s CC strukturou, komplementarni vlakno', fontsize=round(lenbN/65))
            plt.ylabel('Pravdepodobnost CC struktury (%)', fontsize=round(lenbN/65))
            #plt.title('Souhrn ORFu s CC strukturou, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxNccp, dyNccp in zip(liste_xNccp, liste_yNccp):
                if (sxNccp == 0):
                    plt.annotate(' '+str(dyNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95))
                    plt.annotate('     -  '+str(sxNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95), rotation=90)
                elif (sxNccp > 0) and (sxNccp < lenbN):
                    for sxNccpX, dyNccpX in zip(liste_xNccp[1:-1], liste_yNccp[1:-1]):
                        plt.annotate(''+str(dyNccpX), xy = (sxNccpX, dyNccpX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifNccX = liste_xNccp[1:-1].index(sxNccpX)
                        plt.annotate('     - ORF:'+str(sxNccpX - 3*ind_MrkcN_coords[ifNccX] - 3)+' - CC:'+str(sxNccpX), xy = (sxNccpX, dyNccpX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxNccp == lenbN):
                    plt.annotate(' '+str(dyNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95))
                    plt.annotate('     -  '+str(sxNccp), xy = (sxNccp, dyNccp), fontsize=round(lenbN/95), rotation=90)
                else:
                    None
                #plt.scatter(liste_xNccp, liste_stringNccp, s = lenbN/1.5, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xNccp, liste_stringNccp, width = 10, align='center', color='black') # width = lenbN/120
                    
            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbN/80))
            plt.yticks(size=round(lenbN/80))

            #plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno-BIG.svg')
            plt.savefig('figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno-BIG.png', dpi=32) # , dpi=16
            plt.close()

        # pro BIG verzi - roztazeni obr. pres celou str.
        # dava souradnice ORFu a CC struktury
        # pro obr. souhrnu ORFu s CC strukturou, PRAVDEPODOBNOSTI ,komplementarni vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - zacatek:
        # upravene - ukaze ORFy a pozici zacatku TM struktury

        if (len(listTMsP_pozice) == 0):
            print('Zadny ORF s TM sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xPtmA = []
            yPtmA = []

            for iPtmA in range(len(listTMsP_pozice)):

                xPtmA += [count_M_P[iPtmA][0] + 3*apM_P_coords[iPtmA][0] + 3]
                yPtmA += [count_M_P[iPtmA][1]]

            liste_xPtmA = xPtmA
            liste_yPtmA = yPtmA
            
            liste_xPtmA = [0] + liste_xPtmA
            liste_xPtmA = liste_xPtmA + [lenbP]
            
            liste_yPtmA = [0] + liste_yPtmA
            liste_yPtmA = liste_yPtmA + [0]
            
            liste_stringPtmA = liste_yPtmA

            fig = plt.figure()
            fig, ax = plt.subplots()

            #ax.set_xlim([0,lenbP])
            #fig.add_axes([0,0,lenbP,lenbP])
            
            plt.figure(figsize=(round(lenbP/120),round(lenbP/75)))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxPtmA, dyPtmA in zip(liste_xPtmA, liste_yPtmA):
             
                #plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                if (sxPtmA == 0):
                    plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                elif (sxPtmA > 0) and (sxPtmA < lenbP):
                    for sxPtmAX, dyPtmAX in zip(liste_xPtmA[1:-1], liste_yPtmA[1:-1]):
                        plt.annotate(''+str(dyPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifPX = liste_xPtmA[1:-1].index(sxPtmAX)
                        plt.annotate('     - ORF:'+str(sxPtmAX - 3*apM_P_coords[ifPX][0] - 3)+' - TM:'+str(sxPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxPtmA == lenbP):
                    plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                    
                #plt.scatter(liste_xPtmA, liste_stringPtmA, s = lenbP/1.8, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xPtmA, liste_stringPtmA, width = 10, align='center', color='black') # width = lenbP/250,

            plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, zadane vlakno', fontsize=round(lenbP/65))
            plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbP/65))
            #plt.title('Souhrn ORFu s TM strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno.svg')
            plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno.png', dpi=32) # , dpi=16
            plt.close()

        # upravene - ukaze ORFy a pozici zacatku TM struktury
        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - zacatek:
        # upravene - ukaze ORFy a pozici zacatku TM struktury
        # pro BIG verzi - roztazeni obr. pres celou str.

        if (len(listTMsP_pozice) == 0):
            print('Zadny ORF s TM sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xPtmA = []
            yPtmA = []

            for iPtmA in range(len(listTMsP_pozice)):

                xPtmA += [count_M_P[iPtmA][0] + 3*apM_P_coords[iPtmA][0] + 3]
                yPtmA += [count_M_P[iPtmA][1]]

            liste_xPtmA = xPtmA
            liste_yPtmA = yPtmA
            
            liste_xPtmA = [0] + liste_xPtmA
            liste_xPtmA = liste_xPtmA + [lenbP]
            
            liste_yPtmA = [0] + liste_yPtmA
            liste_yPtmA = liste_yPtmA + [0]
            
            liste_stringPtmA = liste_yPtmA

            fig = plt.figure()
            fig, ax = plt.subplots()

            #ax.set_xlim([0,lenbP])
            #fig.add_axes([0,0,lenbP,lenbP])
            
            plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxPtmA, dyPtmA in zip(liste_xPtmA, liste_yPtmA):
             
                #plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                if (sxPtmA == 0):
                    plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                elif (sxPtmA > 0) and (sxPtmA < lenbP):
                    for sxPtmAX, dyPtmAX in zip(liste_xPtmA[1:-1], liste_yPtmA[1:-1]):
                        plt.annotate(''+str(dyPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifPX = liste_xPtmA[1:-1].index(sxPtmAX)
                        plt.annotate('     - ORF:'+str(sxPtmAX - 3*apM_P_coords[ifPX][0] - 3)+' - TM:'+str(sxPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxPtmA == lenbP):
                    plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                    
                #plt.scatter(liste_xPtmA, liste_stringPtmA, s = lenbP/1.8, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xPtmA, liste_stringPtmA, width = 10, align='center', color='black') # width = lenbP/250,

            plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, zadane vlakno', fontsize=round(lenbP/65))
            plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbP/65))
            #plt.title('Souhrn ORFu s TM strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG.svg')
            plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG.png', dpi=32) # , dpi=16
            plt.close()

        # pro BIG verzi - roztazeni obr. pres celou str.
        # upravene - ukaze ORFy a pozici zacatku TM struktury
        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - zacatek:
        # upravene - ukaze ORFy a pozici zacatku TM struktury
        # pro BIG verzi - roztazeni obr. pres celou str.
        # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro TM

        if (len(listTMsP_pozice) == 0):
            print('Zadny ORF s TM sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:

            ## lidsky lokus
            ## popis: >hg38_dna range=chr19:17389648-17406217 5'pad=0 3'pad=0 strand=+
            #TMs = [[705,774],[12101,12170]]

            ## kureci lokus
            ## popis: >galGal6_dna range=chr28:3531163-3538110 5'pad=0 3'pad=0 strand=-
            #TMs = [[1990,2059],[3303,3372]]

            ## mysi lokus
            ## popis: >mm10_dna range=chr8:71519771-71538962 5'pad=0 3'pad=0 strand=-
            #TMs = [[1619,1688],[13293,13362]]

            ## luskoun ostrovni
            ## popis:
            ## >ref|NW_023435982.1|:679097-704212 Manis javanica isolate MJ74 unplaced genomic scaffold, YNU_ManJav_2.0 scaffold_88, whole genome shotgun sequence
            #TMs = [[5864,5933],[15151,15220]]

            ## kalon vabivy
            ## popis:
            ## >ref|NW_006429864.1|:923553-939337 Pteropus alecto unplaced genomic scaffold, ASM32557v1 scaffold160, whole genome shotgun sequence
            TMs = [[4047,4116],[10555,10624]]
             
            xPtmA = []
            yPtmA = []

            for iPtmA in range(len(listTMsP_pozice)):

                xPtmA += [count_M_P[iPtmA][0] + 3*apM_P_coords[iPtmA][0] + 3]
                yPtmA += [count_M_P[iPtmA][1]]

            liste_xPtmA = xPtmA
            liste_yPtmA = yPtmA
            
            liste_xPtmA = [0] + liste_xPtmA
            liste_xPtmA = liste_xPtmA + [lenbP]
            
            liste_yPtmA = [0] + liste_yPtmA
            liste_yPtmA = liste_yPtmA + [0]
            
            liste_stringPtmA = liste_yPtmA

            fig = plt.figure()
            fig, ax = plt.subplots()

            #ax.set_xlim([0,lenbP])
            #fig.add_axes([0,0,lenbP,lenbP])
            
            plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxPtmA, dyPtmA in zip(liste_xPtmA, liste_yPtmA):
             
                #plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                if (sxPtmA == 0):
                    plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                elif (sxPtmA > 0) and (sxPtmA < lenbP):
                    for sxPtmAX, dyPtmAX in zip(liste_xPtmA[1:-1], liste_yPtmA[1:-1]):
                        plt.annotate(''+str(dyPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifPX = liste_xPtmA[1:-1].index(sxPtmAX)
                        plt.annotate('     - ORF:'+str(sxPtmAX - 3*apM_P_coords[ifPX][0] - 3)+' - TM:'+str(sxPtmAX), xy = (sxPtmAX, dyPtmAX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxPtmA == lenbP):
                    plt.annotate(''+str(dyPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPtmA), xy = (sxPtmA, dyPtmA), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                    
                #plt.scatter(liste_xPtmA, liste_stringPtmA, s = lenbP/1.8, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xPtmA, liste_stringPtmA, width = 10, align='center', color='black') # width = lenbP/250,

            for iTMs in range(len(TMs)):
                plt.hlines(y=0,xmin=(TMs[iTMs][0]),xmax=(TMs[iTMs][1]),color='red',linestyle='-',lw=1200)

            plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, zadane vlakno', fontsize=round(lenbP/65))
            plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbP/65))
            #plt.title('Souhrn ORFu s TM strukturou, zadane vlakno\n\n', fontsize=round(lenbP/65))

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG-L.svg')
            plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG-L.png', dpi=32) # , dpi=16
            plt.close()

        # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro TM
        # pro BIG verzi - roztazeni obr. pres celou str.
        # upravene - ukaze ORFy a pozici zacatku TM struktury
        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,zadane vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,komplementarni vlakno, vsechny ramce - zacatek:
        # upravene - ukaze ORFy a pozici zacatku TM struktury

        if (len(listTMsN_pozice) == 0):
            print('Zadny ORF s TM sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xNtmA = []
            yNtmA = []

            for iNtmA in range(len(listTMsN_pozice)):

                xNtmA += [count_M_N[iNtmA][0] + 3*apM_N_coords[iNtmA][0] + 3]
                yNtmA += [count_M_N[iNtmA][1]]

            liste_xNtmA = xNtmA
            liste_yNtmA = yNtmA
            
            liste_xNtmA = [0] + liste_xNtmA
            liste_xNtmA = liste_xNtmA + [lenbN]
            
            liste_yNtmA = [0] + liste_yNtmA
            liste_yNtmA = liste_yNtmA + [0]
            
            liste_stringNtmA = liste_yNtmA

            fig = plt.figure()
            fig, ax = plt.subplots()

            #ax.set_xlim([0,lenbN])
            #fig.add_axes([0,0,lenbN,lenbN])
            
            plt.figure(figsize=(round(lenbN/120),round(lenbN/75)))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxNtmA, dyNtmA in zip(liste_xNtmA, liste_yNtmA):
             
                #plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                if (sxNtmA == 0):
                    plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                elif (sxNtmA > 0) and (sxNtmA < lenbN):
                    for sxNtmAX, dyNtmAX in zip(liste_xNtmA[1:-1], liste_yNtmA[1:-1]):
                        plt.annotate(''+str(dyNtmAX), xy = (sxNtmAX, dyNtmAX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifNX = liste_xNtmA[1:-1].index(sxNtmAX)
                        plt.annotate('     - ORF:'+str(sxNtmAX - 3*apM_N_coords[ifNX][0] - 3)+' - TM:'+str(sxNtmAX), xy = (sxNtmAX, dyNtmAX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxNtmA == lenbN):
                    plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                    
                #plt.scatter(liste_xNtmA, liste_stringNtmA, s = lenbN/1.8, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xNtmA, liste_stringNtmA, width = 10, align='center', color='black') # width = lenbN/250,

            plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, komplementarni vlakno', fontsize=round(lenbN/65))
            plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbN/65))
            #plt.title('Souhrn ORFu s TM strukturou, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbN/80))
            plt.yticks(size=round(lenbN/80))

            #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno.svg')
            plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno.png', dpi=32) # , dpi=16
            plt.close()

        # upravene - ukaze ORFy a pozici zacatku TM struktury
        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,komplementarni vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,komplementarni vlakno, vsechny ramce - zacatek:
        # upravene - ukaze ORFy a pozici zacatku TM struktury
        # pro BIG verzi - roztazeni obr. pres celou str.
        
        if (len(listTMsN_pozice) == 0):
            print('Zadny ORF s TM sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xNtmA = []
            yNtmA = []

            for iNtmA in range(len(listTMsN_pozice)):

                xNtmA += [count_M_N[iNtmA][0] + 3*apM_N_coords[iNtmA][0] + 3]
                yNtmA += [count_M_N[iNtmA][1]]

            liste_xNtmA = xNtmA
            liste_yNtmA = yNtmA
            
            liste_xNtmA = [0] + liste_xNtmA
            liste_xNtmA = liste_xNtmA + [lenbN]
            
            liste_yNtmA = [0] + liste_yNtmA
            liste_yNtmA = liste_yNtmA + [0]
            
            liste_stringNtmA = liste_yNtmA

            fig = plt.figure()
            fig, ax = plt.subplots()

            #ax.set_xlim([0,lenbN])
            #fig.add_axes([0,0,lenbN,lenbN])
            
            plt.figure(figsize=(round(lenbN/40),round(lenbN/75)))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxNtmA, dyNtmA in zip(liste_xNtmA, liste_yNtmA):
             
                #plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')

                if (sxNtmA == 0):
                    plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                elif (sxNtmA > 0) and (sxNtmA < lenbN):
                    for sxNtmAX, dyNtmAX in zip(liste_xNtmA[1:-1], liste_yNtmA[1:-1]):
                        plt.annotate(''+str(dyNtmAX), xy = (sxNtmAX, dyNtmAX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifNX = liste_xNtmA[1:-1].index(sxNtmAX)
                        plt.annotate('     - ORF:'+str(sxNtmAX - 3*apM_N_coords[ifNX][0] - 3)+' - TM:'+str(sxNtmAX), xy = (sxNtmAX, dyNtmAX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxNtmA == lenbN):
                    plt.annotate(''+str(dyNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxNtmA), xy = (sxNtmA, dyNtmA), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                    
                #plt.scatter(liste_xNtmA, liste_stringNtmA, s = lenbN/1.8, c ='black') # picker=(True & 10000) # yellow
                plt.bar(liste_xNtmA, liste_stringNtmA, width = 10, align='center', color='black') # width = lenbN/250,

            plt.xlabel('Pozice pocatku TM (v bp) \n\n Souhrn ORFu s TM strukturou, komplementarni vlakno', fontsize=round(lenbN/65))
            plt.ylabel('Pocet AMK v TM strukture (-)', fontsize=round(lenbN/65))
            #plt.title('Souhrn ORFu s TM strukturou, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbN/80))
            plt.yticks(size=round(lenbN/80))

            #plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno-BIG.svg')
            plt.savefig('figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno-BIG.png', dpi=32) # , dpi=16
            plt.close()

        # pro BIG verzi - roztazeni obr. pres celou str.
        # upravene - ukaze ORFy a pozici zacatku TM struktury
        # pro obr. souhrnu ORFu s TM strukturou, pocet AMK ,komplementarni vlakno, vsechny ramce - konec.

        # PRIPRAVA pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane i komplementarni vlakno, vsechny ramce - zacatek:

        import os
        import subprocess

        owd = os.getcwd() # owd = original working directory
        
        os.chdir("/home/oem/Documents/HTMLproTMaCCaGPIrelatall/kohgpi-1.5/")
        #os.system("./kohgpi t")              # neni vubec treba
        #os.system("./mapfill2A >kohgpi.map") # neni vubec treba
        os.system("kohgpi i")
        process_GPI_N = subprocess.Popen(['/home/oem/Documents/HTMLproTMaCCaGPIrelatall/kohgpi-1.5/kohgpi','/home/oem/Documents/HTMLproTMaCCaGPIrelatall/FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa'],
                                   stdout = subprocess.PIPE,
                                   stderr = subprocess.PIPE)
        stdout, stderr = process_GPI_N.communicate()
        stdout

        os.chdir("/home/oem/Documents/HTMLproTMaCCaGPIrelatall/kohgpi-1.5/")
        #os.system("./kohgpi t")              # neni vubec treba
        #os.system("./mapfill2A >kohgpi.map") # neni vubec treba
        os.system("kohgpi i")
        process_GPI_N = subprocess.Popen(['/home/oem/Documents/HTMLproTMaCCaGPIrelatall/kohgpi-1.5/kohgpi','/home/oem/Documents/HTMLproTMaCCaGPIrelatall/FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa'],
                                   stdout = subprocess.PIPE,
                                   stderr = subprocess.PIPE)
        stdout, stderr = process_GPI_N.communicate()
        stdout

        os.chdir(owd) # owd = original working directory

        ###

        from Bio import SeqIO

        if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos')) == False:
            print('Neni *.pos a ...# xcoords.txt soubor, zadane vlakno.')
            GPI_xcoords_P_R = []

        else:

            pos_filename_P = "FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos"

            xcoords_filename_P = "FASTA_output/peptidesP/peptidesP_c00_Standard_Code.xcoords.txt"

            input_handle_P = open(pos_filename_P,"r")

            output_handle_P = open(xcoords_filename_P,"w")

            for seq_record in SeqIO.parse(input_handle_P,"fasta"):
                output_handle_P.write("%s\n" % (
                    seq_record.id))

            output_handle_P.close()
            input_handle_P.close()

            GPI_xcoords_P = open('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.xcoords.txt')
            GPI_xcoords_P_r = GPI_xcoords_P.read()

            count_GPI_P_n = GPI_xcoords_P_r.count('\n')
            GPI_xcoords_P_rr = GPI_xcoords_P_r.split('\n',count_GPI_P_n)
            GPI_xcoords_P_rrr = GPI_xcoords_P_rr[0:-1]


            iGPiL = []

            for iGP in GPI_xcoords_P_rrr:
                iGPi = int(iGP)
                iGPiL += [iGPi]

            GPI_xcoords_P_R = iGPiL

        #print('GPI_xcoords_P_R je: ', GPI_xcoords_P_R)

        ###

        if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos')) == False:
            print('Neni *.pos a ...xcoords.txt soubor, komplementarni vlakno.')
            GPI_xcoords_N_R = []
        else:

            pos_filename_N = "FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos"

            xcoords_filename_N = "FASTA_output/peptidesN/peptidesN_c00_Standard_Code.xcoords.txt"

            input_handle_N = open(pos_filename_N,"r")

            output_handle_N = open(xcoords_filename_N,"w")

            for seq_record in SeqIO.parse(input_handle_N,"fasta"):
                output_handle_N.write("%s\n" % (
                    seq_record.id))

            output_handle_N.close()
            input_handle_N.close()

            GPI_xcoords_N = open('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.xcoords.txt')
            GPI_xcoords_N_r = GPI_xcoords_N.read()

            count_GPI_N_n = GPI_xcoords_N_r.count('\n')
            GPI_xcoords_N_rr = GPI_xcoords_N_r.split('\n',count_GPI_N_n)
            GPI_xcoords_N_rrr = GPI_xcoords_N_rr[0:-1]


            iGNiL = []

            for iGN in GPI_xcoords_N_rrr:
                iGNi = int(iGN)
                iGNiL += [iGNi]

            GPI_xcoords_N_R = iGNiL

        ###

        if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos')) == False:
            print('Neni *.pos a ...ycoords.txt soubor, zadane vlakno.')
        else:

            posy_filename_P = "FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos"

            ycoords_filename_P = "FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt"

            input_handle_Py = open(posy_filename_P,"r")

            output_handle_Py = open(ycoords_filename_P,"w")

            for seq_record in SeqIO.parse(input_handle_Py,"fasta"):
                output_handle_Py.write("%s\n" % (
                    seq_record.description))

            output_handle_Py.close()
            input_handle_Py.close()

        ###

        if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos')) == False:
            print('Neni *.pos a ...ycoords.txt soubor, komplementarni vlakno.')
        else:

            posy_filename_N = "FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos"

            ycoords_filename_N = "FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt"

            input_handle_Ny = open(posy_filename_N,"r")

            output_handle_Ny = open(ycoords_filename_N,"w")

            for seq_record in SeqIO.parse(input_handle_Ny,"fasta"):
                output_handle_Ny.write("%s\n" % (
                    seq_record.description))

            output_handle_Ny.close()
            input_handle_Ny.close()

        ###

        import re

        if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt')) == False:
            print('Neni ....ycoords.txt soubor pro ziskani % kvality omega-mista, zadane vlakno.')
            omega_sum_P = []
        else:

            OY_P = open("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt",'r')
            OYR_P = OY_P.readlines()

            omega_sum_P = []

            for line_P in OYR_P:
                
                ma_P = re.match(r"(\w+) (\w+) (\W+\w+) (\w+) (\w+\W+\w+) (\w+) (\w+)",line_P)

                if ma_P == None:
                    omegaprocento_P = 0
                else:
                    omegaprocento10_P = 10*int(ma_P[7][0])
                    omegaprocento1_P = 1*int(ma_P[7][1])
                    omegaprocento_P = omegaprocento10_P + omegaprocento1_P
                omega_sum_P += [omegaprocento_P]

        #

        if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt')) == False:
            print('Neni ....ycoords.txt soubor pro ziskani % kvality  omega-mista, komplementarni vlakno.')
            omega_sum_N = []
        else:

            OY_N = open("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt",'r')
            OYR_N = OY_N.readlines()

            omega_sum_N = []

            for line_N in OYR_N:
                
                ma_N = re.match(r"(\w+) (\w+) (\W+\w+) (\w+) (\w+\W+\w+) (\w+) (\w+)",line_N)

                if ma_N == None:
                    omegaprocento_N = 0
                else:
                    omegaprocento10_N = 10*int(ma_N[7][0])
                    omegaprocento1_N = 1*int(ma_N[7][1])
                    omegaprocento_N = omegaprocento10_N + omegaprocento1_N
                omega_sum_N += [omegaprocento_N]

        # PRIPRAVA pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane i komplementarni vlakno, vsechny ramce - konec.

        # vytvareni mnozin - zacatek:

        # vytvareni mnozin TM, CC, GPI - zadane vlakno:

        set_listTMsP_pozice2_0 = set(listTMsP_pozice2_0)

        set_CC_all_P_coords_0 = set(CC_all_P_coords_0)

        set_GPI_xcoords_P_R = set(GPI_xcoords_P_R)


        # vytvareni mnozin s jen TM, nebo jen s CC, nebo jen s GPI - zadane vlakno:

        jenTMgpi_P = set_listTMsP_pozice2_0.difference(set_CC_all_P_coords_0)
        set_jenTM_P = jenTMgpi_P.difference(set_GPI_xcoords_P_R)
        list_jenTM_P = list(set_jenTM_P) # vysledek - tj. seznam pozic, zacatku ORFu s jen TM strukturou, v zadanem vlaknu

        jenCCgpi_P = set_CC_all_P_coords_0.difference(set_listTMsP_pozice2_0)
        set_jenCC_P = jenCCgpi_P.difference(set_GPI_xcoords_P_R)
        list_jenCC_P = list(set_jenCC_P) # vysledek - tj. seznam pozic, zacatku ORFu s jen CC strukturou, v zadanem vlaknu

        jenGPIcc_P = set_GPI_xcoords_P_R.difference(set_listTMsP_pozice2_0)
        set_jenGPI_P = jenGPIcc_P.difference(set_CC_all_P_coords_0)
        list_jenGPI_P = list(set_jenGPI_P) # vysledek - tj. seznam pozic, zacatku ORFu s jen GPI strukturou, v zadanem vlaknu
        
        # vytvareni mnozin TM, CC, GPI - komplementarni vlakno:

        set_listTMsN_pozice2_0 = set(listTMsN_pozice2_0)

        set_CC_all_N_coords_0 = set(CC_all_N_coords_0)

        set_GPI_xcoords_N_R = set(GPI_xcoords_N_R)


        # vytvareni mnozin s jen TM, nebo jen s CC, nebo jen s GPI - komplementarni vlakno:

        jenTMgpi_N = set_listTMsN_pozice2_0.difference(set_CC_all_N_coords_0)
        set_jenTM_N = jenTMgpi_N.difference(set_GPI_xcoords_N_R)
        list_jenTM_N = list(set_jenTM_N) # vysledek - tj. seznam pozic, zacatku ORFu s jen TM strukturou, v komplementarnim vlaknu


        jenCCgpi_N = set_CC_all_N_coords_0.difference(set_listTMsN_pozice2_0)
        set_jenCC_N = jenCCgpi_N.difference(set_GPI_xcoords_N_R)
        list_jenCC_N = list(set_jenCC_N) # vysledek - tj. seznam pozic, zacatku ORFu s jen CC strukturou, v komplementarnim vlaknu

        jenGPIcc_N = set_GPI_xcoords_N_R.difference(set_listTMsN_pozice2_0)
        set_jenGPI_N = jenGPIcc_N.difference(set_CC_all_N_coords_0)
        list_jenGPI_N = list(set_jenGPI_N) # vysledek - tj. seznam pozic, zacatku ORFu s jen GPI strukturou, v komplementarnim vlaknu


        # vytvareni pruniku vzdy jen dvou mnozin - zadane vlakno:
        # napr. kdyz chci TM a CC, musim udelat jejich prunik, ale odecist od tohoto pruniku, spolecne prvky s GPI

        set_TMaCCsGPI_P = set_listTMsP_pozice2_0.intersection(set_CC_all_P_coords_0)
        set_TMaCC_P = set_TMaCCsGPI_P.difference(set_GPI_xcoords_P_R)
        list_TMaCC_P = list(set_TMaCC_P)

        set_TMaGPIsCC_P = set_listTMsP_pozice2_0.intersection(set_GPI_xcoords_P_R)
        set_TMaGPI_P = set_TMaGPIsCC_P.difference(set_CC_all_P_coords_0)
        list_TMaGPI_P = list(set_TMaGPI_P)

        set_GPIaCCsTM_P = set_GPI_xcoords_P_R.intersection(set_CC_all_P_coords_0)
        set_GPIaCC_P = set_GPIaCCsTM_P.difference(set_listTMsP_pozice2_0)
        list_GPIaCC_P = list(set_GPIaCC_P)


        # vytvareni pruniku vzdy jen dvou mnozin - komplementarni vlakno:
        # napr. kdyz chci TM a CC, musim udelat jejich prunik, ale odecist od tohoto pruniku, spolecne prvky s GPI

        set_TMaCCsGPI_N = set_listTMsN_pozice2_0.intersection(set_CC_all_N_coords_0)
        set_TMaCC_N = set_TMaCCsGPI_N.difference(set_GPI_xcoords_N_R)
        list_TMaCC_N = list(set_TMaCC_N)

        set_TMaGPIsCC_N = set_listTMsN_pozice2_0.intersection(set_GPI_xcoords_N_R)
        set_TMaGPI_N = set_TMaGPIsCC_N.difference(set_CC_all_N_coords_0)
        list_TMaGPI_N = list(set_TMaGPI_N)

        set_GPIaCCsTM_N = set_GPI_xcoords_N_R.intersection(set_CC_all_N_coords_0)
        set_GPIaCC_N = set_GPIaCCsTM_N.difference(set_listTMsN_pozice2_0)
        list_GPIaCC_N = list(set_GPIaCC_N)

        
        # vytvareni pruniku vsech 3 (TM,CC,GPI) mnozin:
        # zadane vlakno:
        set_TMaCCxx_P = set_listTMsP_pozice2_0.intersection(set_CC_all_P_coords_0)
        set_TMaCCaGPI_P = set_GPI_xcoords_P_R.intersection(set_TMaCCxx_P)
        list_TMaCCaGPI_P = list(set_TMaCCaGPI_P)

        # komplementarni vlakno:
        set_TMaCCxx_N = set_listTMsN_pozice2_0.intersection(set_CC_all_N_coords_0)
        set_TMaCCaGPI_N = set_GPI_xcoords_N_R.intersection(set_TMaCCxx_N)
        list_TMaCCaGPI_N = list(set_TMaCCaGPI_N)
            
        # vytvareni mnozin - konec.

        # ziskavani obou souradnic pro GPI strukturu - zadane i komplementarni vlakno - zacatek:

        GPI_XYcoords_P_R = []

        for iLOP in listORFsP_pozice:
            for iGPP in GPI_xcoords_P_R:
                if iGPP == iLOP[0]:
                    GPI_XYcoords_P_R += [iLOP]
                else:
                    continue

        GPI_XYcoords_N_R = []

        for iLON in listORFsN_pozice:
            for iGPN in GPI_xcoords_N_R:
                if iGPN == iLON[0]:
                    GPI_XYcoords_N_R += [iLON]
                else:
                    continue

        # ziskavani obou souradnic pro GPI strukturu - zadane i komplementarni vlakno - konec.

        # ziskavani souradnic pro ORFy jen s TM, zadane i komplementarni vlakno - zacatek:

        list_jenTM_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLOP in listORFsP_pozice:
            for iTPP in list_jenTM_P:
                if iTPP == iLOP[0]:
                    list_jenTM_P_01 += [iLOP]
                else:
                    continue
        # vysledek: list_jenTM_P_01

        list_jenTM_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLON in listORFsN_pozice:
            for iTPN in list_jenTM_N:
                if iTPN == iLON[0]:
                    list_jenTM_N_01 += [iLON]
                else:
                    continue
        # vysledek: list_jenTM_N_01

        # ziskavani souradnic pro ORFy jen s TM, zadane i komplementarni vlakno - konec.

        # ziskavani souradnic pro ORFy jen s CC, zadane i komplementarni vlakno - zacatek:

        list_jenCC_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLOP in listORFsP_pozice:
            for iCPP in list_jenCC_P:
                if iCPP == iLOP[0]:
                    list_jenCC_P_01 += [iLOP]
                else:
                    continue
        # vysledek: list_jenCC_P_01

        list_jenCC_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLON in listORFsN_pozice:
            for iCPN in list_jenCC_N:
                if iCPN == iLON[0]:
                    list_jenCC_N_01 += [iLON]
                else:
                    continue
        # vysledek: list_jenCC_N_01

        # ziskavani souradnic pro ORFy jen s CC, zadane i komplementarni vlakno - konec.

        # ziskavani souradnic pro ORFy jen s GPI, zadane i komplementarni vlakno - zacatek:

        list_jenGPI_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLOP in listORFsP_pozice:
            for iGPP in list_jenGPI_P:
                if iGPP == iLOP[0]:
                    list_jenGPI_P_01 += [iLOP]
                else:
                    continue
        # vysledek: list_jenGPI_P_01

        list_jenGPI_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLON in listORFsN_pozice:
            for iGPN in list_jenGPI_N:
                if iGPN == iLON[0]:
                    list_jenGPI_N_01 += [iLON]
                else:
                    continue
        # vysledek: list_jenGPI_N_01

        # ziskavani souradnic pro ORFy jen s GPI, zadane i komplementarni vlakno - konec.

        # ziskavani souradnic pro ORFy jen s TM a zaroven s CC, zadane i komplementarni vlakno - zacatek:

        list_TMaCC_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLOP in listORFsP_pozice:
            for iTCPP in list_TMaCC_P:
                if iTCPP == iLOP[0]:
                    list_TMaCC_P_01 += [iLOP]
                else:
                    continue
        # vysledek: list_TMaCC_P_01

        list_TMaCC_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLON in listORFsN_pozice:
            for iTCPN in list_TMaCC_N:
                if iTCPN == iLON[0]:
                    list_TMaCC_N_01 += [iLON]
                else:
                    continue
        # vysledek: list_TMaCC_N_01

        # ziskavani souradnic pro ORFy jen s TM a zaroven s CC, zadane i komplementarni vlakno - konec.
        
        # ziskavani souradnic pro ORFy jen s TM a zaroven s GPI, zadane i komplementarni vlakno - zacatek:

        list_TMaGPI_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLOP in listORFsP_pozice:
            for iTGPP in list_TMaGPI_P:
                if iTGPP == iLOP[0]:
                    list_TMaGPI_P_01 += [iLOP]
                else:
                    continue
        # vysledek: list_TMaGPI_P_01

        list_TMaGPI_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLON in listORFsN_pozice:
            for iTGPN in list_TMaGPI_N:
                if iTGPN == iLON[0]:
                    list_TMaGPI_N_01 += [iLON]
                else:
                    continue
        # vysledek: list_TMaGPI_N_01

        # ziskavani souradnic pro ORFy jen s TM a zaroven s GPI, zadane i komplementarni vlakno - konec.

        # ziskavani souradnic pro ORFy jen s CC a zaroven s GPI, zadane i komplementarni vlakno - zacatek:

        list_GPIaCC_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLOP in listORFsP_pozice:
            for iGCPP in list_GPIaCC_P:
                if iGCPP == iLOP[0]:
                    list_GPIaCC_P_01 += [iLOP]
                else:
                    continue
        # vysledek: list_GPIaCC_P_01

        list_GPIaCC_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLON in listORFsN_pozice:
            for iGCPN in list_GPIaCC_N:
                if iGCPN == iLON[0]:
                    list_GPIaCC_N_01 += [iLON]
                else:
                    continue
        # vysledek: list_GPIaCC_N_01

        # ziskavani souradnic pro ORFy jen s CC a zaroven s GPI, zadane i komplementarni vlakno - konec.

        # ziskavani souradnic pro ORFy s TM a zaroven s CC a zaroven s GPI, zadane i komplementarni vlakno - zacatek:

        list_TMaCCaGPI_P_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLOP in listORFsP_pozice:
            for iTCGPP in list_TMaCCaGPI_P:
                if iTCGPP == iLOP[0]:
                    list_TMaCCaGPI_P_01 += [iLOP]
                else:
                    continue
        # vysledek: list_TMaCCaGPI_P_01

        list_TMaCCaGPI_N_01 = [] # vzdy 0.a 1. pozice v ORFu dohromady, jako 1 polozka v listu=seznamu

        for iLON in listORFsN_pozice:
            for iTCGPN in list_TMaCCaGPI_N:
                if iTCGPN == iLON[0]:
                    list_TMaCCaGPI_N_01 += [iLON]
                else:
                    continue
        # vysledek: list_TMaCCaGPI_N_01

        # ziskavani souradnic pro ORFy s TM a zaroven s CC a zaroven s GPI, zadane i komplementarni vlakno - konec.

        # "linkove" obrazky pro CC (P i N) - zacatek:

        for ilPcc in range(len(CC_all_P_coords)):

            # "linkovy obr.":
            xsP = np.linspace(1,1,lenbP)
            plt.figure(figsize=(lenbP/120,lenbP/4800))
            plt.hlines(y=1,xmin=1,xmax=lenbP,color='blue',linestyle='-',lw=20)
            plt.hlines(y=1,xmin=(CC_all_P_coords[ilPcc][0]),xmax=(CC_all_P_coords[ilPcc][1]),color='r',linestyle='-',lw=30)
            plt.yticks([])
            plt.xticks([])
            '''
            if ((listORFsP_pozice[irpP_c00_4fA][0])%3 == 0):
                frameP_ORFall_lin = 1
            elif ((listORFsP_pozice[irpP_c00_4fA][0] + 1)%3 == 0):
                frameP_ORFall_lin = 3
            elif ((listORFsP_pozice[irpP_c00_4fA][0] + 2)%3 == 0):
                frameP_ORFall_lin = 2
            else:
                continue
            '''
            #plt.annotate(str(listORFsP_pozice[irpP_c00_4fA][0]),(listORFsP_pozice[irpP_c00_4fA][0],1),fontsize=16)
            #plt.annotate(str(listORFsP_pozice[irpP_c00_4fA][1]),(listORFsP_pozice[irpP_c00_4fA][1],1),fontsize=16)
            #plt.annotate(str(lenbP),(lenbP,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
            #plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
            '''
            plt.annotate('Fr.:'+str(frameP_ORFall_lin),((listORFsP_pozice[irpP_c00_4fA][0]+listORFsP_pozice[irpP_c00_4fA][1])/2,1),fontsize=16,color='green',verticalalignment='top')
            '''
            #plt.savefig('/home/oem/Documents/_html/figures/NA string - fig No. ' + str([irpP_c00_4fA]) + ' provided_string.png')
            plt.savefig('figures/given/NA string - fig No. ' + str(CC_all_P_coords[ilPcc]) + ' provided_string.png')
            #plt.show()
            plt.close()#

        for ilNcc in range(len(CC_all_N_coords)):

            # "linkovy obr.":
            xsN = np.linspace(1,1,lenbN)
            plt.figure(figsize=(lenbN/120,lenbN/4800))
            plt.hlines(y=1,xmin=1,xmax=lenbN,color='blue',linestyle='-',lw=20)
            plt.hlines(y=1,xmin=(CC_all_N_coords[ilNcc][0]),xmax=(CC_all_N_coords[ilNcc][1]),color='r',linestyle='-',lw=30)
            plt.yticks([])
            plt.xticks([])
            '''
            if ((listORFsN_pozice[irpN_c00_4fA][0])%3 == 0):
                frameN_ORFall_lin = 1
            elif ((listORFsN_pozice[irpN_c00_4fA][0] + 1)%3 == 0):
                frameN_ORFall_lin = 3
            elif ((listORFsN_pozice[irpN_c00_4fA][0] + 2)%3 == 0):
                frameN_ORFall_lin = 2
            else:
                continue
            '''
            #plt.annotate(str(listORFsN_pozice[irpN_c00_4fA][0]),(listORFsN_pozice[irpN_c00_4fA][0],1),fontsize=16)
            #plt.annotate(str(listORFsN_pozice[irpN_c00_4fA][1]),(listORFsN_pozice[irpN_c00_4fA][1],1),fontsize=16)
            #plt.annotate(str(lenbN),(lenbN,1),fontsize=12) # pridan popisek delky vlakna na konec modre cary v obr.
            #plt.annotate(str(1),(1,1),fontsize=12) #  pridan popisek zacatku vlakna na zacatek modre cary v obr. = 1. bp, (cislo bp.=1)
            '''
            plt.annotate('Fr.:'+str(frameN_ORFall_lin),((listORFsN_pozice[irpN_c00_4fA][0]+listORFsN_pozice[irpN_c00_4fA][1])/2,1),fontsize=16,color='green',verticalalignment='top')
            '''
            #plt.savefig('/home/oem/Documents/_html/figures/NA string - fig No. ' + str([irpN_c00_4fA]) + ' complementary_string.png')
            plt.savefig('figures/complementary/NA string - fig No. ' + str(CC_all_N_coords[ilNcc]) + ' complementary_string.png')
            #plt.show()
            plt.close()#

        # "linkove" obrazky pro CC (P i N) - konec.

        # "linkove" obrazky pro GPI (P i N) - zacatek:

        for ilPgpi in range(len(GPI_XYcoords_P_R)):

            # "linkovy obr.":
            xsP = np.linspace(1,1,lenbP)
            plt.figure(figsize=(lenbP/120,lenbP/4800))
            plt.hlines(y=1,xmin=1,xmax=lenbP,color='blue',linestyle='-',lw=20)
            plt.hlines(y=1,xmin=(GPI_XYcoords_P_R[ilPgpi][0]),xmax=(GPI_XYcoords_P_R[ilPgpi][1]),color='r',linestyle='-',lw=30)
            
            plt.yticks([])
            plt.xticks([])

            plt.savefig('figures/given/NA string - fig No. ' + str(GPI_XYcoords_P_R[ilPgpi]) + ' provided_string.png')
            #plt.show()
            plt.close()

        for ilNgpi in range(len(GPI_XYcoords_N_R)):

            # "linkovy obr.":
            xsN = np.linspace(1,1,lenbN)
            plt.figure(figsize=(lenbN/120,lenbN/4800))
            plt.hlines(y=1,xmin=1,xmax=lenbN,color='blue',linestyle='-',lw=20)
            plt.hlines(y=1,xmin=(GPI_XYcoords_N_R[ilNgpi][0]),xmax=(GPI_XYcoords_N_R[ilNgpi][1]),color='r',linestyle='-',lw=30)
            
            plt.yticks([])
            plt.xticks([])

            plt.savefig('figures/complementary/NA string - fig No. ' + str(GPI_XYcoords_N_R[ilNgpi]) + ' complementary_string.png')
            #plt.show()
            plt.close()

        # "linkove" obrazky pro GPI (P i N) - konec.

        # vytvareni slovniku pro vykresleni popisu GPI mista v tabulce v HTML - zacatek: - NEPOUZITO v HTML

        if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos')) == False:
            print('Neni *.pos soubor pro tvorbu slovniku pro GPI strukturu, zadane vlakno. - stejne je to NEpouzito v HTML')
        else:
            record_dict_P = SeqIO.to_dict(SeqIO.parse("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.faa.pos","fasta"))

        if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos')) == False:
            print('Neni *.pos soubor pro tvorbu slovniku pro GPI strukturu, komplementarni vlakno.- stejne je to NEpouzito v HTML')
        else:
            record_dict_N = SeqIO.to_dict(SeqIO.parse("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.faa.pos","fasta"))

        # vytvareni slovniku pro vykresleni popisu GPI mista v tabulce v HTML - konec. - NEPOUZITO v HTML

        # vytvareni obrazku pro GPI sestrihove misto - zadane vlakno, vsechny ramce - zacatek:

        import re
        from Bio import SeqIO

        if (os.path.isfile('FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt')) == False:
            print('Neni ...ycoords.txt soubor pro tvorbu obrazku pro GPI strukturu, zadane vlakno, vsechny ramce.')
        else:

            oGPIy_P = open("FASTA_output/peptidesP/peptidesP_c00_Standard_Code.ycoords.txt",'r')
            GPIy_P = oGPIy_P.readlines()

            distances_AAs_P = []
            cleavages_P = []
            XYcoords_P = []

            for line_gpi_P in GPIy_P:

                ma_XYcoords_P = re.match(r"(\w+) (\w+)",line_gpi_P).groups()
                ma_Xcoord_P = int(ma_XYcoords_P[0])
                ma_Ycoord_P = int(ma_XYcoords_P[1])

                XYcoords_P += [[ma_Xcoord_P, ma_Ycoord_P]]

                distance_AAs_P = round((ma_Ycoord_P - ma_Xcoord_P)/3)

                distances_AAs_P += [distance_AAs_P]

                ma_Pc = re.match(r"(\w+) (\w+) (\W+\w+) (\w+) (\w+\W+\w+) (\w+) (\w+)",line_gpi_P)

                if ma_Pc == None:
                    cleavage_int_P = distance_AAs_P
                    cleavages_P += [cleavage_int_P]
                else:
                    cleavage_P = ma_Pc[5][2:]
                    cleavage_int_P = round(int(cleavage_P))
                    cleavages_P += [cleavage_int_P]
                        
            import matplotlib.pyplot as plt

            for iGPIobrP in range(len(omega_sum_P)):

                fig = plt.figure()
                fig, ax = plt.subplots()
                #plt.figure(figsize=(lenbP/200,lenbP/256))

                ax.set_xlim([0,distances_AAs_P[iGPIobrP]])
                ax.set_ylim([0,100])

                if ( distances_AAs_P[iGPIobrP] == cleavages_P[iGPIobrP] ):

                    plt.bar(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], omega_sum_P[iGPIobrP], width = 1, align='center', color='green')
                    plt.annotate(omega_sum_P[iGPIobrP], xy=(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], omega_sum_P[iGPIobrP]) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'bottom')
                    #plt.annotate('None', xy=(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                    plt.annotate('Neurceno omega-misto', xy=( (distances_AAs_P[iGPIobrP])/2,0 ) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                    #plt.annotate('None', xy=(distances_AAs_P[iGPIobrP],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                    
                else:
                    
                    plt.bar(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], omega_sum_P[iGPIobrP], width = 1, align='center', color='green')
                    plt.annotate(omega_sum_P[iGPIobrP], xy=(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], omega_sum_P[iGPIobrP]) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP], xy=(distances_AAs_P[iGPIobrP]-cleavages_P[iGPIobrP],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')

                plt.xlabel('Poradi AMK (N-C konec) a omega-misto (-)', fontsize=16)
                plt.ylabel('Kvalita omega-mista (%)', fontsize=16)
                plt.title('ORF ' + str(XYcoords_P[iGPIobrP]) +' s GPI strukturou, zadane vlakno.\n', fontsize=14)

                plt.savefig('figures/GPI/P' + str(XYcoords_P[iGPIobrP]) + '.png')

                #plt.show()
                plt.close()

        # vytvareni obrazku pro GPI sestrihove misto - zadane vlakno, vsechny ramce - konec.

        # vytvareni obrazku pro GPI sestrihove misto - komplementarni vlakno, vsechny ramce - zacatek:

        #import re
        #from Bio import SeqIO

        if (os.path.isfile('FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt')) == False:
            print('Neni ...ycoords.txt soubor pro tvorbu obrazku pro GPI strukturu, komplementarni vlakno, vsechny ramce.')
        else:

            oGPIy_N = open("FASTA_output/peptidesN/peptidesN_c00_Standard_Code.ycoords.txt",'r')
            GPIy_N = oGPIy_N.readlines()

            distances_AAs_N = []
            cleavages_N = []
            XYcoords_N = []

            for line_gpi_N in GPIy_N:

                ma_XYcoords_N = re.match(r"(\w+) (\w+)",line_gpi_N).groups()
                ma_Xcoord_N = int(ma_XYcoords_N[0])
                ma_Ycoord_N = int(ma_XYcoords_N[1])

                XYcoords_N += [[ma_Xcoord_N, ma_Ycoord_N]]

                distance_AAs_N = round((ma_Ycoord_N - ma_Xcoord_N)/3)

                distances_AAs_N += [distance_AAs_N]

                ma_Nc = re.match(r"(\w+) (\w+) (\W+\w+) (\w+) (\w+\W+\w+) (\w+) (\w+)",line_gpi_N)

                if ma_Nc == None:
                    cleavage_int_N = distance_AAs_N
                    cleavages_N += [cleavage_int_N]
                else:
                    cleavage_N = ma_Nc[5][2:]
                    cleavage_int_N = round(int(cleavage_N))
                    cleavages_N += [cleavage_int_N]
                        
            #import matplotlib.pyplot as plt

            for iGPIobrN in range(len(omega_sum_N)):

                fig = plt.figure()
                fig, ax = plt.subplots()
                #plt.figure(figsize=(lenbN/200,lenbN/256))

                ax.set_xlim([0,distances_AAs_N[iGPIobrN]])
                ax.set_ylim([0,100])

                if ( distances_AAs_N[iGPIobrN] == cleavages_N[iGPIobrN] ):

                    plt.bar(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], omega_sum_N[iGPIobrN], width = 1, align='center', color='green')
                    plt.annotate(omega_sum_N[iGPIobrN], xy=(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], omega_sum_N[iGPIobrN]) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'bottom')
                    #plt.annotate('None', xy=(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                    plt.annotate('Neurceno omega-misto', xy=( (distances_AAs_N[iGPIobrN])/2,0 ) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                    #plt.annotate('None', xy=(distances_AAs_N[iGPIobrN],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')
                    
                else:
                    
                    plt.bar(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], omega_sum_N[iGPIobrN], width = 1, align='center', color='green')
                    plt.annotate(omega_sum_N[iGPIobrN], xy=(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], omega_sum_N[iGPIobrN]) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN], xy=(distances_AAs_N[iGPIobrN]-cleavages_N[iGPIobrN],0) ,fontsize=12, horizontalalignment = 'center', verticalalignment = 'top', color='red')

                plt.xlabel('Poradi AMK (N-C konec) a omega-misto (-)', fontsize=16)
                plt.ylabel('Kvalita omega-mista (%)', fontsize=16)
                plt.title('ORF ' + str(XYcoords_N[iGPIobrN]) +' s GPI strukturou, zadane vlakno.\n', fontsize=14)

                plt.savefig('figures/GPI/N/' + str(XYcoords_N[iGPIobrN]) + '.png')

                #plt.show()
                plt.close()

        # vytvareni obrazku pro GPI sestrihove misto - komplementarni vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - zacatek:
        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist

        if (len(omega_sum_P) == 0):
            print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xPgpi = []
            for ixGPI_P in range(len(GPI_xcoords_P_R)):
                xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P])]
                '''
                if (distances_AAs_P[ixGPI_P] == cleavages_P[ixGPI_P]):
                    xPgpi += [GPI_xcoords_P_R[ixGPI_P]]
                elif (distances_AAs_P[ixGPI_P] > cleavages_P[ixGPI_P]):
                    xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P]) + 3]
                else:
                    None
                '''
            xPgpi
            yPgpi = omega_sum_P
            
            liste_xPgpi = xPgpi
            liste_yPgpi = yPgpi

            liste_xPgpi = [0] + liste_xPgpi
            liste_xPgpi = liste_xPgpi + [lenbP]

            liste_yPgpi = [0] + liste_yPgpi
            liste_yPgpi = liste_yPgpi + [0]

            liste_stringPgpi = liste_yPgpi

            fig = plt.figure()
            
            fig, ax = plt.subplots()

            plt.figure(figsize=(round(lenbP/120),round(lenbP/75)))
            plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, zadane vlakno', fontsize=round(lenbP/65))
            plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbP/65))
            #plt.title('Souhrn ORFu s omega-misty, zadane vlakno\n\n', fontsize=round(lenbP/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxPgpi, dyPgpi in zip(liste_xPgpi, liste_yPgpi):
             
                if (sxPgpi == 0):
                    plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                elif ((sxPgpi > 0) and (sxPgpi < lenbP)):
                    for sxPgpiX, dyPgpiX in zip(liste_xPgpi[1:-1], liste_yPgpi[1:-1]):
                        plt.annotate(''+str(dyPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifPgpiX = liste_xPgpi[1:-1].index(sxPgpiX)
                        plt.annotate('     - ORF:'+str(sxPgpiX - 3*(distances_AAs_P[ifPgpiX] - cleavages_P[ifPgpiX]))+' - GPI:'+str(sxPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxPgpi == lenbP):
                    plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                
                #plt.scatter(liste_xPgpi, liste_stringPgpi, s = lenbP/1.8, c ='black')
                plt.bar(liste_xPgpi, liste_stringPgpi, width = 10,  align='center', color='black') # width = lenbP/250, 

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno.svg')
            plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno.png', dpi=32) # , dpi=16
            plt.close()

        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - zacatek:
        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
        # pro BIG verzi - roztazeni obr. pres celou str.
        
        if (len(omega_sum_P) == 0):
            print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xPgpi = []
            for ixGPI_P in range(len(GPI_xcoords_P_R)):
                xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P])]
                '''
                if (distances_AAs_P[ixGPI_P] == cleavages_P[ixGPI_P]):
                    xPgpi += [GPI_xcoords_P_R[ixGPI_P]]
                elif (distances_AAs_P[ixGPI_P] > cleavages_P[ixGPI_P]):
                    xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P]) + 3]
                else:
                    None
                '''
            xPgpi
            yPgpi = omega_sum_P
            
            liste_xPgpi = xPgpi
            liste_yPgpi = yPgpi

            liste_xPgpi = [0] + liste_xPgpi
            liste_xPgpi = liste_xPgpi + [lenbP]

            liste_yPgpi = [0] + liste_yPgpi
            liste_yPgpi = liste_yPgpi + [0]

            liste_stringPgpi = liste_yPgpi

            fig = plt.figure()
            
            fig, ax = plt.subplots()

            plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))
            plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, zadane vlakno', fontsize=round(lenbP/65))
            plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbP/65))
            #plt.title('Souhrn ORFu s omega-misty, zadane vlakno\n\n', fontsize=round(lenbP/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxPgpi, dyPgpi in zip(liste_xPgpi, liste_yPgpi):
             
                if (sxPgpi == 0):
                    plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                elif ((sxPgpi > 0) and (sxPgpi < lenbP)):
                    for sxPgpiX, dyPgpiX in zip(liste_xPgpi[1:-1], liste_yPgpi[1:-1]):
                        plt.annotate(''+str(dyPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifPgpiX = liste_xPgpi[1:-1].index(sxPgpiX)
                        plt.annotate('     - ORF:'+str(sxPgpiX - 3*(distances_AAs_P[ifPgpiX] - cleavages_P[ifPgpiX]))+' - GPI:'+str(sxPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxPgpi == lenbP):
                    plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                
                #plt.scatter(liste_xPgpi, liste_stringPgpi, s = lenbP/1.8, c ='black')
                plt.bar(liste_xPgpi, liste_stringPgpi, width = 10,  align='center', color='black') # width = lenbP/250, 

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG.svg')
            plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG.png', dpi=32) # , dpi=16
            plt.close()

        # pro BIG verzi - roztazeni obr. pres celou str.
        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - zacatek:
        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
        # pro BIG verzi - roztazeni obr. pres celou str.
        # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro GPI
        
        if (len(omega_sum_P) == 0):
            print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:

            ## lidsky lokus
            ## popis: >hg38_dna range=chr19:17389648-17406217 5'pad=0 3'pad=0 strand=+
            #GPI = 2463

            ## kureci lokus
            ## popis: >galGal6_dna range=chr28:3531163-3538110 5'pad=0 3'pad=0 strand=-
            #GPI = 2817

            ## mysi lokus
            ## popis: >mm10_dna range=chr8:71519771-71538962 5'pad=0 3'pad=0 strand=-
            #GPI = 4210

            ## luskoun ostrovni
            ## popis:
            ## >ref|NW_023435982.1|:679097-704212 Manis javanica isolate MJ74 unplaced genomic scaffold, YNU_ManJav_2.0 scaffold_88, whole genome shotgun sequence
            #GPI = 7606
            
            ## kalon vabivy
            ## popis:
            ## >ref|NW_006429864.1|:923553-939337 Pteropus alecto unplaced genomic scaffold, ASM32557v1 scaffold160, whole genome shotgun sequence
            GPI = 5720

            xPgpi = []
            for ixGPI_P in range(len(GPI_xcoords_P_R)):
                xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P])]
                '''
                if (distances_AAs_P[ixGPI_P] == cleavages_P[ixGPI_P]):
                    xPgpi += [GPI_xcoords_P_R[ixGPI_P]]
                elif (distances_AAs_P[ixGPI_P] > cleavages_P[ixGPI_P]):
                    xPgpi += [GPI_xcoords_P_R[ixGPI_P] + 3*(distances_AAs_P[ixGPI_P] - cleavages_P[ixGPI_P]) + 3]
                else:
                    None
                '''
            xPgpi
            yPgpi = omega_sum_P
            
            liste_xPgpi = xPgpi
            liste_yPgpi = yPgpi

            liste_xPgpi = [0] + liste_xPgpi
            liste_xPgpi = liste_xPgpi + [lenbP]

            liste_yPgpi = [0] + liste_yPgpi
            liste_yPgpi = liste_yPgpi + [0]

            liste_stringPgpi = liste_yPgpi

            fig = plt.figure()
            
            fig, ax = plt.subplots()

            plt.figure(figsize=(round(lenbP/40),round(lenbP/75)))
            plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, zadane vlakno', fontsize=round(lenbP/65))
            plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbP/65))
            #plt.title('Souhrn ORFu s omega-misty, zadane vlakno\n\n', fontsize=round(lenbP/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxPgpi, dyPgpi in zip(liste_xPgpi, liste_yPgpi):
             
                if (sxPgpi == 0):
                    plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                elif ((sxPgpi > 0) and (sxPgpi < lenbP)):
                    for sxPgpiX, dyPgpiX in zip(liste_xPgpi[1:-1], liste_yPgpi[1:-1]):
                        plt.annotate(''+str(dyPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifPgpiX = liste_xPgpi[1:-1].index(sxPgpiX)
                        plt.annotate('     - ORF:'+str(sxPgpiX - 3*(distances_AAs_P[ifPgpiX] - cleavages_P[ifPgpiX]))+' - GPI:'+str(sxPgpiX), xy = (sxPgpiX, dyPgpiX), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxPgpi == lenbP):
                    plt.annotate(''+str(dyPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxPgpi), xy = (sxPgpi, dyPgpi), fontsize=round(lenbP/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                
                #plt.scatter(liste_xPgpi, liste_stringPgpi, s = lenbP/1.8, c ='black')
                plt.bar(liste_xPgpi, liste_stringPgpi, width = 10,  align='center', color='black') # width = lenbP/250, 

            wth = 400 # tj. width
            hht = 4 # tj. height
            xy = (GPI-0.5*wth,0-0.5*hht)
            GPI_rectangle = plt.Rectangle(xy,wth,hht,angle=0.0,color='yellow')
            fig = plt.gcf()
            ax = fig.gca()
            ax.add_patch(GPI_rectangle)

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbP/80))
            plt.yticks(size=round(lenbP/80))

            #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG-L.svg')
            plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG-L.png', dpi=32) # , dpi=16
            plt.close()

        # pro kap. 5, Vysledky - "linkovy obr." do celkoveho obr. pro GPI
        # pro BIG verzi - roztazeni obr. pres celou str.
        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,zadane vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,komplementarni vlakno, vsechny ramce - zacatek:
        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist

        if (len(omega_sum_N) == 0):
            print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xNgpi = []
            for ixGPI_N in range(len(GPI_xcoords_N_R)):
                xNgpi += [GPI_xcoords_N_R[ixGPI_N] + 3*(distances_AAs_N[ixGPI_N] - cleavages_N[ixGPI_N])]
                '''
                if (distances_AAs_N[ixGPI_N] == cleavages_N[ixGPI_N]):
                    xNgpi += [GPI_xcoords_N_R[ixGPI_N]]
                elif (distances_AAs_N[ixGPI_N] > cleavages_N[ixGPI_N]):
                    xNgpi += [GPI_xcoords_N_R[ixGPI_N] + 3*(distances_AAs_N[ixGPI_N] - cleavages_N[ixGPI_N]) + 3]
                else:
                    None
                '''
            xNgpi
            yNgpi = omega_sum_N
            
            liste_xNgpi = xNgpi
            liste_yNgpi = yNgpi

            liste_xNgpi = [0] + liste_xNgpi
            liste_xNgpi = liste_xNgpi + [lenbN]

            liste_yNgpi = [0] + liste_yNgpi
            liste_yNgpi = liste_yNgpi + [0]

            liste_stringNgpi = liste_yNgpi

            fig = plt.figure()
            
            fig, ax = plt.subplots()

            plt.figure(figsize=(round(lenbN/120),round(lenbN/75)))
            plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, komplementarni vlakno', fontsize=round(lenbN/65))
            plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbN/65))
            #plt.title('Souhrn ORFu s omega-misty, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxNgpi, dyNgpi in zip(liste_xNgpi, liste_yNgpi):
             
                if (sxNgpi == 0):
                    plt.annotate(''+str(dyNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                elif ((sxNgpi > 0) and (sxNgpi < lenbN)):
                    for sxNgpiX, dyNgpiX in zip(liste_xNgpi[1:-1], liste_yNgpi[1:-1]):
                        plt.annotate(''+str(dyNgpiX), xy = (sxNgpiX, dyNgpiX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifNgpiX = liste_xNgpi[1:-1].index(sxNgpiX)
                        plt.annotate('     - ORF:'+str(sxNgpiX - 3*(distances_AAs_N[ifNgpiX] - cleavages_N[ifNgpiX]))+' - GPI:'+str(sxNgpiX), xy = (sxNgpiX, dyNgpiX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxNgpi == lenbN):
                    plt.annotate(''+str(dyNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                
                #plt.scatter(liste_xNgpi, liste_stringNgpi, s = lenbN/1.8, c ='black')
                plt.bar(liste_xNgpi, liste_stringNgpi, width = 10,  align='center', color='black') # width = lenbN/250, 

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbN/80))
            plt.yticks(size=round(lenbN/80))

            #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno.svg')
            plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno.png', dpi=32) # , dpi=16
            plt.close()

        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,komplementarni vlakno, vsechny ramce - konec.

        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,komplementarni vlakno, vsechny ramce - zacatek:
        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
        # pro BIG verzi - roztazeni obr. pres celou str.
        
        if (len(omega_sum_N) == 0):
            print('Zadny ORF s GPI sekundarni strukturou v zadanem vlaknu.')
            print('\n')
        else:
            xNgpi = []
            for ixGPI_N in range(len(GPI_xcoords_N_R)):
                xNgpi += [GPI_xcoords_N_R[ixGPI_N] + 3*(distances_AAs_N[ixGPI_N] - cleavages_N[ixGPI_N])]
                '''
                if (distances_AAs_N[ixGPI_N] == cleavages_N[ixGPI_N]):
                    xNgpi += [GPI_xcoords_N_R[ixGPI_N]]
                elif (distances_AAs_N[ixGPI_N] > cleavages_N[ixGPI_N]):
                    xNgpi += [GPI_xcoords_N_R[ixGPI_N] + 3*(distances_AAs_N[ixGPI_N] - cleavages_N[ixGPI_N]) + 3]
                else:
                    None
                '''
            xNgpi
            yNgpi = omega_sum_N
            
            liste_xNgpi = xNgpi
            liste_yNgpi = yNgpi

            liste_xNgpi = [0] + liste_xNgpi
            liste_xNgpi = liste_xNgpi + [lenbN]

            liste_yNgpi = [0] + liste_yNgpi
            liste_yNgpi = liste_yNgpi + [0]

            liste_stringNgpi = liste_yNgpi

            fig = plt.figure()
            
            fig, ax = plt.subplots()

            plt.figure(figsize=(round(lenbN/40),round(lenbN/75)))
            plt.xlabel('Pozice omega-mista (v bp) \n\n Souhrn ORFu s omega-misty, komplementarni vlakno', fontsize=round(lenbN/65))
            plt.ylabel('Kvalita omega-mista (%)', fontsize=round(lenbN/65))
            #plt.title('Souhrn ORFu s omega-misty, komplementarni vlakno\n\n', fontsize=round(lenbN/65))

            #ax = plt.gca()
            #fig = plt.gcf()

            for sxNgpi, dyNgpi in zip(liste_xNgpi, liste_yNgpi):
             
                if (sxNgpi == 0):
                    plt.annotate(''+str(dyNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                elif ((sxNgpi > 0) and (sxNgpi < lenbN)):
                    for sxNgpiX, dyNgpiX in zip(liste_xNgpi[1:-1], liste_yNgpi[1:-1]):
                        plt.annotate(''+str(dyNgpiX), xy = (sxNgpiX, dyNgpiX), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                        ifNgpiX = liste_xNgpi[1:-1].index(sxNgpiX)
                        plt.annotate('     - ORF:'+str(sxNgpiX - 3*(distances_AAs_N[ifNgpiX] - cleavages_N[ifNgpiX]))+' - GPI:'+str(sxNgpiX), xy = (sxNgpiX, dyNgpiX), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                    break
                elif (sxNgpi == lenbN):
                    plt.annotate(''+str(dyNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), horizontalalignment = 'center', verticalalignment = 'bottom')
                    plt.annotate('     -  '+str(sxNgpi), xy = (sxNgpi, dyNgpi), fontsize=round(lenbN/95), rotation=90, horizontalalignment = 'center')
                else:
                    None
                
                #plt.scatter(liste_xNgpi, liste_stringNgpi, s = lenbN/1.8, c ='black')
                plt.bar(liste_xNgpi, liste_stringNgpi, width = 10,  align='center', color='black') # width = lenbN/250, 

            plt.xticks(rotation=90)
            plt.xticks(size=round(lenbN/80))
            plt.yticks(size=round(lenbN/80))

            #plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno-BIG.svg')
            plt.savefig('figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno-BIG.png', dpi=32) # , dpi=16
            plt.close()

        # pro BIG verzi - roztazeni obr. pres celou str.
        # zahrnuje souradnice ORFu i souradnice sestrihovych omega-mist
        # pro obr. souhrnu ORFu s GPI strukturou, kvalita omega mista (v %) ,komplementarni vlakno, vsechny ramce - konec.

        
        index23 = open("index23.html","w")

        from airium import Airium

        a = Airium()

        a('<!DOCTYPE html>')
        with a.html(lang='en'):
          with a.head():
              a.meta(charset="utf-8")
              a.title(_t="TM, CC a GPI predikce")
              a.link(href='styly.css', rel='stylesheet')

          with a.body():
              with a.h3(id="id17"):
                a("TM, CC a GPI predikce:")
              with a.p():
                a('Skript použitelný pro: standardní kód, kód archeí a rostlinných plastidů a kvasinkový alternativní jaderný kód.')
              with a.p():
                a('dolní limit délky vstupní sekvence je (bp):')
                with a.b():a(dolni_limit)
              with a.p():
                  a('horní limit délky vstupní sekvence je (bp):')
                  with a.b():a(horni_limit)
              with a.p():
                a('nukleotidů (délka zadané NK) je (bp):')
                with a.b():a(lenbP)
              with a.p():
                a('minimální délka ORFu je (bp):')
                with a.b():a(limitORF)
              with a.p():
                a('počet ORFů v zadaném vláknu je:')
                with a.b():a(len(listORFsP))
              with a.p():
                a('počet ORFů v komplementárním vláknu je:')
                with a.b():a(len(listORFsN))
              with a.p():
                a('celkový počet ORFů v obou vláknech je:')
                with a.b():a(len(listORFsP) + len(listORFsN))
              with a.p():
                a('počet ORFů s TM sekundární strukturou (zadané vlákno): ')
                with a.b():a(count_TM_dvoj_obrazku_P)
              with a.p():
                a('počet ORFů s TM sekundární strukturou (komplementární vlákno): ')
                with a.b():a(count_TM_dvoj_obrazku_N)
              with a.p():
                a('počet ORFů s CC sekundární strukturou (zadané vlákno): ')
                with a.b():a(CC_count_P)
              with a.p():
                a('počet ORFů s CC sekundární strukturou (komplementární vlákno): ')
                with a.b():a(CC_count_N)
              with a.p():
                a('práh pro detekci CC sekundární struktury - zadané vlákno (%): ')
                with a.b():a(CC_threshold_P*100)
              with a.p():
                a('práh pro detekci CC sekundární struktury - komplementární vlákno (%): ')
                with a.b():a(CC_threshold_N*100)
              with a.p():
                a('počet ORFů s GPI sekundární strukturou (zadané vlákno): ')
                with a.b():a(len(GPI_xcoords_P_R))
              with a.p():
                a('počet ORFů s GPI sekundární strukturou (komplementární vlákno): ')
                with a.b():a(len(GPI_xcoords_N_R))

              with a.p():
                a('počet ORFů jen s TM sekundární strukturou (zadané vlákno): ')
                with a.b():a(len(list_jenTM_P))
              with a.p():
                a('počet ORFů jen s TM sekundární strukturou (komplementární vlákno): ')
                with a.b():a(len(list_jenTM_N))

              with a.p():
                a('počet ORFů jen s CC sekundární strukturou (zadané vlákno): ')
                with a.b():a(len(list_jenCC_P))
              with a.p():
                a('počet ORFů jen s CC sekundární strukturou (komplementární vlákno): ')
                with a.b():a(len(list_jenCC_N))

              with a.p():
                a('počet ORFů jen s GPI sekundární strukturou (zadané vlákno): ')
                with a.b():a(len(list_jenGPI_P))
              with a.p():
                a('počet ORFů jen s GPI sekundární strukturou (komplementární vlákno): ')
                with a.b():a(len(list_jenGPI_N))
                
              with a.p():
                a('počet ORFů jen s TM a CC sekundární strukturou (zadané vlákno): ')
                with a.b():a(len(list_TMaCC_P))
              with a.p():
                a('počet ORFů jen s TM a CC sekundární strukturou (komplementární vlákno): ')
                with a.b():a(len(list_TMaCC_N))

              with a.p():
                a('počet ORFů jen s TM a GPI sekundární strukturou (zadané vlákno): ')
                with a.b():a(len(list_TMaGPI_P))
              with a.p():
                a('počet ORFů jen s TM a GPI sekundární strukturou (komplementární vlákno): ')
                with a.b():a(len(list_TMaGPI_N))

              with a.p():
                a('počet ORFů jen s CC a GPI sekundární strukturou (zadané vlákno): ')
                with a.b():a(len(list_GPIaCC_P))
              with a.p():
                a('počet ORFů jen s CC a GPI sekundární strukturou (komplementární vlákno): ')
                with a.b():a(len(list_GPIaCC_N))

              with a.p():
                a('počet ORFů s TM, CC a GPI sekundární strukturou současně (zadané vlákno): ')
                with a.b():a(len(list_TMaCCaGPI_P))
              with a.p():
                a('počet ORFů s TM, CC a GPI sekundární strukturou současně (komplementární vlákno): ')
                with a.b():a(len(list_TMaCCaGPI_N))

              #with a.p(align='center'):a.a(href='#(konec stránky)', _t='(přejít na konec stránky)')

              a.hr()
              a.hr()
              with a.h2(id="id8", klass='main_header'):a('ORFy - zadané vlákno DNA:')
              with a.table():
                  with a.tr():
                      with a.td():a.img(src="figures/NA-ORFs-Pframe1-given_string-2D.png", alt="")
                      with a.td():a.img(src="figures/NA-ORFs-Pframe2-given_string-2D.png", alt="")
                      with a.td():a.img(src="figures/NA-ORFs-Pframe3-given_string-2D.png", alt="")
                  with a.tr():
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno.png", alt="")
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno.png", alt="")
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno.png", alt="")

              with a.table():
                  with a.tr():
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_TM_strukturou,AMK,zadane_vlakno-BIG.png", alt="")
                  with a.tr():
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,zadane_vlakno-BIG.png", alt="")
                  with a.tr():
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_GPI_strukturou,zadane_vlakno-BIG.png", alt="")

                      
              with a.table(cellpadding='5', border='2', bordercolor='black'): #
                   with a.tr():
                      with a.th(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s TM strukturou - zadané vlákno')
                      with a.th(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s CC strukturou - zadané vlákno')
                      with a.th(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s GPI strukturou - zadané vlákno')
                   with a.tr(): #
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('1. čtecí rámec:')
                              for itmP1all in listTMsP_pozice:
                                      if ((itmP1all[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmP1all), _t='ORF-pozice: '+str(itmP1all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('1. čtecí rámec:')
                              for iccP1all in CC_all_P_coords:
                                      if ((iccP1all[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccP1all), _t='ORF-pozice: '+str(iccP1all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('1. čtecí rámec:')
                              for igpiP1all in GPI_XYcoords_P_R:
                                      if ((igpiP1all[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiP1all), _t='ORF-pozice: '+str(igpiP1all)+'\n')
                                      else:
                                          None
                   with a.tr(): ##
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('2. čtecí rámec:')
                              for itmP2all in listTMsP_pozice:
                                      if ((itmP2all[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmP2all), _t='ORF-pozice: '+str(itmP2all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('2. čtecí rámec:')
                              for iccP2all in CC_all_P_coords:
                                      if ((iccP2all[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccP2all), _t='ORF-pozice: '+str(iccP2all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('2. čtecí rámec:')
                              for igpiP2all in GPI_XYcoords_P_R:
                                      if ((igpiP2all[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiP2all), _t='ORF-pozice: '+str(igpiP2all)+'\n')
                                      else:
                                          None
                   with a.tr(): ###
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('3. čtecí rámec:')
                              for itmP3all in listTMsP_pozice:
                                      if ((itmP3all[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmP3all), _t='ORF-pozice: '+str(itmP3all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('3. čtecí rámec:')
                              for iccP3all in CC_all_P_coords:
                                      if ((iccP3all[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccP3all), _t='ORF-pozice: '+str(iccP3all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('3. čtecí rámec:')
                              for igpiP3all in GPI_XYcoords_P_R:
                                      if ((igpiP3all[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiP3all), _t='ORF-pozice: '+str(igpiP3all)+'\n')
                                      else:
                                          None
              with a.table(cellpadding='5', border='2', bordercolor='black'): # width='99%'
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM strukturou - zadané vlákno - 1. čtecí rámec:')
                              for itmP1 in list_jenTM_P_01:
                                  if ((itmP1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmP1), _t='ORF-pozice: '+str(itmP1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM strukturou - zadané vlákno - 2. čtecí rámec:')
                              for itmP2 in list_jenTM_P_01:
                                  if ((itmP2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmP2), _t='ORF-pozice: '+str(itmP2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM strukturou - zadané vlákno - 3. čtecí rámec:')
                              for itmP3 in list_jenTM_P_01:
                                  if ((itmP3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmP3), _t='ORF-pozice: '+str(itmP3)+'\n')
                                  else:
                                      None # # #
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC strukturou - zadané vlákno - 1. čtecí rámec:')
                              for iccP1 in list_jenCC_P_01:
                                  if ((iccP1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccP1), _t='ORF-pozice: '+str(iccP1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC strukturou - zadané vlákno - 2. čtecí rámec:')
                              for iccP2 in list_jenCC_P_01:
                                  if ((iccP2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccP2), _t='ORF-pozice: '+str(iccP2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC strukturou - zadané vlákno - 3. čtecí rámec:')
                              for iccP3 in list_jenCC_P_01:
                                  if ((iccP3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccP3), _t='ORF-pozice: '+str(iccP3)+'\n')
                                  else:
                                      None # # # 
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s GPI strukturou - zadané vlákno - 1. čtecí rámec:')
                              for igpiP1 in list_jenGPI_P_01:
                                  if ((igpiP1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(igpiP1), _t='ORF-pozice: '+str(igpiP1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s GPI strukturou - zadané vlákno - 2. čtecí rámec:')
                              for igpiP2 in list_jenGPI_P_01:
                                  if ((igpiP2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(igpiP2), _t='ORF-pozice: '+str(igpiP2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s GPI strukturou - zadané vlákno - 3. čtecí rámec:')
                              for igpiP3 in list_jenGPI_P_01:
                                  if ((igpiP3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(igpiP3), _t='ORF-pozice: '+str(igpiP3)+'\n')
                                  else:
                                      None
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a CC strukturou - zadané vlákno - 1. čtecí rámec:')
                              for itmccP1 in list_TMaCC_P_01:
                                  if ((itmccP1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccP1), _t='ORF-pozice: '+str(itmccP1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a CC strukturou - zadané vlákno - 2. čtecí rámec:')
                              for itmccP2 in list_TMaCC_P_01:
                                  if ((itmccP2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccP2), _t='ORF-pozice: '+str(itmccP2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a CC strukturou - zadané vlákno - 3. čtecí rámec:')
                              for itmccP3 in list_TMaCC_P_01:
                                  if ((itmccP3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccP3), _t='ORF-pozice: '+str(itmccP3)+'\n')
                                  else:
                                      None # # #
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a GPI strukturou - zadané vlákno - 1. čtecí rámec:')
                              for itmgpiP1 in list_TMaGPI_P_01:
                                  if ((itmgpiP1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmgpiP1), _t='ORF-pozice: '+str(itmgpiP1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a GPI strukturou - zadané vlákno - 2. čtecí rámec:')
                              for itmgpiP2 in list_TMaGPI_P_01:
                                  if ((itmgpiP2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmgpiP2), _t='ORF-pozice: '+str(itmgpiP2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a GPI strukturou - zadané vlákno - 3. čtecí rámec:')
                              for itmgpiP3 in list_TMaGPI_P_01:
                                  if ((itmgpiP3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmgpiP3), _t='ORF-pozice: '+str(itmgpiP3)+'\n')
                                  else:
                                      None # # #
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC a GPI strukturou - zadané vlákno - 1. čtecí rámec:')
                              for iccgpiP1 in list_GPIaCC_P_01:
                                  if ((iccgpiP1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccgpiP1), _t='ORF-pozice: '+str(iccgpiP1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC a GPI strukturou - zadané vlákno - 2. čtecí rámec:')
                              for iccgpiP2 in list_GPIaCC_P_01:
                                  if ((iccgpiP2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccgpiP2), _t='ORF-pozice: '+str(iccgpiP2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC a GPI strukturou - zadané vlákno - 3. čtecí rámec:')
                              for iccgpiP3 in list_GPIaCC_P_01:
                                  if ((iccgpiP3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccgpiP3), _t='ORF-pozice: '+str(iccgpiP3)+'\n')
                                  else:
                                      None # # #
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s TM, CC a GPI strukturou - zadané vlákno - 1. čtecí rámec:')
                              for itmccgpiP1 in list_TMaCCaGPI_P_01:
                                  if ((itmccgpiP1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccgpiP1), _t='ORF-pozice: '+str(itmccgpiP1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s TM, CC a GPI strukturou - zadané vlákno - 2. čtecí rámec:')
                              for itmccgpiP2 in list_TMaCCaGPI_P_01:
                                  if ((itmccgpiP2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccgpiP2), _t='ORF-pozice: '+str(itmccgpiP2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s TM, CC a GPI strukturou - zadané vlákno - 3. čtecí rámec:')
                              for itmccgpiP3 in list_TMaCCaGPI_P_01:
                                  if ((itmccgpiP3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccgpiP3), _t='ORF-pozice: '+str(itmccgpiP3)+'\n')
                                  else:
                                      None # # #
              a.hr()
              with a.h2(id="id9", klass='main_header'):a('ORFy - komplementární vlákno DNA:')
              with a.table():
                  with a.tr():
                      with a.td():a.img(src="figures/NA-ORFs-Nframe1-complementary_string-2D.png", alt="")
                      with a.td():a.img(src="figures/NA-ORFs-Nframe2-complementary_string-2D.png", alt="")
                      with a.td():a.img(src="figures/NA-ORFs-Nframe3-complementary_string-2D.png", alt="")
                  with a.tr():
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno.png", alt="")
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno.png", alt="")
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno.png", alt="")

              with a.table():
                  with a.tr():
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_TM_strukturou,AMK,komplementarni_vlakno-BIG.png", alt="")
                  with a.tr():
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_CC_strukturou,pravdepodobnost,komplementarni_vlakno-BIG.png", alt="")
                  with a.tr():
                      with a.td():a.img(src="figures/Souhrn_ORFu_s_GPI_strukturou,komplementarni_vlakno-BIG.png", alt="")
                      
              with a.table(cellpadding='5', border='2', bordercolor='black'): #
                   with a.tr():
                      with a.th(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s TM strukturou - komplementární vlákno')
                      with a.th(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s CC strukturou - komplementární vlákno')
                      with a.th(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s GPI strukturou - komplementární vláno')
                   with a.tr(): #
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('1. čtecí rámec:')
                              for itmN1all in listTMsN_pozice:
                                      if ((itmN1all[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmN1all), _t='ORF-pozice: '+str(itmN1all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('1. čtecí rámec:')
                              for iccN1all in CC_all_N_coords:
                                      if ((iccN1all[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccN1all), _t='ORF-pozice: '+str(iccN1all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('1. čtecí rámec:')
                              for igpiN1all in GPI_XYcoords_N_R:
                                      if ((igpiN1all[0])%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiN1all), _t='ORF-pozice: '+str(igpiN1all)+'\n')
                                      else:
                                          None
                   with a.tr(): ##
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('2. čtecí rámec:')
                              for itmN2all in listTMsN_pozice:
                                      if ((itmN2all[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmN2all), _t='ORF-pozice: '+str(itmN2all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('2. čtecí rámec:')
                              for iccN2all in CC_all_N_coords:
                                      if ((iccN2all[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccN2all), _t='ORF-pozice: '+str(iccN2all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('2. čtecí rámec:')
                              for igpiN2all in GPI_XYcoords_N_R:
                                      if ((igpiN2all[0]+2)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiN2all), _t='ORF-pozice: '+str(igpiN2all)+'\n')
                                      else:
                                          None
                   with a.tr(): ###
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('3. čtecí rámec:')
                              for itmN3all in listTMsN_pozice:
                                      if ((itmN3all[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(itmN3all), _t='ORF-pozice: '+str(itmN3all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('3. čtecí rámec:')
                              for iccN3all in CC_all_N_coords:
                                      if ((iccN3all[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(iccN3all), _t='ORF-pozice: '+str(iccN3all)+'\n')
                                      else:
                                          None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('3. čtecí rámec:')
                              for igpiN3all in GPI_XYcoords_N_R:
                                      if ((igpiN3all[0]+1)%3 == 0)==True:
                                          with a.p(align='center'):a.a(href='#'+str(igpiN3all), _t='ORF-pozice: '+str(igpiN3all)+'\n')
                                      else:
                                          None

              with a.table(cellpadding='5', border='2', bordercolor='black'): # width='99%'
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM strukturou - komplementární vlákno - 1. čtecí rámec:')
                              for itmN1 in list_jenTM_N_01:
                                  if ((itmN1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmN1), _t='ORF-pozice: '+str(itmN1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM strukturou - komplementární vlákno - 2. čtecí rámec:')
                              for itmN2 in list_jenTM_N_01:
                                  if ((itmN2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmN2), _t='ORF-pozice: '+str(itmN2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM strukturou - komplementární vlákno - 3. čtecí rámec:')
                              for itmN3 in list_jenTM_N_01:
                                  if ((itmN3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmN3), _t='ORF-pozice: '+str(itmN3)+'\n')
                                  else:
                                      None # # #
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC strukturou - komplementární vlákno - 1. čtecí rámec:')
                              for iccN1 in list_jenCC_N_01:
                                  if ((iccN1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccN1), _t='ORF-pozice: '+str(iccN1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC strukturou - komplementární vlákno - 2. čtecí rámec:')
                              for iccN2 in list_jenCC_N_01:
                                  if ((iccN2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccN2), _t='ORF-pozice: '+str(iccN2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC strukturou - komplementární vlákno - 3. čtecí rámec:')
                              for iccN3 in list_jenCC_N_01:
                                  if ((iccN3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccN3), _t='ORF-pozice: '+str(iccN3)+'\n')
                                  else:
                                      None # # # 
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s GPI strukturou - komplementární vlákno - 1. čtecí rámec:')
                              for igpiN1 in list_jenGPI_N_01:
                                  if ((igpiN1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(igpiN1), _t='ORF-pozice: '+str(igpiN1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s GPI strukturou - komplementární vlákno - 2. čtecí rámec:')
                              for igpiN2 in list_jenGPI_N_01:
                                  if ((igpiN2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(igpiN2), _t='ORF-pozice: '+str(igpiN2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s GPI strukturou - komplementární vlákno - 3. čtecí rámec:')
                              for igpiN3 in list_jenGPI_N_01:
                                  if ((igpiN3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(igpiN3), _t='ORF-pozice: '+str(igpiN3)+'\n')
                                  else:
                                      None
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a CC strukturou - komplementární vlákno - 1. čtecí rámec:')
                              for itmccN1 in list_TMaCC_N_01:
                                  if ((itmccN1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccN1), _t='ORF-pozice: '+str(itmccN1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a CC strukturou - komplementární vlákno - 2. čtecí rámec:')
                              for itmccN2 in list_TMaCC_N_01:
                                  if ((itmccN2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccN2), _t='ORF-pozice: '+str(itmccN2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a CC strukturou - komplementární vlákno - 3. čtecí rámec:')
                              for itmccN3 in list_TMaCC_N_01:
                                  if ((itmccN3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccN3), _t='ORF-pozice: '+str(itmccN3)+'\n')
                                  else:
                                      None # # #
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a GPI strukturou - komplementární vlákno - 1. čtecí rámec:')
                              for itmgpiN1 in list_TMaGPI_N_01:
                                  if ((itmgpiN1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmgpiN1), _t='ORF-pozice: '+str(itmgpiN1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a GPI strukturou - komplementární vlákno - 2. čtecí rámec:')
                              for itmgpiN2 in list_TMaGPI_N_01:
                                  if ((itmgpiN2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmgpiN2), _t='ORF-pozice: '+str(itmgpiN2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s TM a GPI strukturou - komplementární vlákno - 3. čtecí rámec:')
                              for itmgpiN3 in list_TMaGPI_N_01:
                                  if ((itmgpiN3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmgpiN3), _t='ORF-pozice: '+str(itmgpiN3)+'\n')
                                  else:
                                      None # # #
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC a GPI strukturou - komplementární vlákno - 1. čtecí rámec:')
                              for iccgpiN1 in list_GPIaCC_N_01:
                                  if ((iccgpiN1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccgpiN1), _t='ORF-pozice: '+str(iccgpiN1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC a GPI strukturou - komplementární vlákno - 2. čtecí rámec:')
                              for iccgpiN2 in list_GPIaCC_N_01:
                                  if ((iccgpiN2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccgpiN2), _t='ORF-pozice: '+str(iccgpiN2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy jen s CC a GPI strukturou - komplementární vlákno - 3. čtecí rámec:')
                              for iccgpiN3 in list_GPIaCC_N_01:
                                  if ((iccgpiN3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(iccgpiN3), _t='ORF-pozice: '+str(iccgpiN3)+'\n')
                                  else:
                                      None # # #
                   with a.tr():
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s TM, CC a GPI strukturou - komplementární vlákno - 1. čtecí rámec:')
                              for itmccgpiN1 in list_TMaCCaGPI_N_01:
                                  if ((itmccgpiN1[0])%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccgpiN1), _t='ORF-pozice: '+str(itmccgpiN1)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s TM, CC a GPI strukturou - komplementární vlákno - 2. čtecí rámec:')
                              for itmccgpiN2 in list_TMaCCaGPI_N_01:
                                  if ((itmccgpiN2[0]+2)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccgpiN2), _t='ORF-pozice: '+str(itmccgpiN2)+'\n')
                                  else:
                                      None
                      with a.td(width='33%'):
                          with a.p(align='center'):
                              with a.b():a('ORFy s TM, CC a GPI strukturou - komplementární vlákno - 3. čtecí rámec:')
                              for itmccgpiN3 in list_TMaCCaGPI_N_01:
                                  if ((itmccgpiN3[0]+1)%3 == 0)==True:
                                      with a.p(align='center'):a.a(href='#'+str(itmccgpiN3), _t='ORF-pozice: '+str(itmccgpiN3)+'\n')
                                  else:
                                      None # # #

              a.hr()
              with a.h2(id="id10", klass='main_header'):
                a("jednotlivé ORFy (zadané vlákno):")
                if (len(listORFsP) == 0):
                    with a.p():a('Žádný ORF a peptid v zadaném vláknu.')
                else:
                    for iximgP1 in listORFsP_pozice:
                          if ((iximgP1[0])%3 == 0)==True:
                            if (os.path.isfile('figures/given/NA string - fig No. ' + str(iximgP1) + ' provided_string.png'))==True:
                                  with a.p():a('ORF, zadané vlákno, 1. čtecí rámec - pozice: ')
                                  with a.b(id=str(iximgP1)):a(str(iximgP1) + ' z: ' + str(lenbP) + ':')
                                  a.img(src="figures/given/NA string - fig No. " + str(iximgP1) + " provided_string.png", width = "100%", alt="")
                                  with a.table(cellpadding='2', border='2', bordercolor='black'):
                                      with a.tr():
                                          with a.td():
                                              if (os.path.isfile('figures/given/Fig No. ' + str(iximgP1) + ' provided_string.png'))==True:
                                                  a.img(src="figures/given/Fig No. " + str(iximgP1) + " provided_string.png", alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td():
                                              if (os.path.isfile('CC_plots/plotsP-realCC/' + str(iximgP1[0]) + '.png'))==True:
                                                  a.img(src='CC_plots/plotsP-realCC/' + str(iximgP1[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td(align='center'): # bgcolor='green', bgcolor='white'
                                              if (os.path.isfile('figures/GPI/P' + str(iximgP1) + '.png'))==True:
                                                  a.img(src='figures/GPI/P' + str(iximgP1) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None

                                  with a.p(align='center'):
                                    a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                    #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>ORF_'+str(iximgP1[0])+'-'+str(iximgP1[1]))
                                          a.br()
                                          a(str(listORFsP[listORFsP_pozice.index(iximgP1)]))
                                      #with a.td(width='100%'):a([listORFsP[listORFsP_pozice.index(iximgP1)][iSekP1:iSekP1+90] for iSekP1 in range(0,len(listORFsP[listORFsP_pozice.index(iximgP1)]),90)])

                                  with a.p(align='center'):
                                    a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>trORF_'+str(iximgP1[0])+'-'+str(iximgP1[1]))
                                          a.br()
                                          a(str(PepsP_c00[listORFsP_pozice.index(iximgP1)]))

                                  a.hr()
                                
                            #else:
                                #with a.p():a('Neni ORF a peptid v 1. ctecim ramci v zadanem vlaknu.')
                    with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                    a.hr()

                if (len(listORFsP) == 0):
                    with a.p():a('Zadny ORF a peptid v zadanem vlaknu.')
                else:
                    for iximgP2 in listORFsP_pozice:
                          if ((iximgP2[0]+2)%3 == 0)==True:
                            if (os.path.isfile('figures/given/NA string - fig No. ' + str(iximgP2) + ' provided_string.png'))==True:
                                  with a.p():a('ORF, zadané vlákno, 2. čtecí rámec - pozice: ')
                                  with a.b(id=str(iximgP2)):a(str(iximgP2) + ' z: ' + str(lenbP) + ':')
                                  a.img(src="figures/given/NA string - fig No. " + str(iximgP2) + " provided_string.png", width = "100%", alt="")
                                  with a.table(cellpadding='2', border='2', bordercolor='black'):
                                      with a.tr():
                                          with a.td():
                                              if (os.path.isfile('figures/given/Fig No. ' + str(iximgP2) + ' provided_string.png'))==True:
                                                  a.img(src="figures/given/Fig No. " + str(iximgP2) + " provided_string.png", alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td():
                                              if (os.path.isfile('CC_plots/plotsP-realCC/' + str(iximgP2[0]) + '.png'))==True:
                                                  a.img(src='CC_plots/plotsP-realCC/' + str(iximgP2[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td(align='center'): # bgcolor='green'
                                              if (os.path.isfile('figures/GPI/P' + str(iximgP2) + '.png'))==True:
                                                  a.img(src='figures/GPI/P' + str(iximgP2) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                                  
                                  with a.p(align='center'):
                                    a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                    #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>ORF_'+str(iximgP2[0])+'-'+str(iximgP2[1]))
                                          a.br()
                                          a(str(listORFsP[listORFsP_pozice.index(iximgP2)]))
                                      #with a.td(width='100%'):a([listORFsP[listORFsP_pozice.index(iximgP2)][iSekP2:iSekP2+90] for iSekP2 in range(0,len(listORFsP[listORFsP_pozice.index(iximgP2)]),90)])

                                  with a.p(align='center'):
                                    a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>trORF_'+str(iximgP2[0])+'-'+str(iximgP2[1]))
                                          a.br()
                                          a(str(PepsP_c00[listORFsP_pozice.index(iximgP2)]))

                                  a.hr() 
                            #else:
                                #with a.p():a('Neni ORF a peptid ve 2. ctecim ramci v zadanem vlaknu.')
                    with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                    a.hr()

                if (len(listORFsP) == 0):
                    with a.p():a('Zadny ORF a peptid v zadanem vlaknu.')
                else:
                    for iximgP3 in listORFsP_pozice:
                          if ((iximgP3[0]+1)%3 == 0)==True:
                            if (os.path.isfile('figures/given/NA string - fig No. ' + str(iximgP3) + ' provided_string.png'))==True:
                                  with a.p():a('ORF, zadané vlákno, 3. čtecí rámec - pozice: ')
                                  with a.b(id=str(iximgP3)):a(str(iximgP3) + ' z: ' + str(lenbP) + ':')
                                  a.img(src="figures/given/NA string - fig No. " + str(iximgP3) + " provided_string.png", width = "100%", alt="")
                                  with a.table(align="center", cellpadding='2', border='2', bordercolor='black'):
                                      with a.tr():
                                          with a.td():
                                              if (os.path.isfile('figures/given/Fig No. ' + str(iximgP3) + ' provided_string.png'))==True:
                                                  a.img(src="figures/given/Fig No. " + str(iximgP3) + " provided_string.png", alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td():
                                              if (os.path.isfile('CC_plots/plotsP-realCC/' + str(iximgP3[0]) + '.png'))==True:
                                                  a.img(src='CC_plots/plotsP-realCC/' + str(iximgP3[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td(align='center'): # bgcolor='green'
                                              if (os.path.isfile('figures/GPI/P' + str(iximgP3) + '.png'))==True:
                                                  a.img(src='figures/GPI/P' + str(iximgP3) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None

                                  with a.p(align='center'):
                                    a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                    #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>ORF_'+str(iximgP3[0])+'-'+str(iximgP3[1]))
                                          a.br()
                                          a(str(listORFsP[listORFsP_pozice.index(iximgP3)]))
                                      #with a.td(width='100%'):a([listORFsP[listORFsP_pozice.index(iximgP3)][iSekP3:iSekP3+90] for iSekP3 in range(0,len(listORFsP[listORFsP_pozice.index(iximgP3)]),90)])


                                  with a.p(align='center'):
                                    a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>trORF_'+str(iximgP3[0])+'-'+str(iximgP3[1]))
                                          a.br()
                                          a(str(PepsP_c00[listORFsP_pozice.index(iximgP3)]))

                                                  
                                  a.hr()
                            #else:
                                #with a.p():a('Neni ORF a peptid ve 3. ctecim ramci v zadanem vlaknu.')
                    with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                    a.hr()

              a.hr()
              a.hr()
              with a.h2(id="id11", klass='main_header'):
                a("jednotlivé ORFy (komplementární vlákno):")
                if (len(listORFsN) == 0):
                    with a.p():a('Žádný ORF a peptid v komplementárním vláknu.')
                else:
                    for iximgN1 in listORFsN_pozice:
                          if ((iximgN1[0])%3 == 0)==True:
                            if (os.path.isfile('figures/complementary/NA string - fig No. ' + str(iximgN1) + ' complementary_string.png'))==True:
                                  with a.p():a('ORF, komplementární vlákno, 1. čtecí rámec - pozice: ')
                                  with a.b(id=str(iximgN1)):a(str(iximgN1) + ' z: ' + str(lenbN) + ':')
                                  a.img(src="figures/complementary/NA string - fig No. " + str(iximgN1) + " complementary_string.png", width = "100%", alt="")
                                  with a.table(cellpadding='2', border='2', bordercolor='black'):
                                      with a.tr():
                                          with a.td():
                                              if (os.path.isfile('figures/complementary/Fig No. ' + str(iximgN1) + ' complementary_string.png'))==True:
                                                  a.img(src="figures/complementary/Fig No. " + str(iximgN1) + " complementary_string.png", alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td():
                                              if (os.path.isfile('CC_plots/plotsN-realCC/' + str(iximgN1[0]) + '.png'))==True:
                                                  a.img(src='CC_plots/plotsN-realCC/' + str(iximgN1[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td(align='center'): # bgcolor='green' , bgcolor='white'
                                              if (os.path.isfile('figures/GPI/N/' + str(iximgN1) + '.png'))==True:
                                                  a.img(src='figures/GPI/N/' + str(iximgN1) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None

                                  with a.p(align='center'):
                                    a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                    #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>ORF_'+str(iximgN1[0])+'-'+str(iximgN1[1]))
                                          a.br()
                                          a(str(listORFsN[listORFsN_pozice.index(iximgN1)]))
                                      #with a.td(width='100%'):a([listORFsN[listORFsN_pozice.index(iximgN1)][iSekN1:iSekN1+90] for iSekN1 in range(0,len(listORFsN[listORFsN_pozice.index(iximgN1)]),90)])

                                  with a.p(align='center'):
                                    a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>trORF_'+str(iximgN1[0])+'-'+str(iximgN1[1]))
                                          a.br()
                                          a(str(PepsN_c00[listORFsN_pozice.index(iximgN1)]))
                                                  
                                  a.hr()
                            #else:
                                #with a.p():a('Neni ORF a peptid v 1. ctecim ramci v komplementarnim vlaknu.')
                    with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                    a.hr()

                if (len(listORFsN) == 0):
                    with a.p():a('Žádný ORF a peptid v komplementárním vláknu.')
                else:
                    for iximgN2 in listORFsN_pozice:
                          if ((iximgN2[0]+2)%3 == 0)==True:
                            if (os.path.isfile('figures/complementary/NA string - fig No. ' + str(iximgN2) + ' complementary_string.png'))==True:
                                  with a.p():a('ORF, komplementární vlákno, 2. čtecí rámec - pozice: ')
                                  with a.b(id=str(iximgN2)):a(str(iximgN2) + ' z: ' + str(lenbN) + ':')
                                  a.img(src="figures/complementary/NA string - fig No. " + str(iximgN2) + " complementary_string.png", width = "100%", alt="")
                                  with a.table(cellpadding='2', border='2', bordercolor='black'):
                                      with a.tr():
                                          with a.td():
                                              if (os.path.isfile('figures/complementary/Fig No. ' + str(iximgN2) + ' complementary_string.png'))==True:
                                                  a.img(src="figures/complementary/Fig No. " + str(iximgN2) + " complementary_string.png", alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td():
                                              if (os.path.isfile('CC_plots/plotsN-realCC/' + str(iximgN2[0]) + '.png'))==True:
                                                  a.img(src='CC_plots/plotsN-realCC/' + str(iximgN2[0]) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td(align='center'): # bgcolor='green' , bgcolor='white'
                                              if (os.path.isfile('figures/GPI/N/' + str(iximgN2) + '.png'))==True:
                                                  a.img(src='figures/GPI/N/' + str(iximgN2) + '.png', alt="", width="550", height="450", align="absmiddle")
                                              else:
                                                  None

                                  with a.p(align='center'):
                                    a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                    #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>ORF_'+str(iximgN2[0])+'-'+str(iximgN2[1]))
                                          a.br()
                                          a(str(listORFsN[listORFsN_pozice.index(iximgN2)]))
                                      #with a.td(width='100%'):a([listORFsN[listORFsN_pozice.index(iximgN2)][iSekN2:iSekN2+90] for iSekN2 in range(0,len(listORFsN[listORFsN_pozice.index(iximgN2)]),90)])

                                  with a.p(align='center'):
                                    a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>trORF_'+str(iximgN2[0])+'-'+str(iximgN2[1]))
                                          a.br()
                                          a(str(PepsN_c00[listORFsN_pozice.index(iximgN2)]))

                                  a.hr()
                            #else:
                                #with a.p():a('Neni ORF a peptid ve 2. ctecim ramci v komplementarnim vlaknu.')
                    with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                    a.hr()

                if (len(listORFsN) == 0):
                    with a.p():a('Žádný ORF a peptid v komplementárním vláknu.')
                else:
                    for iximgN3 in listORFsN_pozice:
                          if ((iximgN3[0]+1)%3 == 0)==True:
                            if (os.path.isfile('figures/complementary/NA string - fig No. ' + str(iximgN3) + ' complementary_string.png'))==True:
                                  with a.p():a('ORF, komplementární vlákno, 3. čtecí rámec - pozice: ')
                                  with a.b(id=str(iximgN3)):a(str(iximgN3) + ' z: ' + str(lenbN) + ':')
                                  a.img(src="figures/complementary/NA string - fig No. " + str(iximgN3) + " complementary_string.png", width = "100%", alt="")
                                  with a.table(cellpadding='2', border='2', bordercolor='black'):
                                      with a.tr(width = "100%"):
                                          with a.td():
                                              if (os.path.isfile('figures/complementary/Fig No. ' + str(iximgN3) + ' complementary_string.png'))==True:
                                                  a.img(src="figures/complementary/Fig No. " + str(iximgN3) + " complementary_string.png", alt="", width = "550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td():
                                              if (os.path.isfile('CC_plots/plotsN-realCC/' + str(iximgN3[0]) + '.png'))==True:
                                                  a.img(src='CC_plots/plotsN-realCC/' + str(iximgN3[0]) + '.png', alt="", width = "550", height="450", align="absmiddle")
                                              else:
                                                  None
                                          with a.td(align='center'): # bgcolor='white'
                                              if (os.path.isfile('figures/GPI/N/' + str(iximgN3) + '.png'))==True:
                                                  a.img(src='figures/GPI/N/' + str(iximgN3) + '.png', alt="", width = "550", height="450", align="absmiddle")
                                              else:
                                                  None

                                  with a.p(align='center'):
                                    a('Sekvence nukleotidů v rámci ORFu (5` --> 3`):')
                                    #a('Sekvence nukleotidů (děleno po 90ti nt) v rámci ORFu (5` --> 3`):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>ORF_'+str(iximgN3[0])+'-'+str(iximgN3[1]))
                                          a.br()
                                          a(str(listORFsN[listORFsN_pozice.index(iximgN3)]))
                                          #with a.pre():a(str(listORFsN[listORFsN_pozice.index(iximgN3)]))
                                      #with a.td(width='100%'):a([listORFsN[listORFsN_pozice.index(iximgN3)][iSekN3:iSekN3+90] for iSekN3 in range(0,len(listORFsN[listORFsN_pozice.index(iximgN3)]),90)])

                                  with a.p(align='center'):
                                    a('Sekvence aminokyselin v rámci ORFu (N --> C):')
                                  with a.table(cellpadding='5', border='2', bordercolor='black', width='100%'):
                                    with a.tr(width='100%'):
                                      with a.td(width='100%'):
                                          a('>trORF_'+str(iximgN3[0])+'-'+str(iximgN3[1]))
                                          a.br()
                                          a(str(PepsN_c00[listORFsN_pozice.index(iximgN3)]))
                                          #with a.pre():a(str(PepsN_c00[listORFsN_pozice.index(iximgN3)]))
                
                                  a.hr()
                            #else:
                                #with a.p():a('Neni ORF a peptid ve 3. ctecim ramci v komplementarnim vlaknu.')
                    with a.p():a.a(href='#top', _t='Zpět na začátek stránky.')
                    a.hr()
                    
              #with a.p(align="center"):
              #a('(konec stránky)')

              a.hr()
              a.hr()

        index23_str_a = str(a)

        index23.write(index23_str_a)
        index23.close()

        print('\n')
        print('Vysledek byl ulozen do HTML souboru: "index23.html".\n')

elapsed = time.time() - t
print('Elapsed time (seconds): ',elapsed)
