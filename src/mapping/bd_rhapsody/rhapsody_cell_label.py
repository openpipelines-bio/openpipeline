#!/usr/bin/env python

# copied from https://bd-rhapsody-public.s3.amazonaws.com/CellLabel/rhapsody_cell_label.py.txt
# documented at https://bd-rhapsody-bioinfo-docs.genomics.bd.com/steps/steps_cell_label.html

"""
Rhapsody cell label structure
Information on the cell label is captured by the combination of bases in three cell label sections (CLS1, CLS2, CLS3).
Two common linker sequences (L1, L2) separate the three CLS.

--CLS1---|-L1-|--CLS2---|-L2-|--CL3---|--UMI---|-CaptureSequence-


Each cell label section has a whitelist of 96 or 384 possible 9 base sequences.
All the capture oligos from a single bead will have the same cell label.

----------------

V1 beads:

[A96_cell_key1] + [v1_linker1] + [A96_cell_key2] + [v1_linker2] + [A96_cell_key3] + [8 random base UMI] + [18 base polyT capture]


----------------

Enhanced beads:
Enhanced beads contain two different capture oligo types, polyT and 5prime.  On any one bead, the two different capture oligo types have the same cell label sequences.
Compared to the V1 bead, enhanced beads have shorter linker sequences, longer polyT, and 0-3 diversity insert bases at the beginning of the sequence.
The cell label sections use the same 3 sequence whitelists as V1 beads.

polyT capture oligo:
[Enh_insert 0-3 bases] + [A96_cell_key1] + [Enh_linker1] + [A96_cell_key2] + [Enh_linker2] + [A96_cell_key3] + [8 random base UMI] + [25 base polyT capture]

5prime capture oligo:
[Enh_5p_primer] + [A96_cell_key1] + [Enh_5p_linker1] + [A96_cell_key2] + [Enh_5p_linker2] + [A96_cell_key3] + [8 random base UMI] + [Tso_capture_seq]


----------------

Enhanced V2/V3 beads:
Enhanced V2/V3 beads have the same structure as Enhanced beads, but the cell label sections have been updated with increased diversity


polyT capture oligo:
[Enh_insert 0-3 bases] + [B384_cell_key1] + [Enh_linker1] + [B384_cell_key2] + [Enh_linker2] + [B384_cell_key3] + [8 random base UMI] + [25 base polyT capture]

5prime capture oligo:
[Enh_5p_primer] + [B384_cell_key1] + [Enh_5p_linker1] + [B384_cell_key2] + [Enh_5p_linker2] + [B384_cell_key3] + [8 random base UMI] + [Tso_capture_seq]


The only difference between Enh V2 and Enh V3 beads is a different Tso_capture_seq.

----------------

The Rhapsody Sequence Analysis Pipeline will convert each cell label into a single integer representing a unique cell label sequence - which is used in the output files as the 'Cell_index'.
This cell index integer is deterministic and derived from the 3 part cell label as follows: 

- Get the 1-based index for each cell label section from the python sets of sequences below
- Apply this equation:
    (CLS1index - 1) * 384 * 384 + (CLS2index - 1) * 384 + CLS3index

(See label_sections_to_index() function below)


Example: Enhanced bead sequence:
ACACATTGCAGTGAAGATAGTTCGACACTCAAGACA

Each part identified:
A                CACATTGCA          GTGA      AGATAGTTC          GACA      CTCAAGACA
DiversityInsert  A96_cell_key1-33   Linker1   A96_cell_key2-78   Linker2   A96_cell_key3-21

33-78-21
(33 - 1) * 384 * 384 + (78 - 1) * 384 + 21
=4748181


The original sequences of cell label can be determined from the cell index integer by reversing this conversion.
See index_to_label_sections() and index_to_sequence() functions below.

"""

v1_linker1 = 'ACTGGCCTGCGA'
v1_linker2 = 'GGTAGCGGTGACA'

Enh_linker1 = 'GTGA'
Enh_linker2 = 'GACA'

Enh_5p_primer = "ACAGGAAACTCATGGTGCGT"

Enh_5p_linker1 = "AATG"
Enh_5p_linker2 = "CCAC"

Enh_inserts = ["", "A", "GT", "TCA"]

Tso_capture_seq_Enh_EnhV2 = "TATGCGTAGTAGGTATG"
Tso_capture_seq_EnhV3 = "GTGGAGTCGTGATTATA"

A96_cell_key1 = ("GTCGCTATA","CTTGTACTA","CTTCACATA","ACACGCCGG","CGGTCCAGG","AATCGAATG","CCTAGTATA","ATTGGCTAA","AAGACATGC","AAGGCGATC",
                 "GTGTCCTTA","GGATTAGGA","ATGGATCCA","ACATAAGCG","AACTGTATT","ACCTTGCGG","CAGGTGTAG","AGGAGATTA","GCGATTACA","ACCGGATAG",
                 "CCACTTGGA","AGAGAAGTT","TAAGTTCGA","ACGGATATT","TGGCTCAGA","GAATCTGTA","ACCAAGGAC","AGTATCTGT","CACACACTA","ATTAAGTGC",
                 "AAGTAACCC","AAATCCTGT","CACATTGCA","GCACTGTCA","ATACTTAGG","GCAATCCGA","ACGCAATCA","GAGTATTAG","GACGGATTA","CAGCTGACA",
                 "CAACATATT","AACTTCTCC","CTATGAAAT","ATTATTACC","TACCGAGCA","TCTCTTCAA","TAAGCGTTA","GCCTTACAA","AGCACACAG","ACAGTTCCG",
                 "AGTAAAGCC","CAGTTTCAC","CGTTACTAA","TTGTTCCAA","AGAAGCACT","CAGCAAGAT","CAAACCGCC","CTAACTCGC","AATATTGGG","AGAACTTCC",
                 "CAAAGGCAC","AAGCTCAAC","TCCAGTCGA","AGCCATCAC","AACGAGAAG","CTACAGAAC","AGAGCTATG","GAGGATGGA","TGTACCTTA","ACACACAAA",
                 "TCAGGAGGA","GAGGTGCTA","ACCCTGACC","ACAAGGATC","ATCCCGGAG","TATGTGGCA","GCTGCCAAT","ATCAGAGCT","TCGAAGTGA","ATAGACGAG",
                 "AGCCCAATC","CAGAATCGT","ATCTCCACA","ACGAAAGGT","TAGCTTGTA","ACACGAGAT","AACCGCCTC","ATTTAGATG","CAAGCAAGC","CAAAGTGTG",
                 "GGCAAGCAA","GAGCCAATA","ATGTAATGG","CCTGAGCAA","GAGTACATT","TGCGATCTA"
                 )

A96_cell_key2 = ("TACAGGATA","CACCAGGTA","TGTGAAGAA","GATTCATCA","CACCCAAAG","CACAAAGGC","GTGTGTCGA","CTAGGTCCT","ACAGTGGTA","TCGTTAGCA",
                 "AGCGACACC","AAGCTACTT","TGTTCTCCA","ACGCGAAGC","CAGAAATCG","ACCAAAATG","AGTGTTGTC","TAGGGATAC","AGGGCTGGT","TCATCCTAA",
                 "AATCCTGAA","ATCCTAGGA","ACGACCACC","TTCCATTGA","TAGTCTTGA","ACTGTTAGA","ATTCATCGT","ACTTCGAGC","TTGCGTACA","CAGTGCCCG",
                 "GACACTTAA","AGGAGGCGC","GCCTGTTCA","GTACATCTA","AATCAGTTT","ACGATGAAT","TGACAGACA","ATTAGGCAT","GGAGTCTAA","TAGAACACA",
                 "AAATAAATA","CCGACAAGA","CACCTACCC","AAGAGTAGA","TCATTGAGA","GACCTTAGA","CAAGACCTA","GGAATGATA","AAACGTACC","ACTATCCTC",
                 "CCGTATCTA","ACACATGTC","TTGGTATGA","GTGCAGTAA","AGGATTCAA","AGAATGGAG","CTCTCTCAA","GCTAACTCA","ATCAACCGA","ATGAGTTAC",
                 "ACTTGATGA","ACTTTAACT","TTGGAGGTA","GCCAATGTA","ATCCAACCG","GATGAACTG","CCATGCACA","TAGTGACTA","AAACTGCGC","ATTACCAAG",
                 "CACTCGAGA","AACTCATTG","CTTGCTTCA","ACCTGAGTC","AGGTTCGCT","AAGGACTAT","CGTTCGGTA","AGATAGTTC","CAATTGATC","GCATGGCTA",
                 "ACCAGGTGT","AGCTGCCGT","TATAGCCCT","AGAGGACCA","ACAATATGG","CAGCACTTC","CACTTATGT","AGTGAAAGG","AACCCTCGG","AGGCAGCTA",
                 "AACCAAAGT","GAGTGCGAA","CGCTAAGCA","AATTATAAC","TACTAGTCA","CAACAACGG"
                 )

A96_cell_key3 = ("AAGCCTTCT","ATCATTCTG","CACAAGTAT","ACACCTTAG","GAACGACAA","AGTCTGTAC","AAATTACAG","GGCTACAGA","AATGTATCG","CAAGTAGAA",
                 "GATCTCTTA","AACAACGCG","GGTGAGTTA","CAGGGAGGG","TCCGTCTTA","TGCATAGTA","ACTTACGAT","TGTATGCGA","GCTCCTTGA","GGCACAACA",
                 "CTCAAGACA","ACGCTGTTG","ATATTGTAA","AAGTTTACG","CAGCCTGGC","CTATTAGCC","CAAACGTGG","AAAGTCATT","GTCTTGGCA","GATCAGCGA",
                 "ACATTCGGC","AGTAATTAG","TGAAGCCAA","TCTACGACA","CATAACGTT","ATGGGACTC","GATAGAGGA","CTACATGCG","CAACGATCT","GTTAGCCTA",
                 "AGTTGCATC","AAGGGAACT","ACTACATAT","CTAAGCTTC","ACGAACCAG","TACTTCGGA","AACATCCAT","AGCCTGGTT","CAAGTTTCC","CAGGCATTT",
                 "ACGTGGGAG","TCTCACGGA","GCAACATTA","ATGGTCCGT","CTATCATGA","CAATACAAG","AAAGAGGCC","GTAGAAGCA","GCTATGGAA","ACTCCAGGG",
                 "ACAAGTGCA","GATGGTCCA","TCCTCAATA","AATAAACAA","CTGTACGGA","CTAGATAGA","AGCTATGTG","AAATGGAGG","AGCCGCAAG","ACAGTAAAC",
                 "AACGTGTGA","ACTGAATTC","AAGGGTCAG","TGTCTATCA","TCAGATTCA","CACGATCCG","AACAGAAAC","CATGAATGA","CGTACTACG","TTCAGCTCA",
                 "AAGGCCGCA","GGTTGGACA","CGTCTAGGT","AATTCGGCG","CAACCTCCA","CAATAGGGT","ACAGGCTCC","ACAACTAGT","AGTTGTTCT","AATTACCGG",
                 "ACAAACTTT","TCTCGGTTA","ACTAGACCG","ACTCATACG","ATCGAGTCT","CATAGGTCA"
                 )

B384_cell_key1 = ("TGTGTTCGC","TGTGGCGCC","TGTCTAGCG","TGGTTGTCC","TGGTTCCTC","TGGTGTGCT","TGGCGACCG","TGCTGTGGC","TGCTGGCAC","TGCTCTTCC",
                  "TGCCTCACC","TGCCATTAT","TGATGTCTC","TGATGGCCT","TGATGCTTG","TGAAGGACC","TCTGTCTCC","TCTGATTAT","TCTGAGGTT","TCTCGTTCT",
                  "TCTCATCCG","TCCTGGATT","TCAGCATTC","TCACGCCTT","TATGTGCAC","TATGCGGCC","TATGACGAG","TATCTCGTG","TATATGACC","TAGGCTGTG",
                  "TACTGCGTT","TACGTGTCC","TAATCACAT","GTTGTGTTG","GTTGTGGCT","GTTGTCTGT","GTTGTCGAG","GTTGTCCTC","GTTGTATCC","GTTGGTTCT",
                  "GTTGGCGTT","GTTGGAGCG","GTTGCTGCC","GTTGCGCAT","GTTGCAGGT","GTTGCACTG","GTTGATGAT","GTTGATACG","GTTGAAGTC","GTTCTGTGC",
                  "GTTCTCTCG","GTTCTATAT","GTTCGTATG","GTTCGGCCT","GTTCGCGGC","GTTCGATTC","GTTCCGGTT","GTTCCGACG","GTTCACGCT","GTTATCACC",
                  "GTTAGTCCG","GTTAGGTGT","GTTAGAGAC","GTTAGACTT","GTTACCTCT","GTTAATTCC","GTTAAGCGC","GTGTTGCTT","GTGTTCGGT","GTGTTCCAG",
                  "GTGTTCATC","GTGTCACAC","GTGTCAAGT","GTGTACTGC","GTGGTTAGT","GTGGTACCG","GTGGCGATC","GTGCTTCTG","GTGCGTTCC","GTGCGGTAT",
                  "GTGCGCCTT","GTGCGAACT","GTGCAGCCG","GTGCAATTG","GTGCAAGGC","GTCTTGCGC","GTCTGGCCG","GTCTGAGGC","GTCTCAGAT","GTCTCAACC",
                  "GTCTATCGT","GTCGGTGTG","GTCGGAATC","GTCGCTCCG","GTCCTCGCC","GTCCTACCT","GTCCGCTTG","GTCCATTCT","GTCCAATAC","GTCATGTAT",

                  "GTCAGTGGT","GTCAGATAG","GTATTAACT","GTATCAGTC","GTATAGCCT","GTATACTTG","GTATAAGGT","GTAGCATCG","GTACCGTCC","GTACACCTC",
                  "GTAAGTGCC","GTAACAGAG","GGTTGTGTC","GGTTGGCTG","GGTTGACGC","GGTTCGTCG","GGTTCAGTT","GGTTATATT","GGTTAATAC","GGTGTACGT",
                  "GGTGCCGCT","GGTGCATGC","GGTCGTTGC","GGTCGAGGT","GGTAGGCAC","GGTAGCTTG","GGTACATAG","GGTAATCTG","GGCTTGGCC","GGCTTCACG",
                  "GGCTTATGT","GGCTTACTC","GGCTGTCTT","GGCTCTGTG","GGCTCCGGT","GGCTCACCT","GGCGTTGAG","GGCGTGTAC","GGCGTGCTG","GGCGTATCG",
                  "GGCGCTCGT","GGCGCTACC","GGCGAGCCT","GGCGAGATC","GGCGACTTG","GGCCTCTTC","GGCCTACAG","GGCCAGCGC","GGCCAACTT","GGCATTCCT",
                  "GGCATCCGC","GGCATAACC","GGCAACGAT","GGATGTCCG","GGATGAGAG","GGATCTGGC","GGATCCATG","GGATAGGTT","GGAGTCGTG","GGAGAAGGC",
                  "GGACTCCTT","GGACTAGTC","GGACCGTTG","GGAATTAGT","GGAATCTCT","GGAATCGAC","GGAAGCCTC","GCTTGTAGC","GCTTGACCG","GCTTCGGAC",
                  "GCTTCACAT","GCTTAGTCT","GCTGGATAT","GCTGGAACC","GCTGCGATG","GCTGATCAG","GCTGAGCGT","GCTCTTGTC","GCTCTCCTG","GCTCGGTCC",
                  "GCTCCAATT","GCTATTCGC","GCTATGAGT","GCTAGTGTT","GCTAGGATC","GCTAGCACT","GCTACGTAT","GCTAACCTT","GCGTTCCGC","GCGTGTGCC",
                  "GCGTGCATT","GCGTCGGTT","GCGTATGTG","GCGTATACT","GCGGTTCAC","GCGGTCTTG","GCGGCGTCG","GCGGCACCT","GCGCTGGAC","GCGCTCTCC",

                  "GCGCGGCAG","GCGCGATAC","GCGCCGACC","GCGAGCGAG","GCGAGAGGT","GCGAATTAC","GCCTTGCAT","GCCTGCGCT","GCCTAACTG","GCCGTCCGT",
                  "GCCGCTGTC","GCCATGCCG","GCCAGCTAT","GCCAACCAG","GCATGGTTG","GCATCGACG","GCAGGCTAG","GCAGGACGC","GCAGCCATC","GCAGATACC",
                  "GCAGACGTT","GCACTATGT","GCACACGAG","GATTGTCAT","GATTGGTAG","GATTGCACC","GATTCTACT","GATTCGCTT","GATTAGGCC","GATTACGGT",
                  "GATGTTGGC","GATGTTATG","GATGGCCAG","GATCGTTCG","GATCGGAGC","GATCGCCTC","GATCCTCTG","GATCCAGCG","GATACACGC","GAGTTACCT",
                  "GAGTCGTAT","GAGTCGCCG","GAGGTGTAG","GAGGCATTG","GAGCGGACG","GAGCCTGAG","GAGATCTGT","GAGATAATT","GAGACGGCT","GACTTCGTG",
                  "GACTGTTCT","GACTCTTAG","GACCGCATT","GAATTGAGC","GAATATTGC","GAAGGCTCT","GAAGAGACT","GAACTGCCG","GAACGCGTG","CTTGTGTAT",
                  "CTTGTGCGC","CTTGTCATG","CTTGGTCTT","CTTGGTACC","CTTGGATGT","CTTGCTCAC","CTTGCAATC","CTTGAGGCC","CTTGACGGT","CTTCTGATC",
                  "CTTCTCGTT","CTTCTAGGC","CTTCGTTAG","CTTATGTCC","CTTATGCTT","CTTATATAG","CTTAGGTTG","CTTAGGAGC","CTTACTTAT","CTGTTCTCG",
                  "CTGTGCCTC","CTGTCGCAT","CTGTCGAGC","CTGTAGCTG","CTGTACGTT","CTGCTTGCC","CTGCGTAGT","CTGCACACC","CTGATGGAT","CTGAGTCAT",
                  "CTGACGCCG","CTGAACGAG","CTCTTGTAG","CTCTTAGTT","CTCTTACCG","CTCTGCACC","CTCTCGTCC","CTCGTATTG","CTCGACTAT","CTCCTGACG",

                  "CTCACTAGC","CTATACGGC","CGTTCGCTC","CGTTCACCG","CGTATAGTT","CGGTGTTCC","CGGTGTCAG","CGGTCCTGC","CGGCGACTC","CGGCACGGT",
                  "CGGATAGCC","CGGAGAGAT","CGCTAATAG","CGCGTTGGC","CGCGCAGAG","CGCACTGCC","CCTTGTCTC","CCTTGGCGT","CCTTCTGAG","CCTTCTCCT",
                  "CCTTCGACC","CCTTACTTG","CCTGTTCGT","CCTGTATGC","CCTCGGCCG","CCGTTAATT","CCATGTGCG","CCAGTGGTT","CCAGGCATT","CCAGGATCC",
                  "CCAGCGTTG","CATTCCGAT","CATTATACC","CATGTTGAG","ATTGCGTGT","ATTGCGGAC","ATTGCGCCG","ATTGACTTG","ATTCGGCTG","ATTCGCGAG",
                  "ATTCCAAGT","ATTATCTTC","ATTACTGTT","ATTACACTC","ATGTTCTAT","ATGTTACGC","ATGTGTATC","ATGTGGCAG","ATGTCTGTG","ATGGTGCAT",
                  "ATGCTTACT","ATGCTGTCC","ATGCTCGGC","ATGAGGTTC","ATGAGAGTG","ATCTTGGCT","ATCTGTGCG","ATCGGTTCC","ATCATGCTC","ATCATCACT",
                  "ATATCTTAT","ATAGGCGCC","AGTTGGTAT","AGTTGAGCC","AGTGCGACC","AGGTGCTAC","AGGCTTGCG","AGGCCTTCC","AGGCACCTT","AGGAATATG",
                  "AGCGGCCAG","AGCCTGGTC","AGCCTGACT","AGCAATCCG","AGAGATGTT","AGAGAATTC","ACTCGCTTG","ACTCGACCT","ACGTACACC","ACGGATGGT",
                  "ACCAGTCTG","ACATTCGGC","ACATGAGGT","ACACTAATT"
                  )

B384_cell_key2 = ("TTGTGTTGT","TTGTGGTAG","TTGTGCGGA","TTGTCTGTT","TTGTCTAAG","TTGTCATAT","TTGTCACGA","TTGTATGAA","TTGTACAGT","TTGGTTAAT",
                  "TTGGTGCAA","TTGGTCGAG","TTGGTATTA","TTGGCACAG","TTGGATACA","TTGGAAGTG","TTGCGGTTA","TTGCCATTG","TTGCACGCG","TTGCAAGGT",
                  "TTGATGTAT","TTGATAATT","TTGAGACGT","TTGACTACT","TTGACCGAA","TTCTGGTCT","TTCTGCACA","TTCTCCTTA","TTCTCCGCT","TTCTAGGTA",
                  "TTCTAATCG","TTCGTCGTA","TTCGTAGAT","TTCGGCTTG","TTCGGAATA","TTCGCCAGA","TTCGATTGT","TTCGATCAG","TTCCTCGGT","TTCCGGCAG",
                  "TTCCGCATT","TTCCAATTA","TTCATTGAA","TTCATGCTG","TTCAGGAGT","TTCACTATA","TTCAACTCT","TTCAACGTG","TTATGCGTT","TTATGATTG",
                  "TTATCCTGT","TTATCCGAG","TTATATTAT","TTAGGCGCG","TTACTGGAA","TTACTAGTT","TTACGTGGT","TTACGATAT","TTACCTAGA","TTACATGAG",
                  "TTACAGCGT","TTACACGGA","TTACACACT","TTAATCAGT","TTAATAGGA","TTAAGTGTG","TTAACCTTG","TTAACACAA","TGTTCACTT","TGTTCAAGA",
                  "TGTTAAGTG","TGTGTTATG","TGTGTCCAA","TGTGGAGCG","TGTCAGTTA","TGTCAGAAG","TGGTTAGTT","TGGTTACAA","TGGCGTTAT","TGGCGCCAA",
                  "TGGAGTCTT","TGCGTATTG","TGATAGAGA","TGAGGTATT","TGAGAATCT","TCTTGGTAA","TCTTCATAG","TCTGTCCTT","TCTGGAATT","TCTACCGCG",
                  "TCGTTCGAA","TCGTCAGTG","TCGACGAGA","TCATGGCTT","TCACACTTA","TATTCCGAA","TATTATGGT","TATGCTATT","TATCAAGGA","TAGTTCAAT",

                  "TAGCTGCTT","TAGAGGAAG","TACCTGTTA","TACACCTGT","GTTGTGCGT","GTTGGCTAT","GTTGCCAAG","GTTGACCTT","GTTCTGCTA","GTTCTGAAT",
                  "GTTCTATCA","GTTCGCGTG","GTTCCTTAT","GTTAGCAGT","GTTACTGTG","GTTACTCAA","GTTAAGAGA","GTTAACTTA","GTGTCGGCA","GTGTCCATT",
                  "GTGCTTGAG","GTGCTCGTT","GTGCTCACA","GTGCCTGGA","GTCTTGTCG","GTCTTGATT","GTCTTCCGT","GTCTTAAGA","GTCTCATCT","GTCTACGAG",
                  "GTCGTTGCT","GTCGTGTTA","GTCGGTAAT","GTCGGATGT","GTCGAGCTG","GTCCGGACT","GTCCAACAT","GTCAGACGA","GTCAGAATT","GTCACTCTT",
                  "GTCAAGGAA","GTATGTCTT","GTATGTACA","GTATCGGTT","GTATATGTA","GTATACAAT","GTAGTTAAG","GTAGTCGAT","GTAGCCTTA","GTAGATACT",
                  "GTACGATTA","GTACAGTCT","GTAATTCGT","GCTTGGCAG","GCTTGCTTG","GCTTGAGGA","GCTTCATTA","GCTTATGCG","GCTGTGTAG","GCTGTCATG",
                  "GCTGGTTGT","GCTGGACTG","GCTGCCTAA","GCTGATATT","GCTCTTAGT","GCTCTATTG","GCTCGCCGT","GCTCCGCTG","GCTATTCTG","GCTATACGA",
                  "GCTACTAAG","GCTACATGT","GCTAACTCT","GCGTTGTAA","GCGTTCTCT","GCGTGCGTA","GCGTCTTGA","GCGTCCGAT","GCGTAAGAG","GCGCTTACG",
                  "GCGCGGATT","GCGCCATAT","GCGCATGAA","GCGATCAAT","GCGAGCCTT","GCGAGATTG","GCGAGAACA","GCCTTGGTA","GCCTTCTAG","GCCTTCACA",
                  "GCCTGAGTG","GCCTCACGT","GCCGGCGAA","GCCGCACAA","GCCATGCTT","GCCATATAT","GCCAATTCG","GCATTCGTT","GCATGATGT","GCAGTTGGA",

                  "GCAGTGTCT","GCACTTGTG","GCAATCTGT","GCAACACTT","GATTGTATT","GATTGCGAG","GATTCCAGT","GATTCATAT","GATTATCAG","GATTAGGTT",
                  "GATGTTGCG","GATGGATCT","GATGCTGAT","GATGCCTTG","GATCTCCTT","GATCGCTTA","GATATTGAA","GATATTACT","GAGTGTTAT","GAGCTCAGT",
                  "GAGCGTGCT","GAGCGTCGA","GAGCGGTTG","GAGCGACTT","GAGCCGAAT","GAGATAGAT","GAGACCTAT","GACGGTCGT","GACGCAGGT","GACGATATG",
                  "GACCTATCT","GAATTAGGA","GAATCAGCT","GAAGTTCAT","GAAGTGGTT","GAAGTATTG","GAAGGCATT","GAACGCTGT","CTTGTCCAG","CTTGGATTG",
                  "CTTGCTGAA","CTTGCCGTG","CTTGATTCT","CTTCTGTCG","CTTCGGCGT","CTTATGAGT","CTTACCGAT","CTGTTAGGT","CTGTCGTCT","CTGTATAAT",
                  "CTGGCTCAT","CTGGATGCG","CTGCGTGTG","CTGCGCGGT","CTGCCGATT","CTGCATTGT","CTGATTAAG","CTGAGATAT","CTGACCTGT","CTCGTATCT",
                  "CTCGGCAAG","CTCGCAATT","CTCCTGCTT","CTCCTAAGT","CTCCGGATG","CTCCGAGCG","CTCACAGGT","CTATTCTAT","CTATTAGTG","CTATGAATT",
                  "CTACATATT","CGTGGCATT","CGTCTTAAT","CGTCTGGTT","CGTCACTGT","CGTAGGTCT","CGGTTCGAG","CGGTTCATT","CGGTGCTCT","CGGTAATTG",
                  "CGGCCTGAT","CGGATATAG","CGGAATATT","CGCTCCAAT","CGCGTTCGT","CGCAGGTTG","CGAGGATGT","CGAGCTGTT","CGACGGCTT","CCTTGTGTG",
                  "CCTGTCTCA","CCTGACTAT","CCTACCTTG","CCGTAGATT","CCGGCTGGT","CATCGGACG","CATCGATAA","CATCCTTCT","CAGTTCTGT","CAGTGCCAG",

                  "CAGGCACTG","CAGCCTCTT","CACTTATAT","CACTGGTCG","CACTGCATG","CACGCGTTG","CACGATGTT","CACCATCTG","CACAGGCGT","ATTGTACAA",
                  "ATTGGTATG","ATTGCTAAT","ATTGCATAG","ATTGCAGTT","ATTCTGCAG","ATTCTACGT","ATTCGGATT","ATTCCGTTG","ATTCATCAA","ATTCAAGAG",
                  "ATTAGCCTT","ATTAATATT","ATGTTAGAG","ATGTTAACT","ATGTAGTCG","ATGGTGTAG","ATGGATTAT","ATCTTGAAG","ATCTGATAT","ATCTCAGAA",
                  "ATCGCTCAA","ATCGCGTCG","ATCCATGGT","ATCATGAGA","ATCATAGTT","ATCAGCGAG","ATCACCATT","ATAGTAATT","ATAGCTGTG","ATACTCTCG",
                  "ATACCTCAT","AGTTGCGCG","AGTTGAATT","AGTTATGAT","AGTGTCCGT","AGTGGCTTG","AGTGCTTCT","AGTATCATT","AGTACACAA","AGGTATGCG",
                  "AGGTATAGT","AGGCTACTT","AGGCCAGGT","AGGAGCGAT","AGCTTATAG","AGCTCTAGA","AGCGTGTAT","AGCGTCACA","AGCCTTCAT","AGCCTGTCG",
                  "AGCCTCGAG","AGCACTGAA","AGATGTACG","AGAGTTAAT","AGACCTCTG","ACTTCTATA","ACTGTCGAG","ACTGTATGT","ACTCTGTAA","ACTCGCGAA",
                  "ACTAGATCT","ACTAACGTT","ACGTTACTG","ACGTGGAAT","ACGGACTCT","ACGCCTAAT","ACGCCGTTA","ACGACGTGT","ACCTCGCAT","ACCATCATA",
                  "ACATATATT","ACAGGCACA","ACACCTGAG","ACACATTCT"
                  )

B384_cell_key3 = ("TTGTGGCTG","TTGTGGAGT","TTGTGCGAC","TTGTCTTCA","TTGTAAGAT","TTGGTTCTG","TTGGTGCGT","TTGGTCTAC","TTGGTAACT","TTGGCGTGC",
                  "TTGGATTAG","TTGGAGACG","TTGGAATCA","TTGCGGCGA","TTGCGCTCG","TTGCCTTAC","TTGCCGGAT","TTGCATGCT","TTGCACGTC","TTGCACCAT",
                  "TTGAACCTG","TTCTCGCGT","TTCTCAACT","TTCTACTCA","TTCGTCCAT","TTCGGATAC","TTCGGACGT","TTCGCAATC","TTCCGGTGC","TTCCGACTG",
                  "TTCATTATG","TTCATGGAT","TTCAGCGCA","TTCACCTCG","TTCAAGCAG","TTCAACTAC","TTATGCCAG","TTATGCATC","TTATCGTAC","TTATACCTA",
                  "TTATAATAG","TTATAAGTC","TTAGTTAGC","TTAGCTCAT","TTAGCACTA","TTAGATATG","TTACTACGA","TTACCGTCA","TTACAGAGC","TTAATTGCA",
                  "TTAACAGAT","TGTTGGCTA","TGTTGATGA","TGTTAAGCT","TGTGGCCGA","TGTGCTAGC","TGTGCGTCA","TGTCGCAGT","TGTCGAGCA","TGTACAACG",
                  "TGGTTCCGA","TGGTTCACT","TGGTCAAGT","TGGCTTGTA","TGGCTGTCG","TGGCGTATG","TGGCGCGCT","TGGATGTAC","TGGACTTGC","TGGAATACT",
                  "TGCTAGCGA","TGCGTTGCT","TGCGGTCTG","TGCGCTTAG","TGCGCGACG","TGCCTGCAT","TGCCTAGAC","TGCACGAGT","TGAGTGTGC","TGAGGCTCG",
                  "TCTTCCGTC","TCTTATAGT","TCTTACCAT","TCTGTTGTC","TCTGTTACT","TCTGGCTAG","TCTCAGATC","TCTAGTTGA","TCTAGTACG","TCGTACTAC",
                  "TCGGTGTAG","TCGGCTGCT","TCGCTACTG","TCGATCACG","TCGAGGCAT","TCCGGCGTC","TCCGGAGCT","TCCGCTCGT","TCCGAGTAC","TCCATTCAT",

                  "TCCATGGTC","TCCAAGTCG","TCATTACGT","TCATGCACT","TCAGGTTGC","TCAGACCGT","TCACTCAGT","TCAAGCTCA","TATTGCGCA","TATTCGGCT",
                  "TATTCCAGC","TATTCATCA","TATGTTCAG","TATGGTATG","TATGCAAGT","TATCTGGTC","TATCTGACT","TATCCAGAT","TATCAGTCG","TATCACGCT",
                  "TAGGCGCGA","TAGGCACAT","TAGGATCGT","TAGCATTGC","TAGAGTTAC","TAGACTGAT","TACTTGTCG","TACGTCCGA","TACCGTACT","TACCGCGAT",
                  "TACCAGGAC","TACAGAAGT","TAAGTGCAT","TAAGCTACT","GTTGACCGA","GTTCTCGAC","GTTCCTGCT","GTTATGATG","GTGCTTGCA","GTGCCGCGT",
                  "GTATTGCTG","GTATTCCGA","GTATTAAGC","GTATGACGT","GTAGTTGTC","GTAGTACAT","GTAGCTCGA","GGTTGCTCA","GGTTGAGTA","GGTTAACGT",
                  "GGTGTGGCA","GGTCTTCAG","GGTCGTCTA","GGTCGGCGT","GGTCCGACT","GGTCATGTC","GGTCACATG","GGTAGTGCT","GGTAGCGTC","GGTACCAGT",
                  "GGTAAGGAT","GGCTTGTGC","GGCTTGACT","GGCTTACGA","GGCTGTAGT","GGCTGGCAG","GGCTCCATC","GGCGTGGAT","GGCGTAATC","GGCGCAAGT",
                  "GGCGAGTAG","GGCGACCGT","GGCCTGTCA","GGCCATTGC","GGCACTCTG","GGATGTCAT","GGAGTAACT","GGAGAACGA","GGACTGGCT","GGACGTTCA",
                  "GGAACGTGC","GCTGTCCAT","GCTGGTTCA","GCTGCAACT","GCTCGTTAC","GCTATAGAT","GCTAGTCGT","GCTACCATG","GCGTTCTGA","GCGTGTTAG",
                  "GCGGTATCG","GCGGAGCAT","GCGCGGTGC","GCGCCTAGT","GCGCCGGCT","GCCTTCATG","GCCATACTG","GCATGTTGA","GCATGCTAC","GCAGTATAC",

                  "GCAGGTACT","GCAGCGCGT","GCACCTCAT","GCAATTCGA","GATTGCCGT","GATGAACAT","GATCTTCGA","GATCTGCAT","GAGTGGCAT","GAGTCGGAC",
                  "GAGTATGAT","GAGGCGAGT","GAGGCAACG","GAGCGCACT","GAATAGGCT","ATTGTCACT","ATTGTATCA","ATTGGTCAG","ATTGGCGAT","ATTGATCGT",
                  "ATTCGTAGT","ATTCATACG","ATTCAGGAC","ATTACTTCA","ATTAATTAG","ATTAAGCAT","ATGTCTCTA","ATGTAGCGT","ATGGCATAC","ATGGAGATC",
                  "ATGGACTCG","ATGGAACGA","ATGCTTCAT","ATGCTCGCT","ATGCGACGT","ATGCCGTAG","ATGAGTTCG","ATGACTATC","ATGACCGAC","ATCTTATGC",
                  "ATCTTACTA","ATCTATCAG","ATCGTGTAC","ATCGTCTGA","ATCGGCATG","ATCGCGAGC","ATCGCAACG","ATCGATGCT","ATCGAATAG","ATCCTTCTG",
                  "ATCCTGCGT","ATCCGCACT","ATCCATTAC","ATCCAAGCA","ATCAGATCA","ATCACACAT","ATCAACGTC","ATCAACCGA","ATATTGAGT","ATATTCGTC",
                  "ATATTACAG","ATATCTTGA","ATATCGCAT","ATATCAATC","ATAGTCCTG","ATAGGTCTA","ATAGCTGAC","ATAGCGGTA","AGTTCGCTG","AGTTACAGC",
                  "AGTTAACTA","AGTGCAATC","AGTCTGGTA","AGTCTGAGC","AGTCTACAT","AGTCGAACT","AGTCCATCG","AGTCATTCA","AGTATCCAG","AGTAGACTG",
                  "AGTAATCGA","AGTAAGTGC","AGGTTGGCT","AGGTTCTAG","AGGTGTTCA","AGGTGCCAT","AGGTCTGAT","AGGTCGTAC","AGGTCAGCA","AGGCTTATC",
                  "AGGCTATGA","AGGCCGACG","AGGCCAAGC","AGGCAGGTC","AGGCAAGAT","AGGAGCAGT","AGGACCGCT","AGGAATTAC","AGCTTGGAC","AGCTTAAGT",

                  "AGCTACACG","AGCGTTACG","AGCGGTGCA","AGCGGAGTC","AGCGGACGA","AGCGCGCTA","AGCGATAGC","AGCGACTCA","AGCCTCTAC","AGCCGTCGT",
                  "AGCATGATC","AGCACTTCG","AGCACGGCA","AGATTCTGA","AGATTAGAT","AGATGATAG","AGATATGTA","AGATACCGT","AGAGTGCGT","AGAGCCGAT",
                  "AGACTCACT","ACTTGCCTA","ACTTGAGCA","ACTTCTAGC","ACTTCGACT","ACTTAGTAC","ACTGTTGAT","ACTGTAACG","ACTGGTATC","ACTGACGTC",
                  "ACTGAAGCT","ACTCTGATG","ACTCCTGAC","ACTCCGCTA","ACTCAACTG","ACTATTGCA","ACTAGGCAG","ACTACGCGT","ACTAATACT","ACGTTCGTA",
                  "ACGTGTGCT","ACGTGTATG","ACGTGGAGC","ACGTCTTCG","ACGTCAGTC","ACGGTCTCA","ACGGTCCGT","ACGGTACAG","ACGGCGCTG","ACGCTGCGA",
                  "ACGCGTGTA","ACGCGCCAG","ACGATGTCG","ACGATGGAT","ACGATCTAC","ACGAGCTGA","ACGAGCATC","ACGAATCGT","ACGAACGCA","ACCTTGTAG",
                  "ACCTGTTGC","ACCTGTCAT","ACCTCGATC","ACCTAGGTA","ACCTACTGA","ACCTAATCG","ACCGTAGCA","ACCGGTAGT","ACCGGCTAC","ACCGCTTCA",
                  "ACATTGTGC","ACATTCTCG","ACATGGCTG","ACATGACGA","ACATATGAT","ACATATACG","ACAGCGTAC","ACACTTGCT","ACACTATCA","ACACGCATG",
                  "ACACCAGTA","ACACCAACT","ACACATAGT","ACACACCTA"
                  )


def label_sections_to_index(label):
    """ 
    Return the cell_index integer based on input 3 part cell label string
    
    """

    cl1, cl2, cl3 = [int(n) for n in label.split('-')]
    return (cl1 - 1) * 384 * 384 + (cl2 - 1) * 384 + (cl3 - 1) + 1


# print(label_sections_to_index('1-1-1'))
# print(label_sections_to_index('33-78-21'))
# print(label_sections_to_index('43-12-77'))
# print(label_sections_to_index('96-96-96'))
# print(label_sections_to_index('135-43-344'))
# print(label_sections_to_index('384-384-384'))
# print('-')

#----------------------------------


def index_to_label_sections(index):

    zerobased = int(index) - 1

    cl1 = (int((zerobased) / 384 / 384) % 384) + 1
    cl2 = (int((zerobased) / 384) % 384) + 1
    cl3 = (zerobased % 384) + 1

    return f'{cl1}-{cl2}-{cl3}'


# print(index_to_label_sections(1))
# print(index_to_label_sections(4748181))
# print(index_to_label_sections(6197453))
# print(index_to_label_sections(14044896))
# print(index_to_label_sections(19775576))
# print(index_to_label_sections(56623104))
# print('-')
#----------------------------------


def index_to_sequence(index, bead_version):

    zerobased = int(index) - 1

    cl1 = (int((zerobased) / 384 / 384) % 384) + 1
    cl2 = (int((zerobased) / 384) % 384) + 1
    cl3 = (zerobased % 384) + 1

    if bead_version == 'v1':
        cls1_sequence = A96_cell_key1[cl1-1]
        cls2_sequence = A96_cell_key2[cl2-1]
        cls3_sequence = A96_cell_key3[cl3-1]

        return f'{cls1_sequence}{v1_linker1}{cls2_sequence}{v1_linker2}{cls3_sequence}'

    elif bead_version == 'Enh':

        diversityInsert = ''

        if 1 <= cl1 <= 24:
            diversityInsert = ''
        elif 25 <= cl1 <= 48:
            diversityInsert = 'A'
        elif 49 <= cl1 <= 72:
            diversityInsert = 'GT'
        else: # 73 <= cl1 <= 96:
            diversityInsert = 'TCA'

        cls1_sequence = A96_cell_key1[cl1-1]
        cls2_sequence = A96_cell_key2[cl2-1]
        cls3_sequence = A96_cell_key3[cl3-1]

        return f'{diversityInsert}{cls1_sequence}{Enh_linker1}{cls2_sequence}{Enh_linker2}{cls3_sequence}'

    elif bead_version == 'EnhV2':

        diversityInsert = ''
        subIndex = ((cl1-1) % 96) + 1

        if 1 <= subIndex <= 24:
            diversityInsert = ''
        elif 25 <= subIndex <= 48:
            diversityInsert = 'A'
        elif 49 <= subIndex <= 72:
            diversityInsert = 'GT'
        else: # 73 <= subIndex <= 96:
            diversityInsert = 'TCA'

        cls1_sequence = B384_cell_key1[cl1-1]
        cls2_sequence = B384_cell_key2[cl2-1]
        cls3_sequence = B384_cell_key3[cl3-1]

        return f'{diversityInsert}{cls1_sequence}{Enh_linker1}{cls2_sequence}{Enh_linker2}{cls3_sequence}'


# print(index_to_sequence(4748181, 'Enh'))
# print(index_to_sequence(52923177, 'EnhV2'))

#----------------------------------


def create_cell_index_fasta_V1():
    with open('Rhapsody_cellBarcodeV1_IndexToSequence.fasta', 'w') as f:
        for cl1 in range(1, 96+1):
            for cl2 in range(1, 96+1):
                for cl3 in range(1, 96+1):
                    index = label_sections_to_index(f'{cl1}-{cl2}-{cl3}')
                    sequence = index_to_sequence(index, 'v1')
                    f.write(f'>{index}\n')
                    f.write(f'{sequence}\n')

#create_cell_index_fasta_V1()


def create_cell_index_fasta_Enh():
    with open('Rhapsody_cellBarcodeEnh_IndexToSequence.fasta', 'w') as f:
        for cl1 in range(1, 96+1):
            for cl2 in range(1, 96+1):
                for cl3 in range(1, 96+1):
                    index = label_sections_to_index(f'{cl1}-{cl2}-{cl3}')
                    sequence = index_to_sequence(index, 'Enh')
                    f.write(f'>{index}\n')
                    f.write(f'{sequence}\n')

#create_cell_index_fasta_Enh()

def create_cell_index_fasta_EnhV2():
    with open('Rhapsody_cellBarcodeEnhV2_IndexToSequence.fasta', 'w') as f:
        for cl1 in range(1, 384+1):
            for cl2 in range(1, 384+1):
                for cl3 in range(1, 384+1):
                    index = label_sections_to_index(f'{cl1}-{cl2}-{cl3}')
                    sequence = index_to_sequence(index, 'EnhV2')
                    f.write(f'>{index}\n')
                    f.write(f'{sequence}\n')

#create_cell_index_fasta_EnhV2()
