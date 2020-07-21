import re,sys,argparse
from collections import Counter

"""
FOREST.py ver.1 
Written by Kaoru R. Komatsu

"""
parser = argparse.ArgumentParser(
    prog = "FOREST",
    usage = "RNA motif extraction and library design",
    description=u"Output file: Multiple FASTA file of the extracted RNA motifs.",
    add_help=True)

parser.add_argument('filename',help=u"Input file: Multiple FASTA file with RNA secondary structure in dot-bracket format.")

parser.add_argument("-L", type=int, default=100,help=u"This option limits the maximum length of the extracted motif (Default: 100).")
parser.add_argument("-bn", type=int, default=5, help=u"Barcode-number")
parser.add_argument("-b","--barcodes", help=u"Barcode-number")
parser.add_argument("-lib","--library",  action="store_true",help=u"The attachements above")
parser.add_argument("-t","--templates", action="store_true",help=u"DNAtemplate")
parser.add_argument("-a","--array", action="store_true",help=u"DNAtemplate")

args = parser.parse_args()
p = re.compile('\(+\.+\)+')

def bracket_divider(seq):
    i = 0
    while seq[i] == "." or seq[i] ==")" :
        i += 1
        if i == len(seq):
            break
    return [seq[0:i],seq[i:]]

def revcom(seq):
    seq = seq.upper()
    seq = seq.replace("A","t")
    seq = seq.replace("U","a")
    seq = seq.replace("T","a")
    seq = seq.replace("G","c")
    seq = seq.replace("C","g")
    return seq[::-1].upper()

def hitpoint_brew(hitpoint,seq):
    i = 0
    while hitpoint != 0:
        if seq[i] == ")" or seq[i] == "(":
            hitpoint -= 1
        i += 1
    return seq[0:i]

def hitpoint_counter(hitpoint,seq,bracket):
    partition = ""
    i = 0
    if bracket == ")":
        while hitpoint != partition.count(bracket):
            partition = seq[0:i]
            i += 1
        return partition
    elif bracket == "(":
        seq_rev = seq[::-1]
        while hitpoint != partition.count(bracket):
            partition = seq_rev[0:i]
            i += 1
        return partition[::-1]

def loop_brew(seq):
    output = []
    main = p.findall(seq) 
    part = p.split(seq)

    for m in range(len(main)):

        n = m+1
        left = part[m]
        right = part[n]

        if left == "":
            left = "."
        if right == "":
            right = "."

        part_left = bracket_divider(left)[1] 
        part_right = bracket_divider(right)[0] 
        whole = ''.join([part_left,main[m],part_right])
        hitpoint_right_g = hitpoint_left_g = min(whole.count("("),whole.count(")"))
        
        main_c = Counter(main[m])
          
        hitpoint_left = hitpoint_left_g - int(main_c["("])
        hitpoint_right = hitpoint_right_g - int(main_c[')'])

        if hitpoint_right == 0 and hitpoint_left == 0:
            output.append(main[m].strip('.'))
            
        elif hitpoint_right < 0 or hitpoint_left < 0:
            refseq = main[m]

            if hitpoint_right < 0:
                refseq = hitpoint_counter(hitpoint_right_g,refseq,")")
                hitright_left = hitpoint_counter(hitpoint_left,left,"(")
                final_out = (''.join([hitright_left,refseq]))

            elif hitpoint_left < 0:
                refseq = hitpoint_counter(hitpoint_left_g,refseq,"(")
                hitright_left = hitpoint_counter(hitpoint_right,right,")")
                final_out = (''.join([refseq,hitright_left]))


            elif hitpoint_right < 0 and hitpoint_left < 0:
                print("Error")
                break

            output.append(final_out)

        else:
            final_right =  hitpoint_counter(hitpoint_right,part_right,")")
            final_left =  hitpoint_counter(hitpoint_left,part_left,"(")
            output.append(''.join([final_left,main[m],final_right]))
    return output


def dot_replace(seq,looplist):
    for i in looplist:
        num = len(i)
        seq = seq.replace(i,"."*num,1)
    return seq

def packman(seq):
    i = 0
    t = 0
    while seq[i] != ".":
        i += 1
    while seq[len(seq)-1-t] != ".":
        t += 1
    pack_num = max(i,t)
    i = 0
    t = 0
    while seq[:i].count("(") != pack_num:
        i +=1
    while seq[len(seq)-1-t:len(seq)].count(")") != pack_num:
        t +=1
    return seq[i:len(seq)-1-t].strip('.')

def loop_list_judge(Secondary_loop_list):

    if len(Secondary_loop_list) >= 1:
        return "Read"
    elif Secondary_loop_list == []:
        return "Zero"
    else:
        return "Error"

def locus_find(seq,loop):
    left = seq.index(loop)
    right = left + len(loop)
    return [left,right]

def theshold(secondary_seq,loop_structure,limit):
    conserved_area = 0
    output = []
    if len(loop_structure) > limit:
        while len(loop_structure) >= limit:
            loop_structure = packman(loop_structure)
    locus = locus_find(secondary_seq,loop_structure)
    left = locus[0]
    right = locus[1]
    return [left,right]

def loop_replace_extract(seq,seq2,single_structure_result):
    output = []
    for i in single_structure_result:
        left = seq.find(i)
        right = left + len(i)
        rep = '!' * len(i)
        output.append([seq[left:right],seq2[left:right]])
        seq = seq.replace(str(i),rep,1)
    return output

def conjugateStem(seq,n=17):
    forward = (''.join(list("GTGTACGAAGTTTCAGC")[0:n]))
    reverse = (''.join(list("GCTGAAGCTTCGTGCAC")[:n]))
    output = forward + seq + reverse
    return output 

 
LibraryLengthlimit = args.L
def terminalmotif_extraction(name,seq,seq2):
    output = {}

    ##############################################################
    # Single terminal motif extraction#
    # ############################################################

    loop_list =  loop_brew(seq)
    result1 = loop_replace_extract(seq,seq2,loop_list)

    strct = seq
    motifcount = 1
    for t in range(len(result1)):
        identical_name = '_'.join([name,"Motif",str(motifcount)])
        
        LRlocus = theshold(strct,result1[t][0],limit = args.L)
        outseqandstrc = [strct[LRlocus[0]:LRlocus[1]],seq2[LRlocus[0]:LRlocus[1]]]

        if len(outseqandstrc[1]) <= args.L:
            output[identical_name] = outseqandstrc
            motifcount +=1
        rep = '!' * len(result1[t][0])
        strct = strct.replace(str(result1[t][0]),rep,1)

    ##############################################################
    # Multiple terminal motif extraction#
    ##############################################################
    
    def multiterminalmotif_extraction(name,secondary_seq,Secondary_loop_list, LibraryLengthlimit,replacetimes):
        output2 = {}
        if loop_list_judge(Secondary_loop_list) == "Read":
            for k,i in enumerate(Secondary_loop_list):
                loop_rep = i
                identical_name_multi = '_'.join([name,"Multi",str(k+1),"ComplexLevel",str(1+replacetimes)])
                LRlocus = theshold(secondary_seq,loop_rep,LibraryLengthlimit)
                outseqandstrc = [seq[LRlocus[0]:LRlocus[1]],seq2[LRlocus[0]:LRlocus[1]]]
                if len(outseqandstrc[1]) <= LibraryLengthlimit:
                    output2[identical_name_multi] = outseqandstrc
        return output2

    secondary_seq = dot_replace(seq, loop_list)

    replacetimes = 0
    while "("  in secondary_seq:
        Secondary_loop_list = loop_brew(secondary_seq)
        output.update(multiterminalmotif_extraction(name, secondary_seq,Secondary_loop_list,args.L,replacetimes))
        secondary_seq = dot_replace(secondary_seq,Secondary_loop_list)
        replacetimes +=1
    return output


outdict = {}

f = open(args.filename)
for line in f:
        line = line.strip('\n')
        if line.startswith('>'):
            #header#
            tokens = line.split('|')
            name = tokens[0]
        elif line[0].lower() in ["a","t","g","c"]:
            #Sequence
            seq = line
        elif line.startswith('(') or line.startswith('.') or line.startswith(')'):
            #Structure
            tokens3 = line.split(' ')
            structure = tokens3[0]
            structure = structure.replace("&",".")
            stru = structure.replace("]",".").replace("[",".") #Remove the PseudoKnots.
            extractedmotif_info = terminalmotif_extraction(name,stru,seq)
            outdict.update(extractedmotif_info)
f.close()



if args.barcodes:
    filebar = open(args.barcodes)
    barcodelist = []
    for line in filebar:
        if line[0].lower() in ["a","t","g","c"]:
            line = line.strip('\n')
            barcodelist.append(line)
    filebar.close()

    print('# Loading %d barcodes' % len(barcodelist))


if args.library:
    #Add barcodes
    outdict_barcodes = {}
    barcode_list_num = 0
    used_barcode_list = []

    uniquelist = [sq[1] for (st,sq) in outdict.items() if sq[1] != ""]

    print("# The number of barcodes per each RNA structure: %d" % args.bn)
    print("# The number of RNA structures: %d" % len(uniquelist))
    if len(uniquelist)*args.bn > len(barcodelist):
        print("Error: You need more %d barcodes" % (len(uniquelist)*args.bn - len(barcodelist)))
    else:

        for i,(k,v) in enumerate(outdict.items()):
            #print(i)#Number
            #print(k)#Name
            #print(v)#[Strucure,Sequence]
            for barcode_id in range(1,args.bn+1):
                probe_id = k +  "_Barcode_" + str(barcode_id)
                e_barcode = str(barcodelist[barcode_list_num])
                outdict_barcodes[probe_id] = [v[0],"GGG" + e_barcode + conjugateStem(v[1])]
                barcode_list_num +=1
                used_barcode_list.append(e_barcode)

         
if args.library:
    for i,(k,v) in enumerate(outdict_barcodes.items()):
        if v[1] != "" and v[0] != "":
            if args.templates:
                print(k + "_template")
                output_bs = "GCGCTAATACGACTCACTATA" + v[1]
                print(revcom(output_bs.upper().replace("U","T")))
            elif args.array:
                print(k+"_array")
                print(revcom("GGG"+used_barcode_list[i]))

            else:
                print(k)
                print(v[1].upper().replace("T","U"))

else:
    for k,v in outdict.items():
        if v[1] != "" and v[0] != "":
            print(k)
            print(v[1])
            print(v[0])
