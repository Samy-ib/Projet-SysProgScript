from Bio import SeqIO
from collections import Counter
import random
import numpy as np



#=============================================#
# Ce fichier contient les fonctions servant a #
# manipuler nos différents type de séquence   #
#=============================================#

def generate():
    """
        Generation d'une séquence ADN aléatoire de max 1000 nucléotide.
    """
    seq=''
    for i in range(random.randint(1,1000)):
        seq+=random.choice('ATCG')
    return seq

def valide(seq):
    """
        Vérifie la validité d'une séquence ADN.
    """
    return len(seq)==(seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T'))

def adn_to_arn(seq):
    """
       Transforme une séquence ADN en ARN. 
    """
    return seq.replace("T", "U")

def adn_to_proteine(seq):
    """
       Transforme une séquence adn en protéine après avoir vérifier 
       que la taille de cette dernière est divisible par trois.
       Si elle ne l'est pas les dernier élément seront supprimer pour qu'elle 
       puisse être convertis en protéine.
    """
    if len(seq)%3 != 0 :
        seq = seq[:-(len(seq)%3)]
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    for i in range(0, len(seq), 3): 
        codon = seq[i:i + 3] 
        protein+= table[codon] 
    return protein 

def arn_to_proteine(seq):
    """
       Tranforme une séquence ARN en proteine 
       (même raisonement que "adn_ro_protein()") 
    """
    if len(seq)%3 != 0 :
        seq = seq[:-(len(seq)%3)]
    table = {
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C', # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '_', 'UGA': '_', # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '_', 'UGG': 'W', # UxG

    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', # CxG

    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', # AxG

    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    protein ="" 
    for i in range(0, len(seq), 3): 
        codon = seq[i:i + 3] 
        protein+= table[codon] 
    return protein 

def calcul_freq(seq):
    """
       Calcul la fréquence de chaque nucléotide de la séquence ADN .
    """
    return  "Le nombre d'Adenine dans la chaine est :" + str(seq.count('A')) + " (" + str("%.2f" % (seq.count('A')/len(seq)*100)) + "%" + ")" + "\n" + \
            "Le nombre de Cytosine dans la chaine est :" + str(seq.count('C')) + " (" + str("%.2f" % (seq.count('C')/len(seq)*100)) + "%" + ")" + "\n" + \
            "Le nombre de Thymine dans la chaine est :" + str(seq.count('T')) + " (" + str("%.2f" % (seq.count('T')/len(seq)*100)) + "%" + ")"  + "\n" + \
            "Le nombre de Guanine dans la chaine est :" + str(seq.count('G')) + " (" + str("%.2f" % (seq.count('G')/len(seq)*100)) + "%" + ")" 

def comp_inv(seq):
    """
       Complément inverse. 
    """
    for a in range(len(seq)):
        if seq[a]=='A':
            seq[a]='T'
        elif seq[a]=='T':
            seq[a]='A'
        elif seq[a]=='C':
            seq[a]='G'
        elif seq[a]=='G':
            seq[a]='C'
    return "".join(seq[::-1])

def taux_gc(seq):
    """
       Le taux de "GC" dans la séquence ADN. 
    """
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100

def freq_codon(seq):
    """
       Fréquence des codons.
       On transforme d'abbord notre séquence en protéine.
       Cette liste est ensuite donner en paramètre a la fonction 
       "Counter" du module "collection" qui nous renvois un dictionnaire
       de la fréquence de chaque codon.
       Ce dernier est transformer en liste avant d'être transformer en une chaine
       pour être afficher.
    """
    table = {
            'F': 'Phenylalanine',       'I': 'Isoleucine',      'V': 'Valine',      'L': 'Leucine',         
            'M': 'Methionine',          'S': 'Serine',          'P': 'Proline',     'T': 'Threonine',
            'A': 'Alanine',             'Y': 'Tyrosine',        'H': 'Histidine',   'N': 'Asparagine',     
            'D': 'Acide Aspartique',    '_': 'Codon Stop',      'Q': 'Glutamine',   'K': 'Lysine',     
            'E': 'Acide Glutamique',    'C': 'Cysteine',        'R': 'Arginine',    'G': 'Glycine',
            'W': 'Tryptophane'
        }
    if len(seq)%3 != 0 :
        seq = seq[:-(len(seq)%3)]
    freq = ""
    seq=adn_to_proteine(seq)
    count=dict((Counter(list(seq))))
    l=[]
    [l.extend([k,v]) for k,v in count.items()]
    for i in range(0,len(l), 2):
        l[i] = table[l[i]]
    for x in range(0,len(l),2):
        freq += l[x] + " : " + str(l[x+1]) +"\n"
    return freq

def masse_proteique(seq):
    """
       Calcul de la masse protéique a l'aide
       d'un dictionnaire contenant la masse 
       chaque protéine.
    """
    masses = { 
        'A':71.03711, 'C':103.00919, 'D':115.02694, 'E':129.04259, 
        'F':147.06841, 'G':57.02146, 'H':137.05891, 'I':113.08406, 
        'K':128.09496, 'L':113.08406, 'M':131.04049, 'N':114.04293, 
        'P':97.05276, 'Q':128.05858, 'R':156.10111, 'S':87.03203,                  
        'T':101.04768, 'V':99.06841, 'W':99.06841, 'Y':163.06333,
        '_':0,
    } 
    prot = adn_to_proteine(seq)
    return sum(masses[a] for a in prot)

def epissage_adn(seq, intr1, intr2):
    seq = seq.replace(intr1,"").replace(intr2,"")
    return arn_to_proteine(seq)
    
def assem(seq, taille):
    list_seq = [seq[start:start+taille] for start in range(0, len(seq), taille)]
    if (len(list_seq[-1]) != len(list_seq[0])):
        list_seq.pop()
    return search(list_seq)

def my_overlap(s1,s2): 
    """
        Reçoit deux chaines de caractères et nous retourne la plus petite
        chaine contenant s1 et s2. ainsi que le nombre de nombre de caracteres
        communs entre les deux.
    """
    if s2 in s1:
        return [s1, len(s2)]
    elif s1 in s2:
        return [s2, len(s1)]

    for i in range(1, len(s2)):
        if s1.startswith(s2[i:]):
            return [s2[:i] + s1,i]
        if s1.endswith(s2[:-i]):
            return [s1 + s2[-i:],i]
    return [s1+s2,0]

def my_most_overlap(seqs):                            # SEQ NAME, OVERLAPPING STRINGS, SEQ1 ID, SEQ2 ID
    """
        On lui donne une liste de chaine de caracteres(Liste de sequences)
        elle nous retourne la position dans cette liste des 2 première chaînes
        avec le plus de sous chaine comune (most overlaping strings).
    """
    L=[]
    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if j!=i :
                ovl=my_overlap(seqs[i],seqs[j])
                ovl.extend([i,j])
                L.append(ovl)
    L = np.array(L)
    li=np.argmax(L[:,1],axis=0)
    return L[li,2:].tolist()

def search(seqs):
    """
        Concatène les deux chaine les plus "Overlapper" a chaque itterations
        tout en supprimant une des deux chaînes (La seconde reçoit le résultat de la concaténation).
        Ceci s'arrête quand la longueur de la liste est égale a 1.
    """
    while(len(seqs)>1):
        ov=my_most_overlap(seqs)
        seqs[int(ov[0])]=my_overlap(seqs[int(ov[0])], seqs[int(ov[1])])[0]
        del seqs[int(ov[1])]
    return seqs[0]