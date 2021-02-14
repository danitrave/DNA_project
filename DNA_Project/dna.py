rna2amino = {'UUU':'F','UAC':'Y','UUA':'L','UUG':'L','CUU':'L','CUC':'L','CUA':'L','CUG':'L','AUU':'I','AUC':'I','AUA':'I','AUG':'M','GUU':'V','GUC':'V','GUA':'V','GUG':'V','UCU':'S','UCC':'S','UCA':'S','UCG':'S','CCU':'P','CCC':'P','CCA':'P','CCG':'P','ACU':'T','ACC':'T','ACA':'T','ACG':'T','GCU':'A','GCC':'A','GCA':'A','GCG':'A','UAU':'Y','CAU':'H','CAC':'H','CAA':'Q','CAG':'Q','AAU':'N','AAC':'N','AAA':'K','AAG':'K','GAU':'D','GAC':'D','GAA':'E','GAG':'E','UGU':'C','UGC':'C','UGG':'W','CGU':'R','CGC':'R','CGA':'R','CGG':'R','AGU':'S','AGC':'S','AGA':'R','AGG':'R','GGU':'G','GGC':'G','GGA':'G','GGG':'G'}
amb_symb = { "R":"AG", "Y":"CT", "K":"GT", "M":"AC", "S":"CG", "W":"AT", "B":"CGT", "D":"AGT", "H":"ACT", "V":"ACG", "N":"ACGT" }
dinucleotides = ("AA","AT","AC","AG","TT","TA","TC","TG","CC","CG","CT","CA","GG","GC","GA","GT")
start_stop_codons = {"ATG":">>>","TAA":"<<<","TAG":"<<<","TGA":"<<<"}
import random
import turtle
import urllib.request
tortoise = turtle.Turtle()


def gcContent(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    gc = 0
    
    for e in dna:
        if e == "G" or e == "C":
            gc += 1
            
    return (gc/len(dna))

    """Returns the GC content ratio of a DNA sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a real number between 0 and 1
    Example: gcContent("atcgttcaag") = 0.4
    """


def countCodon(dna, codon):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    codon = codon.upper()
    count = 0
    
    for e in range(0, len(dna)-1, 3):
        if dna[e:e+3] == codon:
               count += 1
               
    return count

    """Returns the number of (non-overlapping) occurrences of a codon in a DNA sequence
    Parameters:
        dna: a string object representing a DNA sequence
        codon: a three-letter string object representing the codon to
                search for
    Return value: the integer number of instances of the target codon in dna
    Example: countCodon("aaaaaaaa", "aaa") = 2
    """

def countACG(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    s = ""
    
    for e in dna:
        if e != "T":
            s += e
            
    return len(s)

    """Returns the number of nucleotides that are not T in a DNA sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: the integer number of nucleotides in dna that are not T
    Example: countACG("atcgttcaag") = 7
    """


def printCodons(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    
    Fr1 = tuple()
    sFr1 = ""
 
    for e in range(0, len(dna), 3):
         Fr1 += (dna[e:e+3],)

    for e in Fr1:
         sFr1 += str(e) + " "

    return print(sFr1)

    """Prints the sequence of non-overlapping codons in dna, 
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: None
    Example: printCodons("ggtacactgta") would print: ggt aca ctg
    """


def findCodon(dna, codon):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    codon = codon.upper()
    div_dna = ""

    for i in range(0,len(dna),3):
        div_dna += dna[i:i+3] 

    if div_dna.find(codon) != -1:
        return div_dna.find(codon)
    else:
        return None

    """Returns the index of the first occurrence of codon in dna
    Parameters:
        dna: a string object representing a DNA sequence
        codon: a three-letter string object representing the codon to
                search for
    Return value:
        (if codon found): the integer index of the first occurrence of codon in dna
        (if codon not found): None
    Example: findCodon("ggtacactacgta", "tac") = 2 
    """


def findATG(dna):

    if dna == "":
        return None
    
    dna = dna.upper()
    index_ATG = []
    i = 0

    if dna.find("ATG") == -1:
        return 
    else:
        while dna.find("ATG",i) != -1:
            index_ATG.append(dna.find("ATG",i))
            i = dna.find("ATG", i)+1

    return index_ATG

    """Returns a list of all the positions of the codon ATG in dna
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a list of integer numbers
    Example: findATG("gatgtatgta") = [1,5] 
    """


def printReadingFrames(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()

    Frame = ["","",""]                        #list of all reading frames 

    for x in range(3):                        #consinders offsets of 0, 1 and 2                      
        for i in range(x,len(dna)-x,3):
            if len(dna[i:i+3]) == 3:
                Frame[x] += dna[i:i+3] + " "
    
    print(Frame[0])
    print(Frame[1])
    print(Frame[2])
    return

    """Prints the sequences of non-overlapping codons in the dna
    with the three possible reading frames in separate columns
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: None
    Example: printReadingFrames("aggcctggc") should print
        agg    ggc    gcc
        cct    ctg    tgg
        ggc
    """

def firstSSR(dna, seq):
    
    if dna == "" or seq == "":
        return None
    
    dna = dna.upper()
    seq = seq.upper()
    count = 0

    if dna.find(seq) == -1:
        return 0

    else:
        x = dna.find(seq)
        while dna[x:x+len(seq)] == seq:   #goes on as long as consiquent nucleotides are the same as the sequence
            count += 1
            x += len(seq)
        return count
    
    """Returns the length (number of repeats) of the first SSR in dna
    that repeats the sequence seq
    Parameters:
        dna: a string object representing a DNA sequence
        seq: a string object representing a short sequence of DNA
    Return value:
        (if seq found): the integer number of seq repeats in the first SSR
        (if seq not found): 0
    Example: firstssr("aggcctggcggcggc", "ggc") = 1
    """

def longestSSR(dna, seq):

    if dna == "" or seq == "":
        return None
    
    dna = dna.upper()
    seq = seq.upper()
    length = 0

    if dna.find(seq) == -1:
        return 0
        
    else:
        while dna.find(seq) != -1:                           
            
            if firstSSR(dna,seq) > length:                   #saves index if length SSR is greater than those preaviously found 
                length = firstSSR(dna,seq)
                
            y = dna.find(seq) + firstSSR(dna,seq)*len(seq)   #index at which SSR ends 
            dna = dna[y:]
            
        return length
    
    """Returns the length of the longest SSR in dna that repeats the sequence seq
    Parameters:
        dna: a string object representing a DNA sequence
        seq: a string object representing a short sequence of DNA
    Return value:
        (if seq found): the integer length of longest SSR in dna repeating seq
        (if seq not found): 0
    Example: longestSSR("aggcctggcggcggc", "gcc") = 3
    """

def longestSSRdin(dna):
    
    if len(dna) < 2:
        return None
    
    dna = dna.upper()
    length = 0
    dinucleotide = []

    for din in dinucleotides:             #considers all dinucleotides (see tuple start)
        SSR_din = longestSSR(dna, din)
        
        if SSR_din > length:
            length = SSR_din
            dinucleotide = [din]

        elif SSR_din == length:
            dinucleotide.append(din)    #saves all SSR of same length

    return (dinucleotide , length)
    
    """Finds the longest SSR in dna for all the possible dinucleotides
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: (if len(dna)>1): a pair (din, len)
        din: a two-letter string object representing the dinucleotide with longest SSR in dna
        len: integer representing the length of the longest SSR of din
                  (if len(dna)<2): None
    Example: longestSSRdin("ctctctgcgccacacaca") = ("ca", 4)
    """

def complement(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    dna1 = ""
    
    for e in dna:
        if e == "G":
            dna1 += "C"
        elif e == "C":
            dna1 += "G"
        elif e == "A":
            dna1 += "T"
        else:
            dna1 += "A"
            
    return dna1

    """Returns the complement of a dna sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing the complement of the DNA sequence
    Example: complement("acgtac") = "tgcatg"
    """

def reverseComplement(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    dna2 = complement(dna) 
    dna3 = dna2[::-1]
    
    return dna3

    """Returns the reverse  complement of a dna sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing the reverse complement of the DNA sequence
    Example: reverseComplement("acgtac") = "gtacgt"
    """
    
def palindrome(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    dna1 = reverseComplement(dna)
    
    if dna1 == dna:
        return True
    else:
        return False
    
    """Returns true if dna is the same as its reverse complement
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: Bool
    Example: palindrome("atat") = True
    """

def dna2rna(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    dna1 = dna.replace('T' , 'U')
    
    return dna1
    
    """Returns a copy of dna in which every "t" has been replaced by a "u"
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing an RNA sequence
    Example: dna2rna("actgat") = "acugau"
    """

def transcribe(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    rna = dna2rna(reverseComplement(dna))
    
    return rna
    
    """Returns the RNA equivalent of the reverse complement of dna
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing an RNA sequence
    Example: transcribe("acgtac") = "guacgu"
    """

def clean(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    c_dna = ""
    
    for e in dna:
        if e != "A" and e != "C" and e != "G" and e != "T":
            c_dna += "N"
        else:
            c_dna += e
            
    return c_dna

    """Returns a new DNA string in which every character in dna
    that is not an "a", "c", "g", or "t" is replaced with an "n"
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing a clean DNA sequence
    Example: clean("goat") = "gnat"
    """

def fix(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    new_dna = ""
    
    for e in dna:
        if e in amb_symb:                       #amb_symb defined at the begining of project
            n = len(amb_symb[e])-1                          #changes element selected
            new_dna += amb_symb[e][random.randint(0, n)]
        else:
            new_dna += e
            
    return new_dna

    """Returns a DNA string in which each ambiguous symbol is replaced
    with one of the possible bases it represents, each with equal probability
    Parameter:
        dna: a string object representing a DNA sequence with ambiguous simbols
    Return value: a string object representing a fixed DNA sequence
    Example: fix("arta") returns either "aata" or "agta"
        (with probability 1/2 for each one of these two possible returns)
    """

def fixAll(dna):
    
    if dna == "":
        return None
    
    dna = dna.upper()
    p = ""
    t = ""
    fixed_dna = [p,]
    fixed_dna2 = []
    
    for e in dna:
        
        if e in amb_symb:                 
        
            for x in range(len(fixed_dna)):
                
                for i in range(len(amb_symb[e])):                                       
                    t = fixed_dna[x] + amb_symb[e][i]  #makes new string; each contains one of the possible solution for the ambiguos symbol
                    fixed_dna2.append(t)
                    t = ""
                    
            fixed_dna, fixed_dna2 = fixed_dna2 , fixed_dna
            fixed_dna2.clear()
                   
        else:                                 #if e is not an ambiguos symbol
            for x in range(len(fixed_dna)):
                fixed_dna[x] += e             #add e to every string present in fixed_dna

    return fixed_dna

    """Returns the list containing all possible DNA strings in which each
    ambiguous symbol in dna is replaced with all of the possible bases it
    represents; all combinations of replacements should be considered and
    returned within the list of fixed DNA strings
    Parameter:
        dna: a string object representing a DNA sequence with ambiguous simbols
    Return value: a list of string objects representing fixed DNA sequences
    Examples:   fixAll("arma") = ["aaaa", "aaca", "agaa", "agca"]
                fixAll("agata") = ["agata"]
    """

def readFASTA(filename):
    
    if filename == "":
        return None 
    
    a  = open(filename , 'r')
    lines = a.readlines()
    a.close()

    return str("".join("".join(lines[1:]).split("\n")))

    """Reads a FASTA file and removes the header and all the newline characters
    Returns the DNA sequence contained in the file as a string
    Parameter:
        filename: the name of file containing a DNA sequence in the FASTA format
    Return value: a string object representing the DNA sequence in the file
    """


def getFASTA(id):
    
    genome = []
    read_file = urllib.request.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+id+'&rettype=fasta&retmode=text')
    
    for line in read_file:
        line = line.decode('utf-8')
        genome.append(line)
        
    return str("".join("".join(genome[1:]).split("\n")))

    """Fetches the DNA sequence with the given id from the NCBI database
    and returns it as a string
    Parameter:
        id: a string object representing the identifier (NCBI accession number) of a DNA sequence
    Return value: a string object containing the dna sequence with the given id
    """

def table(n):           #table of logest prefix that is also suffix for each character in n
    
    n = n.upper()
    if n == "":
        return None
    
    i = 1               #index of character
    j = 0               #index of prefix 
    shift = [0]         #list of best shifts for any character in n
    
    while i<len(n) and i!=len(n):
        
        if n[i] == n[j]:     #finds longest prefix that is also a suffix
            shift.append(j+1)
            i += 1
            j += 1

        elif j > 0:
            j = shift[j-1]    #mismatch - prefix index is recomputed using table
            
        else:
            shift.append(j)   #mismatch with prefix len = 0 
            i += 1
            
    return shift


def findseq(dna,seq):      #with KMP algorithm 

    if dna == "":
        return

    dna = dna.upper()
    seq = seq.upper()
    
    shift = table(seq)      #compute "smart-jumps" 
    i = 0                   #index of char in substring                  
    match = list()         
    
    for j in range(len(dna)):
        
        if i == len(seq):
            match.append(j-i)
            i=0
            
        while i>0 and dna[j]!=seq[i]:
            i=shift[i-1]
            
        if dna[j]==seq[i]:
            i+=1
            
    return match

    """Returns the list of indexes (starting positions) of all the
    occurrences of seq in dna
    Parameters:
        dna: a string object representing a DNA sequence
        seq: a string object representing a sequence of DNA
    Return value: a list of integer numbers
    """


def mark(dna):

    if dna == "":
        return None
    
    dna = dna.upper()
    d_dna = tuple()
    x = 0              #check value
    DNA = ""

    for e in range(0, len(dna)-2, 3):
        d_dna += (dna[e:e+3],)

    
    for e in d_dna:
        if e == "ATG":
            if x == 0:                           #will ATG only if we are not already in coding region
                DNA += start_stop_codons["ATG"]
                x = 1                            #changes value of check so to not change following ATG
            else:
                DNA += e                   
        elif e == "TAA":
            DNA += start_stop_codons["TAA"]
            x = 0
        elif e == "TAG":
            DNA += start_stop_codons["TAG"]
            x = 0
        elif e == "TGA":
            DNA += start_stop_codons["TGA"]
            x = 0
        else:
            DNA += e
            
    return DNA

    """Returns a new DNA string in which every start codon (atg) in dna
    is replaced with ">>>" and every stop codon (taa, tag, or tga) is
    replaced with "<<<"
    The function does not consider overlapping codons (but just the
    reading frame starting from offset 0)
    If two subsequent start codons are found (without a stop codon in
    between), only the first one should be replaced with ">>>", since
    codon atg also codes for the Methionine amino-acid
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object represnting a marked sequence of DNA
    Example: mark("ttgatggagatgcattagaag") = "ttg>>>gagatgcat<<<aag"
    """

def proteins(marked_dna):

    if marked_dna == "":
        return None
    
    codons = tuple(rna2amino.keys())
    gene = list()
    proteins = list()
    p = ""
    x = 0

    while marked_dna.find("<<<",x) != -1:     #makes a list of all coding regions 
        gene.append(marked_dna[marked_dna.find(">>>",x)+3:marked_dna.find("<<<",x)])
        x = marked_dna.find("<<<",x)+3

    for e in range(len(gene)):                  
        gene[e] = gene[e].replace(">>>","ATG")  #fixes mistake if dna is not marked correctly
        gene[e] = dna2rna(gene[e])              #transcript of gene: from dna to rna
        
        for i in range(0, len(gene[e]), 3):   
            cod = gene[e][i:i+3]                #selects all codons
            
            for y in range(len(codons)):               
                if cod == codons[y]: 
                    p += rna2amino[codons[y]]   #appends corresponding aminoacid 

        proteins.append(p)
        p = ""

    return proteins


    """Returns the list of proteins traduced from marked_dna
    Proteins are represented as strings of amino-acids
    Proteins are obtained from the RNA sequences obtained from
    the dna enclosed between the markers ">>>" and "<<<"
    Parameter:
        marked_dna: a string object representing a marked sequence of DNA
    Return Value: a list of string objects representing proteins
    Example: proteins("ttg>>>gagcat<<<aagcag>>>aca>>>caccaacag<<<aga") = ["EH","HQQ"]
    """



""" The following lines define the base settings and functions for plotting through the turtle module
You do not need to add any code to the functions plot and bar!
"""

width = 1200		# width of the window
cols = width // 6	# number of columns of text
height = 600		# height of the window
rows = height // 100	# number of rows of text



def plot(tortoise, index, value, window):
	"""Plots GC fraction value for window ending at position index."""
	
	if (index == window) or (index - window + 1) // cols != (index - window) // cols:
		tortoise.up()	
		tortoise.goto((index - window + 1) % cols, \
		              (index - window + 1) // cols + 0.7 + value * 0.25)
		tortoise.down()
	else:
		tortoise.goto((index - window + 1) % cols, \
		              (index - window + 1) // cols + 0.7 + value * 0.25)

		
def bar(tortoise, index, rf):
	"""Draws a colored bar over codon starting at position index in
	   reading frame rf. Puts the turtle's pen up and down to
	   handle line breaks properly."""
	   
	tortoise.up()
	tortoise.goto(index % cols, index // cols + (rf + 1) / 5)
	tortoise.down()
	tortoise.forward(1)
	tortoise.up()
	tortoise.goto((index + 1) % cols, (index + 1) // cols + (rf + 1) / 5)
	tortoise.down()
	tortoise.forward(1)
	tortoise.up()
	tortoise.goto((index + 2) % cols, (index + 2) // cols + (rf + 1) / 5)
	tortoise.down()
	tortoise.forward(1)


"""
Complete the following two functions to accomplish the required tasks
"""

def orf(dna, rf, tortoise):
    """Finds and draws all ORFs in the reading frame rf
    Blue bars begin and end on start and stop codons
    Parameters:
        dna: a string object representing a sequence of DNA
        rf: reading frame offset (it's value can be either 0, 1, or 2)
        tortoise: the drawing turtle
    Return value: None
    """
	
    # to place a bar in the current color over the codon starting at
    # position index in reading frame rf, call
    # bar(tortoise, index, rf)
    
    if dna == "":
        return None
    
    dna = dna.upper()
    marked_dna = mark(dna[rf:])
    index = rf
    d_dna = ()
    
    for e in range(0, len(dna), 3):
        d_dna += (marked_dna[e:e+3],)

    tortoise.pencolor("red")
    
    for e in d_dna:
        if e == ">>>":
            tortoise.pencolor("blue")
            bar(tortoise, index, rf)
        elif e == "<<<":
            bar(tortoise, index, rf)
            tortoise.pencolor("red")
        else:
            bar(tortoise, index, rf)
        index += 3

    return
	

def gcFreq(dna, window, tortoise):
    """Computes and plots the GC frequency in dna over a sliding window
    Parameters:
        dna: a string object representing a sequence of DNA
        window: integer size of the sliding window
        tortoise: the drawing turtle
    Return value: None
    """
    
    if dna == "":
        return None
    
    # draws red lines at 0.5 above the sequence:
	
    tortoise.pencolor('red')
    for index in range(len(dna) // cols + 1):
    	tortoise.up()
    	tortoise.goto(0, index + 0.825)
    	tortoise.down()
    	if index < len(dna) // cols:
    		tortoise.goto(cols - 1, index + 0.825)
    	else:
    		tortoise.goto((len(dna) - window) % cols, index + 0.825)
    tortoise.up()
    tortoise.pencolor('blue')

    # get initial window count
	
    # get subsequent window counts and plot them
    # to plot a fraction for the window ending at position index,
    # call plot(tortoise, index, fraction, window)

    dna = dna.upper()
    
    for i in range(0, len(dna)-window+1):
        dna_parts = dna[i:i+window]
        value = gcContent(dna_parts)
        plot(tortoise, index, value, window)
        index = i + window
        




"""
The viewer functuion calls the functions you have coded to build the final plot!
"""

def viewer(dna):
    """Displays GC content and ORFs in 3 forward reading frames."""

    dna = dna.upper()   # makes everything upper case

    tortoise = turtle.Turtle()
    screen = tortoise.getscreen()
    screen.setup(width, height)     # makes a long, thin window
    screen.setworldcoordinates(0, 0, cols, rows) # scales coord system so 1 char fits at each point
    screen.tracer(100)
    tortoise.hideturtle()
    tortoise.speed(0)
    tortoise.up()

    # Prints the DNA string in the window:
	
    for index in range(len(dna)):
    	tortoise.goto(index % cols, index // cols)
    	tortoise.write(dna[index], font = ('Courier', 9, 'normal'))
		
    # Finds ORFs in forward reading frames 0, 1, 2:
	
    tortoise.width(5)
    for rf in range(3):
    	orf(dna, rf, tortoise)
		
    # Plots GC frequency:
	
    tortoise.width(1)
    gcFreq(dna, 5, tortoise)

    screen.update()


    
