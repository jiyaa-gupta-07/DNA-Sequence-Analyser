import matplotlib.pyplot as plt
import numpy as np

print("=" * 70)
print("                           C the ATG")
print("                  Analyzing DNA sequences through Python.")
print("=" * 70)
print()
def is_valid_dna(sequence):
    if sequence == "":
        return False, "ERROR: DNA sequence can not be empty."

    if not all(base in "ATGC" for base in sequence):
        return False, "ERROR: DNA sequence can only contain the letters A,T,G and C."

    return True, "Your DNA sequence is valid!"

def is_valid_dna_motif(sequence):
    if sequence == "":
        return False, "DNA motif search skipped."

    elif not all(base in "ATGC" for base in sequence):
        return False, "ERROR: DNA sequence can only contain the letters A,T,G and C."

    return True, "DNA motif validated successfully!"


def is_valid_dna_complement_motif(sequence):
    if sequence == "":
        return False, "Complementary DNA motif search skipped."

    elif not all(base in "ATGC" for base in sequence):
        return False, "ERROR: DNA sequence can only contain the letters A,T,G and C."

    return True, "Complementary DNA motif validated successfully!"

def is_valid_mrna_motif(sequence):
    if sequence == "":
        return False, "mRNA motif search skipped."

    elif not all(base in "AUGC" for base in sequence):
        return False, "ERROR: mRNA sequence can only contain the letters A,U,G and C."

    return True, "mRNA motif validated successfully!"


def length_of_dna(sequence):
    return len(sequence)

def complement_of_dna(sequence):
   if sequence =="":
      return "ERROR: DNA sequence is empty "
   complement_sequence = ""
   for i in sequence:
      if i == "A":
         complement_sequence+="T"
      elif i =="T":
         complement_sequence+="A"
      elif i =="G":
         complement_sequence+="C"
      elif i =="C":
         complement_sequence+="G"
      else:
         print("")
   return complement_sequence

def transcription(sequence):
   if sequence =="":
      return "ERROR: DNA Sequence is Empty"
   mRNA_sequence = ""
   for i in sequence:
      if(i=="T"):
         mRNA_sequence+="U"
      else:
         mRNA_sequence+=i
   return mRNA_sequence


def translation(sequence,codon_table):
   protein = []
   for i in range(0,len(sequence)-2,3):
      codon = sequence[i:i+3]
      amino_acid = codon_table.get(codon,"Unknown")
      if (amino_acid == "STOP"):
         break
      protein.append(amino_acid)
   return "-".join(protein)
   
dna = input("Your DNA sequence: ").strip().upper()

valid, message = is_valid_dna(dna)
print(message)
if valid:
 A = dna.count("A")
 T = dna.count("T")
 G = dna.count("G")
 C = dna.count("C")
 PURINES= A+G
 PYRIMIDINES= T+C
 GC_CONTENT = (((G+C)/length_of_dna(dna))*100)
 AT_CONTENT = (((A+T)/length_of_dna(dna))*100)

 if (G + C) != 0:
    GC_SKEW = (G - C) / (G + C)
 else:
    GC_SKEW = None

 if (A + T) != 0:
    AT_SKEW = (A - T) / (A + T)
 else:
    AT_SKEW = None

 print()
 print("------------------------------SEQUENCE LENGTH------------------------------")
 print()
 print("SEQUENCE LENGTH:", length_of_dna(dna), "base pairs")
 print()
 print("----------------------------COMPOSITION OF BASES----------------------------")
 print()
 print("NUCLEOTIDE FREQUENCY:")
 print("A                         :", A)
 print("T                         :", T)
 print("G                         :", G)
 print("C                         :", C)
 plt.bar(["A","T","G","C"],[A,T,G,C],color="grey")
 plt.xlabel("Nucleotide Bases")
 plt.ylabel("Frequency (Count)")
 plt.title("Nucleotide Frequency Distribution",weight="bold")
 plt.show()

 print("PURINES AND PYRIMIDINES: ")
 print("PURINES (A+G)             :", PURINES)
 print("PYRIMIDINES (T+C)         :", PYRIMIDINES)
 print("Percentage of Purines     :",round((PURINES/len(dna)*100),3),"%")
 print("Percentage of Pyrimidines :",round((PYRIMIDINES/len(dna)*100),3),"%")
 print("GC Content                :",round(GC_CONTENT,3))
 print("AT Content                :",round(AT_CONTENT,3))
 
 if GC_SKEW is not None:
    print("GC Skew:", round(GC_SKEW, 3))
 else:
    print("GC Skew could not be calculated since G + C = 0.")

 if AT_SKEW is not None:
    print("AT Skew:", round(AT_SKEW, 3))
 else:
    print("AT Skew could not be calculated since A + T = 0.")

 plt.figure(figsize=(6,6))
 plt.pie(x=[PURINES,PYRIMIDINES],labels=["Purines","Pyrimidines"],colors=["lightgrey","dimgrey"],startangle=90,autopct="%1.1f%%",wedgeprops={"edgecolor": "white","linewidth": 2},textprops={"fontsize":12, "weight":"bold"})
 plt.title("Purines vs Pyrimidines Composition",weight="bold")
 plt.show()
 
 print()
 print("------------------------------SEQUENCE ANALYSIS------------------------------")
 print()
 reverse = dna[::-1]
 print("The reverse of the input DNA strand is: ",reverse)
 print()
 print("Considering the input DNA strand as the coding strand (or) sense strand (or) non-template strand:- ")
 print()
 complement = complement_of_dna(dna)
 print("The complementary non-coding (or) anti-sense (or) template strand is: ",complement)
 search_in_dna= input("Enter the DNA motif to search (e.g., GAATTC) in the input DNA strand (Press Enter to skip): ").strip().upper()
 valid, message = is_valid_dna_motif(search_in_dna)
 print(message)

 if valid:
    positions = []

    for i in range(len(dna) - len(search_in_dna) + 1):
        current_sequence = dna[i:i+len(search_in_dna)]

        if current_sequence == search_in_dna:
            positions.append(i + 1)

    if len(positions) > 0:
        all_positions = ", ".join(map(str, positions))
        print(f"DNA motif '{search_in_dna}' found.")
        print(f"Occurrences : {len(positions)}")
        print(f"Position(s) : {all_positions}")

    else:
        print(f"DNA motif '{search_in_dna}' was not found.")
 print()
 search_in_dna_complement= input("Enter the DNA motif to search (e.g., GAATTC) in the complementary DNA strand (Press Enter to skip): ").strip().upper()
 valid, message = is_valid_dna_complement_motif(search_in_dna_complement)
 print(message)
 if valid:
    positions = []

    for i in range(len(complement) - len(search_in_dna_complement) + 1):
        current_sequence = complement[i:i + len(search_in_dna_complement)]

        if current_sequence == search_in_dna_complement:
            positions.append(i + 1)

    if len(positions) > 0:
        all_positions = ", ".join(map(str, positions))
        print(f"DNA motif '{search_in_dna_complement}' found.")
        print(f"Occurrences : {len(positions)}")
        print(f"Position(s) : {all_positions}")

    else:
        print(f"DNA motif '{search_in_dna_complement}' was not found.")
 print()
 print("------------------------------MOLECULAR WEIGHT------------------------------")
 print()
 MWSSC= 313.21*A + 304.20*T + 329.21*G + 289.18*C
 MWDS= 617.4 * (A+T) + 618.4*(G+C) + 36
 Ac = complement.count("A")
 Tc = complement.count("T")
 Gc = complement.count("G")
 Cc = complement.count("C")
 MWSSNC= 313.21*Ac + 304.20*Tc + 329.21*Gc + 289.18*Cc
 print("Molecular Weight of the Entered Input (CODING) Single Strand             :", round(MWSSC, 3),"Da") 
 print("Molecular Weight of the Complementary (NON CODING) Single Strand         :", round(MWSSNC, 3),"Da") 
 print("Molecular Weight of the Double Stranded linear DNA (NON CODING + CODING) :", round(MWDS, 3),"Da") 
 print()
 print("------------------------------TRANSCRIPTION: FROM DNA TO mRNA!------------------------------")

 mRNA = transcription(dna)
 print("Transcribed mRNA: ",mRNA)
 print()
 search_in_mrna= input("Enter the mRNA motif to search (e.g., AUG, UAA) in the mRNA strand (Press Enter to skip): ").strip().upper()
 valid, message = is_valid_mrna_motif(search_in_mrna)
 print(message)
 if valid:
    positions = []

    for i in range(len(mRNA) - len(search_in_mrna) + 1):
        current_sequence = mRNA[i:i + len(search_in_mrna)]

        if current_sequence == search_in_mrna:
            positions.append(i + 1)

    if len(positions) > 0:
        all_positions = ", ".join(map(str, positions))
        print(f"mRNA motif '{search_in_mrna}' found.")
        print(f"Occurrences : {len(positions)}")
        print(f"Position(s) : {all_positions}")

    else:
        print(f"mRNA motif '{search_in_mrna}' was not found.")
 print()
 print("------------------------------TRANSLATION: FROM mRNA TO PROTEIN!------------------------------")
 print()
 def cleaned_mRNA(mRNA):

    start_index = mRNA.find("AUG")

    if start_index == -1:
        return "", -1, -1, None

    stop_codons = ["UAA", "UAG", "UGA"]

    for i in range(start_index, len(mRNA), 3):

        codon = mRNA[i:i+3]

        if len(codon) < 3:
            break

        if codon in stop_codons:
            clean_mRNA = mRNA[start_index:i]
            end_index = i
            end_codon = codon

            return clean_mRNA, start_index, end_index, end_codon

    return "", start_index, -1, None
 

 clean_mRNA, start_index, end_index, end_codon = cleaned_mRNA(mRNA)

 print("Open Reading Frame (ORF) Detection")

 if start_index == -1:
     print("Start Codon                         : Not Found")
     print("Stop Codon                          : Not Applicable")
     print("ORF Status                          : No Open Reading Frame Detected")
     print()
     print("Translation                         : Aborted")
     print("Reason                              : No start codon (AUG) detected.")

 elif end_index == -1:
     print(f"Start Codon                         : Found (AUG) at Position {start_index + 1}")
     print("Stop Codon                           : Not Found")
     print("ORF Status                           : Incomplete ORF")
     print()
     print("Translation                          : Aborted")
     print("Reason                               : No in-frame stop codon detected.")

 else:
     orf_length = end_index - start_index + 3

     print(f"Start Codon                           : Found (AUG) at Position {start_index + 1}")
     print(f"Stop Codon                            : Found ({end_codon}) at Position {end_index + 1}")
     print(f"ORF Length (including stop)           : {orf_length} Nucleotides")
     print("ORF Status                             : Complete ORF Detected")

 Codon_table = {
   'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
   'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
   'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
   'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
   'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
   'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
   'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
   'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
   'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
   'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
   'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
   'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
   'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
   'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
   'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
   'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
   }

 amino_acids = {
    "Ala": "Alanine",
    "Arg": "Arginine",
    "Asn": "Asparagine",
    "Asp": "Aspartic Acid",
    "Cys": "Cysteine",
    "Gln": "Glutamine",
    "Glu": "Glutamic Acid",
    "Gly": "Glycine",
    "His": "Histidine",
    "Ile": "Isoleucine",
    "Leu": "Leucine",
    "Lys": "Lysine",
    "Met": "Methionine",
    "Phe": "Phenylalanine",
    "Pro": "Proline",
    "Ser": "Serine",
    "Thr": "Threonine",
    "Trp": "Tryptophan",
    "Tyr": "Tyrosine",
    "Val": "Valine"
 }
 
 protein_sequence = translation(clean_mRNA, Codon_table)
 if protein_sequence == "":
     print("No protein sequence available for analysis.")
 else:
     print(f"Translated Amino Acid Sequence      : {protein_sequence}")
     print()
     print("------------------------------PROTEIN ANALYSIS------------------------------")
     print()
     protein_list = protein_sequence.split("-")

     protein_length = len(protein_list)

     frequency = {}
 
     for amino in amino_acids:
      frequency[amino] = 0

     for amino in protein_list:
       if amino in frequency:
        frequency[amino] += 1
       else:
        frequency[amino] = 1

     sorted_frequency = dict(sorted(frequency.items(), key=lambda item: item[1], reverse=True))

     max_frequency = max(frequency.values())
     most_abundant = []

     for amino, count in frequency.items():
      if count == max_frequency:
        most_abundant.append(amino)

     non_zero_frequency = {amino: count for amino, count in frequency.items() if count > 0}

     min_frequency = min(non_zero_frequency.values())
     least_abundant = []
 
     for amino, count in non_zero_frequency.items():
      if count == min_frequency:
        least_abundant.append(amino)

     print(f"Protein Length : {protein_length} amino acids")
     print("Most Abundant Amino acid(s)            :", end=" ")
     for amino in most_abundant:
      print(f"{amino_acids[amino]} ({amino})", end=", ")

     print()

     print("Least Abundant Amino acid(s)           :", end=" ")
     for amino in least_abundant:
      print(f"{amino_acids[amino]} ({amino})", end=", ")
 
     print()

     print("Absent Amino acid(s)                   :", end=" ")
     for amino in frequency:
      if(frequency[amino]==0):
       print(f"{amino_acids[amino]} ({amino})", end=", ")
 
     print()
     print("-" * 45)
     print(f"{'AMINO ACID':<30}{'FREQUENCY':^10}")
     print("-" * 45)
     for code, frequency_count in sorted_frequency.items():
      amino_display = f"{amino_acids[code]} ({code})"
      print(f"{amino_display:<30}{frequency_count:^10}")
      print("-" * 45)

     pareto_data = {aa: count for aa, count in sorted_frequency.items() if count > 0}
     if (len(pareto_data) > 0):
         labels = list(pareto_data.keys())
         counts = list(pareto_data.values())

         cumulative = np.cumsum(counts)
         cumulative_percent = cumulative / cumulative[-1] * 100
 
         fig, ax1 = plt.subplots(figsize=(12,6))

         ax1.bar(labels, counts, color="grey")
         ax1.set_xlabel("Amino Acids", weight="bold")
         ax1.set_ylabel("Frequency", weight="bold")
         ax1.set_title("Amino Acid vs  Frequencies", weight="bold")

         ax2 = ax1.twinx()
         ax2.plot(labels, cumulative_percent, color="red", marker="o", linewidth=2)
         ax2.set_ylabel("Cumulative Percentage (%)", weight="bold")
         ax2.set_ylim(0, 110)
         ax2.axhline(80, color="blue", linestyle="--", linewidth=1.5, label="80%")
         plt.xticks(rotation=45, ha="right")
         plt.tight_layout()
         plt.show()

     else:
         print("No amino acids available to plot.")


