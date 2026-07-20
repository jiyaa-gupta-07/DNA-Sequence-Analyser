# C THE ATG

**(i) Introduction / Overview:**
‘C the ATG’ is a Python-based bioinformatics toolkit for DNA, mRNA and protein analysis. It simulates the Central Dogma In Molecular Biology
through an intuitive command-line interface with graphical tools.

A single change in a DNA sequence can alter an entire protein and ultimately change a person's life—as seen in disorders such as sickle-cell anaemia. Yet, understanding this relationship has often been abstract. ‘C the ATG’ bridges this gap by allowing users to computationally explore the journey from DNA to mRNA to protein while visualizing each step of the process.

**(ii) List Of Features:**
1. DNA Sequence Validation : Ensures that only biologically valid DNA sequences (A, T, G, and C) are accepted for analysis.
2. Sequence Length Analysis : Reports the total number of nucleotides in the input DNA sequence.
3. Base Composition Analysis : Calculates nucleotide frequencies, purine and pyrimidine composition, GC content, AT content, GC skew, and AT skew, accompanied by graphical         visualizations.
4. Sequence Analysis : Generates the reverse strand and complementary DNA strand, and supports motif searching in both.
5. Molecular Weight Estimation : Calculates the molecular weight of the coding strand, complementary strand, and double-stranded DNA.
6. DNA Transcription : Simulates transcription to generate the corresponding mRNA sequence and supports mRNA motif searching.
7. Protein Translation : Detects valid Open Reading Frames (ORFs) and translates coding sequences using the Standard Genetic Code.
8. Protein Analysis : Reports protein length, amino acid composition, Pareto chart visualization, and protein motif searching.

**(iii) Technologies Used:**
1. Python 3
2. Matplotlib

**(iv) Biological Principles Implemented:**
1. DNA Structure and Complementary Base Pairing 
2. GC and AT Content
3. GC and AT Skew
4. Molecular Weight Estimation
5. Coding and Template DNA Strands
6. Central Dogma In Molecular Biology
7. The Standard Genetic Code
8. Open Reading Frames (ORFs)
9. Start and Stop Codons
10. Point Mutations and their Biological Consequences (Case Study)

**(v) Case Study: Human β-globin (HBB) Gene – Normal vs Sickle Cell Variant**

**(vi) Applications & Future scope:**

**Applications:**
1. Demonstrates the Central Dogma of Molecular Biology.
2. Supports introductory bioinformatics and genetics education.
3. Enables mutation analysis by comparing normal and altered gene sequences.
4. Assists in understanding genetic disorders such as sickle-cell anaemia.
5. Serves as a foundation for basic DNA, mRNA, and protein sequence analysis.
   
**Future Scope:**
1. FASTA & Multi-FASTA support
2. Reverse complement translation
3. Six-frame translation
4. Codon usage analysis
5. DNA restriction map generation
6. PDF report generation
7. GUI/Desktop application
8. Integration with public biological databases such as NCBI

**(vii) What drove us to build ‘C the ATG’?**
The project is  a collaboration between Jiya Gupta, a first-year BS Data Science student with a strong interest in biology, and Sathyaseelan U., a third-year B.Tech student at SRM Institute of Science and Technology with experience in Python.
The project combines Jiya's interest in molecular biology with Sathyaseelan's programming expertise, creating a bioinformatics toolkit built entirely from scratch in Python.

**Why "C the ATG"?**
The name is a play on words—"C" sounds like "See", reflecting our goal of helping users visualize the genetic information encoded in the DNA alphabet (A, T, G, and C). It also references ATG, the DNA triplet that corresponds to the AUG start codon in mRNA, marking the beginning of protein synthesis. Just as AUG initiates translation, C the ATG represents the beginning of our journey into computational biology!
