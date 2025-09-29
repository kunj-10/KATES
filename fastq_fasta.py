from Bio import SeqIO

fastq_file = "edna_sequences.fastq"
fasta_file = "edna_sequences.fasta"
# Convert
SeqIO.convert(fastq_file, "fastq", fasta_file, "fasta")
print(f"Converted {fastq_file} to {fasta_file}")
