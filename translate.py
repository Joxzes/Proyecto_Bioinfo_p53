import Bio
import os
import sys
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO 
from Bio.Data import CodonTable

# Leer la secuencia del archivo FASTA
p53 = SeqIO.read('p53.fa', "fasta")
output_file = "p53_prot.fa"

# Obtener la subsecuencia a partir de la posición 142
seq = p53.seq[142:]

# Traducir la subsecuencia hasta el primer codón de parada
translated_seq = seq.translate(to_stop=True)

# Crear un nuevo registro de tipo SeqRecord para la secuencia traducida
translated_record = SeqRecord(
    translated_seq,
    id=p53.id,
    description="Proteína traducida hasta el primer codón de parada"
)

# Guardar el registro traducido en un nuevo archivo FASTA
SeqIO.write(translated_record, output_file, "fasta")

print(f"Secuencia traducida guardada en: {output_file}")
