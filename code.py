import Bio
import os
import sys
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO 
from Bio.Data import CodonTable


p53 = SeqIO.read('p53.fa', "fasta")
pato = SeqIO.read('pato_1.fa', "fasta").seq
seq = p53.seq
posicion = 569
mutation = "T"
location = 142 + (posicion - 1)
nucleotido_actual = p53.seq[location]

print(f"Nucleótido actual en la posición {location}: {nucleotido_actual}")

secuencia_modificada = (
    p53.seq[:location] + mutation + p53.seq[location + 1:]
)

proteina = secuencia_modificada[142:]
translated_seq = proteina.translate(to_stop=True)

p53_modificado = SeqRecord(
    Seq(secuencia_modificada),
    id=p53.id + "_modificado",
    description=f"Secuencia con cambio en posición {location} ({nucleotido_actual}->{mutation})"
)

translated_record = SeqRecord(
    translated_seq,
    id=p53.id,
    description="Proteína traducida hasta el primer codón de parada"
)

output_file = "lpato_5.fa"
output_file = "p53_prot.fa"
SeqIO.write(translated_record, output_file, "fasta")
SeqIO.write(p53_modificado, output_file, "fasta")

print(f"Secuencia modificada guardada en '{output_file}'.")
print(f"Secuencia traducida guardada en: {output_file}")