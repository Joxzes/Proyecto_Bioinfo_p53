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
output_file = "lpato_5.fa"
posicion = 569
mutation = "T"
location = 142 + (posicion - 1)
nucleotido_actual = p53.seq[location]

# print(seq[location])
# print(pato[location])

print(f"Nucleótido actual en la posición {location}: {nucleotido_actual}")

# Realizar el cambio
secuencia_modificada = (
    p53.seq[:location] + mutation + p53.seq[location + 1:]
)

# Crear un nuevo registro FASTA con la secuencia modificada
p53_modificado = SeqRecord(
    Seq(secuencia_modificada),
    id=p53.id + "_modificado",
    description=f"Secuencia con cambio en posición {location} ({nucleotido_actual}->{mutation})"
)

# Guardar la secuencia modificada en un archivo
SeqIO.write(p53_modificado, output_file, "fasta")

print(f"Secuencia modificada guardada en '{output_file}'.")