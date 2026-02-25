#!/usr/bin/env python
import argparse as arg
import subprocess as sp
import os as os


parser = arg.ArgumentParser()
parser.add_argument('-i', '--input', required=True)
parser.add_argument('-o', '--output', required=True) 
args = parser.parse_args()


print('FindCleave v0.0.1')
print('Kevin Kuchinski - https://github.com/KevinKuchinski/')


HPAI_H5_CDS = 'MEKIVLLLAIVSLVKSDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKTHNGKLCDLNGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIVEKASPANDLCYPGDFNDYEE'
HPAI_H5_CDS += 'LKHLLSRTNHFEKIQIIPKSSWSNHDASSGVSSACPYHGRSSFFRNVVWLIKKNSAYPTIKRSYNNTNQEDLLVLWGIHHPNDAAEQTKLYQNPTTYISVGTSTLNQRLVPEIAT'
HPAI_H5_CDS += 'RPKVNGQSGRMEFFWTILKPNDAINFESNGNFIAPEYAYKIVKKGDSAIMKSELEYGNCNTKCQTPMGAINSSMPFHNIHPLTIGECPKYVKSNRLVLATGLRNTPQRERRRKKRGLF'
HPAI_H5_CDS += 'GAIAGFIEGGWQGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGVTNKVNSIIDKMNTQFEAVGREFNNLERRIENLNKQMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKV'
HPAI_H5_CDS += 'RLQLRDNAKELGNGCFEFYHKCDNECMESVKNGTYDYPQYSEEARLNREEISGVKLESMGTYQILSIYSTVASSLALAIMVAGLSLWMCSNGSLQCRICI'
ref_cleave_start = 337
ref_cleave_end = 349
min_cov = 90


print('Aligning sequences to prototypical HPAI H5...')
input_dir, input_file = os.path.split(args.input)
db_path = os.path.join(input_dir, input_file + '_blast_db.fa')
with open(db_path, 'w') as db_file:
    db_file.write('>HPAI_H5\n')
    db_file.write(HPAI_H5_CDS + '\n')
terminal_command = f'makeblastdb -in {db_path} -dbtype prot'
sp.run(terminal_command, stdout=sp.DEVNULL, stderr=sp.DEVNULL, shell=True)
input_path = args.input
blast_path = os.path.join(input_dir, input_file + '_blast_results.tsv')
cols = 'qseqid bitscore sstart send qstart qend qseq'
terminal_command = f'blastx -query {input_path} -db {db_path} -outfmt "6 {cols}" > {blast_path}'
sp.run(terminal_command, stdout=sp.DEVNULL, stderr=sp.DEVNULL, shell=True)
for suffix in ['', '.phr', '.pin', '.psq', '.pdb', '.pot', '.ptf', '.pto']:
    os.remove(db_path + suffix)


print('Identifying cleavage sites...')
best_bitscores = {}
with open(blast_path, 'r') as blast_file:
    for line in blast_file:
        qseqid, bitscore, sstart, send, qstart, qend, qseq = line.strip().split('\t')
        best_bitscores[qseqid] = 0
with open(blast_path, 'r') as blast_file:
    for line in blast_file:
        qseqid, bitscore, sstart, send, qstart, qend, qseq = line.strip().split('\t')
        bitscore = float(bitscore)
        best_bitscores[qseqid] = bitscore if bitscore > best_bitscores[qseqid] else best_bitscores[qseqid]
cleave_seqs = {}
with open(blast_path, 'r') as blast_file:
    for line in blast_file:
        qseqid, bitscore, sstart, send, qstart, qend, qseq = line.strip().split('\t')
        if float(bitscore) == best_bitscores[qseqid] and int(sstart) <= ref_cleave_start and int(send) >= ref_cleave_end:
            best_bitscores[qseqid] = float(bitscore)
            start = ref_cleave_start - int(sstart)
            end = int(send) - ref_cleave_end
            cleave_seq = qseq[start:-end].replace('-', '')
            cleave_start = int(qstart) + (start * 3)
            cleave_end = int(qend) - (end * 3)
            cleave_seqs[qseqid] = (bitscore, cleave_start, cleave_end, cleave_seq)
os.remove(blast_path)


with open(args.output, 'w') as output_file:
    output_file.write('seq_name\tbitscore\tcleave_start\tcleave_end\tcleavage_site_motif\n')
    for seq_name, cleave_info in cleave_seqs.items():
        line = [seq_name] + [str(s) for s in cleave_info]
        output_file.write('\t'.join(line) + '\n')
print('Done.')