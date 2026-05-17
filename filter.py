import sys

desired_chromosomes = {str(i): True for i in range(1, 23)}
desired_chromosomes['X'] = True
desired_chromosomes['Y'] = True

for line in sys.stdin:
    if line.startswith('ensembl_transcript_id'):
        continue
    columns = line.strip('\n').split('\t')
    if len(columns) > 4:
        chromosome_name = columns[4]
        if chromosome_name in desired_chromosomes:
            sys.stdout.write(line)
