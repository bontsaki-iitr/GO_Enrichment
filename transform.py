import sys

for line in sys.stdin:
    columns = line.strip('\n').split('\t')
    if len(columns) < 8:
        continue

    chromosome = columns[4]
    strand_val = columns[5]
    gene_name  = columns[6]
    tss_str    = columns[7]

    try:
        tss = int(tss_str)
    except ValueError:
        continue

    col1 = chromosome
    col2 = str(tss)
    col3 = str(tss + 1)
    col4 = f"{col1}@{col2}-{col3}|{gene_name}"
    col5 = "."
    col6 = "+" if strand_val == "1" else ("-" if strand_val == "-1" else ".")

    sys.stdout.write("\t".join([col1, col2, col3, col4, col5, col6]) + "\n")
