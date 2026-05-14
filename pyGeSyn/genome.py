def extract_sequence(fasta_path, chrom, start, end):
    current_chrom = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_chrom is not None and current_chrom == chrom:
                    seq = ''.join(current_seq)
                    if start < 1:
                        start = 1
                    if end > len(seq):
                        end = len(seq)
                    return seq[start - 1:end]
                header = line[1:].split()[0]
                current_chrom = header
                current_seq = []
            elif current_chrom is not None:
                current_seq.append(line)

    if current_chrom is not None and current_chrom == chrom:
        seq = ''.join(current_seq)
        if start < 1:
            start = 1
        if end > len(seq):
            end = len(seq)
        return seq[start - 1:end]

    raise ValueError(f"Chromosome '{chrom}' not found in {fasta_path}")
