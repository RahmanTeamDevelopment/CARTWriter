from __future__ import division

from operator import itemgetter
import pysam
import datetime
import os

allowed_chroms = map(str, range(1, 24)) + ['X', 'Y', 'MT']

class TranscriptDBWriter(object):
    """Class for creating new transcript database"""

    def __init__(self, fn, source='', build='', columns=[]):
        """Constructor of the TranscriptDBWriter class"""

        self._fn = fn
        self._source = source
        self._build = build
        self._columns = [x.lower() for x in columns]
        self._records = {c: [] for c in allowed_chroms}
        self.idx_chrom = self._columns.index('chrom')
        self.idx_start = self._columns.index('start')
        self.idx_end = self._columns.index('end')

    def add(self, transcript):
        """Add transcript to DB"""

        record = []
        for c in self._columns:
            if c in ['exons', 'cdna_exons']:
                record.append(','.join([str(e.start) + '-' + str(e.end) for e in getattr(transcript, c.lower())]))
            elif c in ['start', 'end', 'coding_start', 'coding_end', 'cdna_coding_start', 'cdna_coding_end']:
                record.append(int(getattr(transcript, c.lower())))
            else:
                record.append(str(getattr(transcript, c.lower())))
        self._records[transcript.chrom].append(record)

    def _sort_records(self):
        """Sort records by chrom, start, end"""

        idx_start = self._columns.index('start')
        idx_end = self._columns.index('end')
        for c in allowed_chroms:
            if c in self._records:
                self._records[c] = sorted(self._records[c], key=itemgetter(idx_start, idx_end))

    def _index_with_tabix(self):
        """Compress and index output file by Tabix"""

        pysam.tabix_compress(self._fn + '_tmp', self._fn + '.gz', force=True)
        pysam.tabix_index(self._fn + '.gz', seq_col=4, start_col=6, end_col=7, meta_char='#', force=True)


    def finalize(self, options):
        """Write to file, compress and index, clean up"""

        # Sort records by CHROM, START, END
        self._sort_records()

        # Initialize file and write header
        out = open(self._fn + '_tmp', 'w')
        out.write('#createdby: ' + self._source + '\n')
        out.write('#date: ' + str(datetime.datetime.today()).split()[0] + '\n')
        out.write('#build: ' + self._build + '\n')
        out.write('#ncbi_source_db: ' + options.ncbi + '\n')
        out.write('#ucsc_source_db: ' + options.ucsc + '\n')
        out.write('#hgnc_biomart_file: ' + options.hgnc + '\n')

        # Write records to file
        for c in allowed_chroms:
            if c in self._records:
                for record in self._records[c]:
                    record = map(str, record)

                    old_record = [
                        record[0],
                        record[2],
                        record[1],
                        record[3],
                        record[5],
                        '1' if record[4] == '+' else '-1',
                        record[6],
                        record[7],
                        record[-2],
                        str(int(record[9]) + 1),
                        str(int(record[10]) + 1)
                    ]

                    for e in record[8].split(','):
                        [start, end] = e.split('-')
                        old_record.append(start)
                        old_record.append(end)

                    out.write('\t'.join(old_record) + '\n')

        out.close()

        # Compress and index by Tabix
        self._index_with_tabix()

        # Remove temporary file
        os.remove(self._fn + '_tmp')


def output_gff3(transcript, outfile):

    attr = ';'.join(['ID=' + transcript.id, 'HGNCID=' + transcript.hgnc_id, 'GENE_SYMBOL=' + transcript.gene_symbol])
    outfile.write('\t'.join([transcript.chrom, '.', 'transcript', str(transcript.start + 1), str(transcript.end), '.', transcript.strand, '.', attr]) + '\n')

    # Exons
    for i in range(len(transcript.exons)):
        exon = transcript.exons[i]
        exon_id = 'EXON' + transcript.id[4:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + exon_id, 'Parent=' + transcript.id])
        outfile.write('\t'.join([transcript.chrom, '.', 'exon', str(exon.start + 1), str(exon.end), '.', transcript.strand, '.', attr]) + '\n')

    # CDS
    cds_regs = transcript.cds_regions()
    cdspos = 0
    for i in range(len(cds_regs)):
        cds_reg = cds_regs[i]

        cds_id = 'CDS' + transcript.id[4:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + cds_id, 'Parent=' + transcript.id])

        if cdspos % 3 == 0:
            phase = 0
        elif cdspos % 3 == 1:
            phase = 2
        else:
            phase = 1

        outfile.write('\t'.join([transcript.chrom, '.', 'CDS', str(cds_reg[0] + 1), str(cds_reg[1]), '.', transcript.strand, str(phase), attr]) + '\n')
        cdspos += cds_reg[1] - cds_reg[0]


def output_genepred(transcript, outfile):

    exons = transcript.exons
    if transcript.strand == '-':
        exons = exons[::-1]

    record = [
                transcript.id,
                'chr'+transcript.chrom,
                transcript.strand,
                str(transcript.start),
                str(transcript.end),
                str(min(transcript.coding_start, transcript.coding_end)),
                str(max(transcript.coding_start, transcript.coding_end)+1),
                str(len(transcript.exons)),
                ''.join([str(e.start)+',' for e in exons]),
                ''.join([str(e.end)+',' for e in exons]),
                '0',
                'HGNC:'+transcript.hgnc_id,
                'cmpl',
                'cmpl',
                ''.join(frame_offsets(transcript))
            ]

    outfile.write('\t'.join(record) + '\n')


def output_gbk(transcript, ref, outfile):

    full_dna = read_full_dna_sequence(transcript, ref)
    gbk_header(transcript, outfile, full_dna)
    gbk_features(transcript, outfile, ref)
    gbk_origin(outfile, full_dna)
    outfile.write('//\n')


def output_fasta(transcript, outfile, ref):

    width = 70
    s = read_mrna_sequence(transcript, ref)
    outfile.write('>{}\n'.format(transcript.id))
    while len(s)>0:
        outfile.write(s[:width]+'\n')
        s = s[width:]


def get_coding_sequence(transcript, ref):

    ret = ''
    for e in transcript.cds_regions():
        s = ref.read_sequence(transcript.chrom, e[0], e[1])
        if transcript.strand == '-':
            s = reverse_complement(s)
        ret += s
    return ret


def get_protein_sequence(transcript, ref):

    return translate(get_coding_sequence(transcript, ref))[:-1]


def translate(dna):

    gencode = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'}

    ret = ''
    index = 0
    while index + 3 <= len(dna):
        codon = dna[index:index + 3].upper()
        if 'N' in codon:
            ret += '?'
            index += 3
            continue
        ret += gencode[codon]
        index += 3
    return ret


def gbk_header(transcript, outfile, full_dna):

    now = datetime.datetime.now()
    date = '{}-{}-{}'.format(now.strftime("%d"),now.strftime("%b").upper(),now.strftime("%Y"))
    outfile.write('LOCUS       {}              {} bp    DNA . {}\n'.format(transcript.id, len(full_dna), date))
    outfile.write('ACCESSION   {}\n'.format(transcript.id))
    outfile.write('VERSION     {}\n'.format(transcript.id))
    outfile.write('KEYWORDS    .\n')
    outfile.write('SOURCE      Homo sapiens (human)\n')
    outfile.write('  ORGANISM  Homo sapiens\n')
    outfile.write('            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;\n')
    outfile.write('            Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;\n')
    outfile.write('            Catarrhini; Hominidae; Homo.\n')


def gbk_features(transcript, outfile, ref):

    outfile.write('FEATURES             Location/Qualifiers\n')
    feature_mrna(transcript, outfile)
    feature_cds(transcript, outfile, ref)


def feature_mrna(transcript, outfile):

    exon_strings = []
    for e in transcript.exons:
        exon_start = e.start - transcript.start + 1
        exon_end = e.end - transcript.start
        exon_string = '{}..{}'.format(exon_start, exon_end)
        if transcript.strand == '+':
            exon_strings.append(exon_string)
        else:
            exon_strings.append('complement({})'.format(exon_string))

    if len(exon_strings) == 1:
        outfile.write('     mRNA            {}\n'.format(exon_strings[0]))
        return

    first = True
    while len(exon_strings) > 0:
        s = ','.join(exon_strings[:2])
        if first:
            s = '     mRNA            join(' + s
            first = False
        else:
            s = '                     ' + s
        exon_strings = exon_strings[2:]
        if len(exon_strings) > 0:
            s += ','
        else:
            s += ')'
        outfile.write('{}\n'.format(s))


def feature_cds(transcript, outfile, ref):

    exon_strings = []
    for e in transcript.cds_regions():
        exon_start = e[0] - transcript.start + 1
        exon_end = e[1] - transcript.start
        exon_string = '{}..{}'.format(exon_start, exon_end)
        if transcript.strand == '+':
            exon_strings.append(exon_string)
        else:
            exon_strings.append('complement({})'.format(exon_string))

    if len(exon_strings) == 1:
        outfile.write('     CDS             {}\n'.format(exon_strings[0]))
        return

    first = True
    while len(exon_strings) > 0:
        s = ','.join(exon_strings[:2])
        if first:
            s = '     CDS             join(' + s
            first = False
        else:
            s = '                     ' + s
        exon_strings = exon_strings[2:]
        if len(exon_strings) > 0:
            s += ','
        else:
            s += ')'
        outfile.write('{}\n'.format(s))

    comment = '/translation=\"{}\"'.format(get_protein_sequence(transcript, ref))
    while len(comment) > 0:
        line = comment[:58]
        outfile.write('                     {}\n'.format(line))
        comment = comment[58:]


def gbk_origin(outfile, full_dna):

    outfile.write('ORIGIN\n')
    s = full_dna
    lineidx = 0
    while len(s) > 0:
        line = s[:60]
        num = str(lineidx * 60 + 1)
        num = ' ' * (9 - len(num)) + num
        outfile.write('{} {} {} {} {} {} {}\n'.format(num, line[:10], line[10:20], line[20:30], line[30:40], line[40:50], line[50:]))
        lineidx += 1
        s = s[60:]


def read_mrna_sequence(transcript, ref):

    ret = ''

    for e in transcript.exons:
        e_seq = ref.read_sequence(transcript.chrom, e.start, e.end)
        if transcript.strand == '-':
            e_seq = reverse_complement(e_seq)
        ret += e_seq

    return ret


def read_full_dna_sequence(transcript, ref):

    return ref.read_sequence(transcript.chrom, transcript.start, transcript.end)


def frame_offsets(transcript):

    ret = []
    cds_sum = 0
    for exon in transcript.exons:
        cds_region = exon.get_cds(transcript.coding_start, transcript.coding_end)
        if cds_region is not None:
            ret.append(str(cds_sum % 3)+',')
            cds_sum += cds_region[1] - cds_region[0]
        else:
            ret.append('-1,')
    if transcript.strand == '+':
        return ret
    else:
        return ret[::-1]


def reverse_complement(dna_seq):

    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}
    ret = ''
    for base in dna_seq[::-1]:
        ret += complement[base]
    return ret