
from tgmi.transcripts import TranscriptDB
import helper
import reference
import gzip


def number_of_input_carts(fn):

    ret = 0
    for line in open(fn):
        line = line.strip()
        if line=='' or line.startswith('#'):
            continue
        ret += 1
    return ret


def read_gene_symbol_file(fn):

    ret = {}
    for line in open(fn):
        line = line.strip()
        if line=='' or line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[1]
    return ret


def read_excluded_list(fn):

    ret = {}
    for line in open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[1]
    return ret


def main(ver, options):

    print 'Input file: {} -> {} CARTs\n'.format(options.input, number_of_input_carts(options.input))

    # Reading gene symbol file
    genes_symbols = read_gene_symbol_file(options.hgnc)
    print 'HGNC BioMart file: {}'.format(options.hgnc)

    # Initializing reference sequence reader
    ref = reference.Reference(options.ref)
    print 'Reference genome file: {}\n'.format(options.ref)

    # Initializing transcript database writer
    columns = ['ID', 'HGNC_ID', 'GENE_SYMBOL', 'INFO', 'STRAND', 'CHROM', 'START', 'END', 'EXONS', 'CODING_START', 'CODING_END', 'CDNA_CODING_START', 'CDNA_CODING_END']
    tdb_writer = helper.TranscriptDBWriter(options.output+'_cava', source='CARTdb ' + ver, build='GRCh37', columns=columns)

    # Reading NCBI database
    db_ncbi = TranscriptDB(options.ncbi)
    db_ncbi.read()
    db_ncbi_excluded = read_excluded_list(options.ncbi[:-3]+'_excluded.txt')
    print 'Transcript database (NCBI mapping): {} -> {} transcripts'.format(options.ncbi,len(db_ncbi._data))

    # Reading UCSC database
    db_ucsc = TranscriptDB(options.ucsc)
    db_ucsc.read()
    db_ucsc_excluded = read_excluded_list(options.ucsc[:-3] + '_excluded.txt')
    print 'Transcript database (UCSC mapping): {} -> {} transcripts\n'.format(options.ucsc, len(db_ucsc._data))

    # Initializing output files
    out_gff = open(options.output + '.gff', 'w')
    out_gff.write('##gff-version 3\n')
    out_source = open(options.output + '_source.txt', 'w')
    out_source.write('#CARTID\trelated_NM\tsource_db\n')
    out_missing = open(options.output + '_missing.txt', 'w')
    out_missing.write('#CARTID\trelated_NM\treason_NCBI_db\treason_UCSC_db\n')
    out_genepred = open(options.output + '.gp', 'w')
    out_fasta = open(options.output + '.fa', 'w')
    if options.gbk:
        out_gbk = gzip.open(options.output + '.gbk.gz', 'w')

    # Iterating through input records
    print 'Processing data ... '
    counter = 0
    counter_ncbi = 0
    counter_ucsc = 0
    counter_missing = 0
    for line in open(options.input):
        line = line.strip()
        if line=='' or line.startswith('#'):
            continue

        cols = line.split()
        cart_id = cols[0]
        nm = cols[1]

        # Create transcript object and write to _source and _missing output files

        if db_ncbi.contains(nm):
            transcript = db_ncbi.by_id(nm)
            out_source.write('{}\t{}\tNCBI\n'.format(cart_id, nm))
            counter_ncbi += 1
        elif db_ucsc.contains(nm):
            transcript = db_ucsc.by_id(nm)
            out_source.write('{}\t{}\tUCSC\n'.format(cart_id, nm))
            counter_ucsc += 1
        else:
            issues = [cart_id, nm]
            if nm in db_ncbi_excluded:
                issues.append(db_ncbi_excluded[nm])
            else:
                issues.append('not_found')
            if nm in db_ucsc_excluded:
                issues.append(db_ucsc_excluded[nm])
            else:
                issues.append('not_found')
            out_missing.write('\t'.join(issues) + '\n')
            counter_missing += 1
            continue

        transcript.id = cart_id

        if transcript.hgnc_id in genes_symbols:
            transcript.gene_symbol = genes_symbols[transcript.hgnc_id]
        else:
            transcript.gene_symbol = '?'
            print '!WARNING: {} ({}) not found in HGNC BioMart file. Gene symbol set to \"?\".'.format(transcript.hgnc_id, cart_id)

        transcript.hgnc_id = transcript.hgnc_id[5:]

        # Adding to transcript database writer
        tdb_writer.add(transcript)

        # Writing to gff file
        helper.output_gff3(transcript, out_gff)

        # Writing to gp file
        helper.output_genepred(transcript, out_genepred)

        # Writing to gbk.gz file
        if options.gbk:
            helper.output_gbk(transcript, ref, out_gbk)

        # Writing to fasta file
        helper.output_fasta(transcript, out_fasta, ref)

        counter += 1

    # Finalizing transcript database writer
    tdb_writer.finalize(options)

    # Closing output files
    out_gff.close()
    out_source.close()
    out_missing.close()
    out_fasta.close()
    out_genepred.close()
    if options.gbk:
        out_gbk.close()

    # Printing out summary info
    print '... done'
    print '\nSummary:'
    print '{} CARTs added to output database ({} with NCBI, {} with UCSC mapping)'.format(counter, counter_ncbi, counter_ucsc)
    print '{} CARTs missing from output database'.format(counter_missing)

    print '\nOutput files:'
    print ' - CAVA database: {}_cava.gz (+.tbi)'.format(options.output)
    print ' - GFF3 file: {}.gff'.format(options.output)
    print ' - GenePred file: {}.gp'.format(options.output)
    if options.gbk:
        print ' - GenBank file: {}.gbk.gz'.format(options.output)
    print ' - FASTA file: {}.fa'.format(options.output)
    print ' - Mapping source file: {}_source.txt'.format(options.output)
    print ' - Missing CARTs file: {}_missing.txt'.format(options.output)