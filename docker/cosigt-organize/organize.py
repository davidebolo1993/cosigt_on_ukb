#!/usr/bin/env python3

import os
import sys
import yaml
import argparse

def validate_graph(graph_path):
    """Validate .og graph file"""
    if not os.path.exists(graph_path):
        print(f'Graph file: {graph_path} does not exist!')
        return False
    if not os.path.isabs(graph_path):
        print(f'Graph file: {graph_path} is not readable!')
        return False
    if not graph_path.endswith('.og'):
        print(f'Graph file: {graph_path} must have .og extension!')
        return False
    return True

def read_graphs_file(graphs_file):
    """Read and validate the graphs specified in the tsv"""
    graphs_out = dict()
    if not os.path.exists(graphs_file):
        print(f'Table with graphs: {graphs_file} does not exist!')
        sys.exit(1)

    with open(graphs_file) as graphs_in:
        for line in graphs_in:
            line = line.rstrip()
            if not line:
                continue
            entries = line.split('\t')
            if len(entries) != 2:
                print(f'Table with graphs: {graphs_file} does not match expected format!')
                sys.exit(1)

            region = entries[0]
            graph = entries[1]

            if region in graphs_out:
                print(f'Table with graphs: {graphs_file} contains duplicate region!')
                sys.exit(1)

            if validate_graph(graph):
                # Keep absolute paths as-is, make relative paths absolute
                graphs_out[region] = graph if os.path.isabs(graph) else os.path.abspath(graph)
            else:
                sys.exit(1)

    print(f'Loaded table with graphs: {graphs_file}!')
    return graphs_out

def validate_alignment(alignment):
    """Validate alignment file"""
    if not os.path.exists(alignment):
        print(f'Alignment file: {alignment} does not exist!')
        return False
    if alignment.endswith('.bam'):
        if not os.path.exists(alignment + '.bai') and not os.path.exists(alignment + '.csi'):
            print(f'Alignment file: {alignment} is not indexed - expected .bai/.csi!')
            return False
    if alignment.endswith('.cram'):
        if not os.path.exists(alignment + '.crai'):
            print(f'Alignment file: {alignment} is not indexed - expected .crai!')
            return False
    return True

def read_alignments_map(alignment_map_file):
    """Read alignments TSV: alignment_path<tab>sample_id"""
    aln_dict = dict()
    if not os.path.exists(alignment_map_file):
        print(f'Alignments table: {alignment_map_file} does not exist!')
        sys.exit(1)

    with open(alignment_map_file, 'r') as map_in:
        for line in map_in:
            line = line.rstrip()
            if not line:
                continue
            try:
                aln_path, aln_id = line.split('\t')
            except ValueError:
                print(f'Alignments table has invalid line: {line}')
                sys.exit(1)

            # Keep absolute paths as-is, make relative paths absolute
            aln_path = aln_path if os.path.isabs(aln_path) else os.path.abspath(aln_path)

            if aln_id in aln_dict:
                print(f'Duplicate sample id: {aln_id}!')
                sys.exit(1)

            if not validate_alignment(aln_path):
                sys.exit(1)

            aln_dict[aln_id] = aln_path

    print(f'Loaded alignments from: {alignment_map_file}!')
    return aln_dict

def read_bed(bed_file, graphs_dict):
    """Read bed file and organize regions"""
    bed_dict = dict()
    if not os.path.exists(bed_file):
        print(f'Bed file: {bed_file} does not exist!')
        sys.exit(1)

    with open(bed_file, 'r') as bedin:
        for line in bedin:
            bed_entry = line.rstrip().split('\t')
            if len(bed_entry) < 3:
                print(f'Invalid entry in {bed_file}: less than 3 columns')
                sys.exit(1)

            chrom = bed_entry[0]
            start = bed_entry[1]
            end = bed_entry[2]
            annot = bed_entry[3] if len(bed_entry) >= 4 else "unknown"
            alt = bed_entry[4] if len(bed_entry) == 5 else None

            region = '_'.join([chrom, start, end])

            if region not in graphs_dict:
                print(f'Region: {region} in bed file is missing in the graphs table!')
                sys.exit(1)

            if region not in bed_dict:
                bed_dict[region] = [(chrom, start, end, annot, alt)]
            else:
                bed_dict[region].append((chrom, start, end, annot, alt))

    print(f'Loaded bed file: {bed_file}!')
    return bed_dict

def read_genome(genome_file):
    """Read and validate reference genome"""
    if not os.path.exists(genome_file):
        print(f'Genome file: {genome_file} does not exist!')
        sys.exit(1)
    if not genome_file.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz', '.fna', '.fna.gz')):
        print(f'Genome file: {genome_file} must be FASTA format!')
        sys.exit(1)
    if not os.path.exists(genome_file + '.fai'):
        print(f'Genome file: {genome_file} is not indexed!')
        sys.exit(1)

    print(f'Loaded reference: {genome_file}!')
    # Keep absolute paths as-is, make relative paths absolute
    return genome_file if os.path.isabs(genome_file) else os.path.abspath(genome_file)

def write_graphs(graphs_dict, bed_dict, RESOURCES):
    """Write graph symlinks to resources/assemblies"""
    graphs_dir = os.path.join(RESOURCES, 'assemblies')
    os.makedirs(graphs_dir, exist_ok=True)

    for region, graph_path in graphs_dict.items():
        if region in bed_dict:
            chrom = region.split('_')[0]
            chrom_folder = os.path.join(graphs_dir, chrom)
            os.makedirs(chrom_folder, exist_ok=True)
            graph_link = os.path.join(chrom_folder, f'{region}.og')
            if os.path.lexists(graph_link):
                os.remove(graph_link)
            os.symlink(graph_path, graph_link)

    print(f'Added graphs to: {graphs_dir}!')

def write_alignments(aln_dict, RESOURCES):
    """Write alignment symlinks to resources/alignments"""
    aln_dir = os.path.join(RESOURCES, 'alignments')
    os.makedirs(aln_dir, exist_ok=True)

    samples = []
    for sample_id, aln_path in aln_dict.items():
        if aln_path.endswith('.bam'):
            aln_link = f'{sample_id}.bam'
            idx_ext = '.bai' if os.path.exists(aln_path + '.bai') else '.csi'
        else:
            aln_link = f'{sample_id}.cram'
            idx_ext = '.crai'

        link_path = os.path.join(aln_dir, aln_link)
        if os.path.lexists(link_path):
            os.remove(link_path)
        os.symlink(aln_path, link_path)

        idx_link = os.path.join(aln_dir, aln_link + idx_ext)
        if os.path.lexists(idx_link):
            os.remove(idx_link)
        os.symlink(aln_path + idx_ext, idx_link)

        samples.append(sample_id)

    print(f'Added alignments to: {aln_dir}!')
    return samples

def write_regions(bed_dict, RESOURCES):
    """Write region bed files to resources/regions"""
    reg_dir = os.path.join(RESOURCES, 'regions')
    os.makedirs(reg_dir, exist_ok=True)

    regions = []
    for region, bed_entries in bed_dict.items():
        bed_dir = os.path.join(reg_dir, region.split('_')[0])
        os.makedirs(bed_dir, exist_ok=True)

        for subr in bed_entries:
            region_name = '_'.join(subr[:-2])
            bed_out = os.path.join(bed_dir, f'{region_name}.bed')
            with open(bed_out, 'w') as b_out:
                b_out.write('\t'.join(subr[:-1]) + '\n')

                if subr[-1] is not None:
                    alts = subr[-1].split(',')
                    for alt in alts:
                        last = alt.rfind(':')
                        chr_alt, rest_alt = alt[:last], alt[last+1:]
                        start_alt, end_alt = rest_alt.split('-')
                        b_out.write('\t'.join([chr_alt, start_alt, end_alt]) + '\n')

            regions.append(region_name)

    print(f'Added regions to: {reg_dir}!')
    return regions

def write_reference(reference_path, RESOURCES):
    """Write reference symlink to resources/reference"""
    ref_dir = os.path.join(RESOURCES, 'reference')
    os.makedirs(ref_dir, exist_ok=True)

    ref_name = os.path.basename(reference_path)
    ref_link = os.path.join(ref_dir, ref_name)
    if os.path.lexists(ref_link):
        os.remove(ref_link)
    os.symlink(reference_path, ref_link)

    fai_link = os.path.join(ref_dir, ref_name + '.fai')
    if os.path.lexists(fai_link):
        os.remove(fai_link)
    os.symlink(reference_path + '.fai', fai_link)

    # Return RELATIVE path for config
    ref_rel = os.path.join('resources', 'reference', ref_name)
    print(f'Added reference to: {ref_dir}!')
    return ref_rel

def write_config(samples, regions, reference, output_dir, pansn_prefix, CONFIG):
    """Write config.yaml"""
    config = {
        'samples': samples,
        'regions': regions,
        'reference': reference,
        'output': os.path.abspath(output_dir),
        'pansn_prefix': pansn_prefix
    }

    config_out = os.path.join(CONFIG, 'config.yaml')
    with open(config_out, 'w') as yml_out:
        yaml.dump(config, yml_out, default_flow_style=False)

    print(f'Wrote config to: {config_out}!')

def setup_arg_parser():
    """Argument parser"""
    parser = argparse.ArgumentParser(
        prog='organize.py',
        description='Organize inputs for COSIGT pipeline'
    )

    parser.add_argument('-g', '--graphs', required=True, help='Graphs TSV')
    parser.add_argument('-r', '--reads', required=True, help='Alignments TSV')
    parser.add_argument('-b', '--bed', required=True, help='Regions BED')
    parser.add_argument('-f', '--reference', required=True, help='Reference FASTA')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--pansn', default='grch38#1#', help='PanSN prefix')

    return parser

def main():
    """Main function"""
    parser = setup_arg_parser()
    args = parser.parse_args()

    BASE = os.getcwd()
    RESOURCES = os.path.join(BASE, 'resources')
    CONFIG = os.path.join(BASE, 'config')

    if os.path.isdir(RESOURCES) and os.listdir(RESOURCES):
        print(f"Directory {RESOURCES} is not empty: clean it and retry!")
        sys.exit(1)

    os.makedirs(RESOURCES, exist_ok=True)
    os.makedirs(CONFIG, exist_ok=True)

    print("=" * 60)
    print("Reading inputs...")
    print("=" * 60)

    graphs_dict = read_graphs_file(args.graphs)
    bed_dict = read_bed(args.bed, graphs_dict)
    aln_dict = read_alignments_map(args.reads)
    reference = read_genome(args.reference)

    print("=" * 60)
    print("Writing resources...")
    print("=" * 60)

    write_graphs(graphs_dict, bed_dict, RESOURCES)
    samples = write_alignments(aln_dict, RESOURCES)
    regions = write_regions(bed_dict, RESOURCES)
    reference_link = write_reference(reference, RESOURCES)

    print("=" * 60)
    print("Writing configuration...")
    print("=" * 60)

    write_config(samples, regions, reference_link, args.output, args.pansn, CONFIG)

    print("=" * 60)
    print("✓ Setup complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()

