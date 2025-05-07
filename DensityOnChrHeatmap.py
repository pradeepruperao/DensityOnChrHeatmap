#!/usr/bin/env python3
# Usage: python3 scriptname.py 
# Input files are given in the last line of this script

import svgwrite
import statistics

def read_data(filename):
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            chrom, start, end, value = line.strip().split()
            if chrom not in data:
                data[chrom] = []
            data[chrom].append((int(start), int(end), float(value)))
    return data

def get_chromosome_length(data):
    # Compute chromosome length up to the last gene's end position only
    return max(end for _, end, _ in data)

def create_chromosome_shape(x, y, width, height):
    radius = height / 2
    path = svgwrite.path.Path(d=f"M {x+radius},{y}")
    path.push(f"L {x+width-radius},{y}")
    path.push(f"A {radius},{radius} 0 0 1 {x+width-radius},{y+height}")
    path.push(f"L {x+radius},{y+height}")
    path.push(f"A {radius},{radius} 0 0 1 {x+radius},{y}")
    path.push("Z")
    return path

def create_visualization(genes_file, snps_file, output_file):
    genes_data = read_data(genes_file)
    snps_data = read_data(snps_file)

    chromosomes = sorted(genes_data.keys(), key=lambda x: int(x[3:]))
    svg_width = 1200
    chr_height = 14
    chr_spacing = 120
    margin_top = 160
    margin_bottom = 100
    margin_side = 100
    svg_height = len(chromosomes) * (chr_height + chr_spacing) + margin_top + margin_bottom
    #print("Chrom === ",len(chromosomes))

    dwg = svgwrite.Drawing(output_file, size=(svg_width, svg_height))
    dwg.add(dwg.rect(insert=(0, 0), size=('100%', '100%'), fill='white'))

    chr_lengths = {}
    for chrom in chromosomes:
        gene_ends = [end for _, end, _ in genes_data[chrom]]
        chr_lengths[chrom] = max(gene_ends)
        #print("chrlen == ",chrom," : ",genes_data[chrom])
    max_chr_len = max(chr_lengths.values())
    gene_values = [value for chrom in genes_data for _, _, value in genes_data[chrom]]
    max_snp_value = max(value for chrom in snps_data for _, _, value in snps_data[chrom])

    mean = statistics.mean(gene_values)
    stdev = statistics.stdev(gene_values)
    min_clamp = mean - stdev
    max_clamp = mean + stdev

    colors = ['#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061']

    for i, chrom in enumerate(chromosomes):
        y = margin_top + i * (chr_height + chr_spacing)
        chr_width = (svg_width - 2 * margin_side) * (chr_lengths[chrom] / max_chr_len)

        chr_shape = create_chromosome_shape(margin_side, y, chr_width, chr_height)
        chr_shape.fill('none').stroke('black', width=1)
        dwg.add(chr_shape)

        clip_path = dwg.defs.add(dwg.clipPath(id=f'clip-{chrom}'))
        clip_path.add(chr_shape)

        heatmap_group = dwg.g(clip_path=f'url(#clip-{chrom})')
        for start, end, value in genes_data[chrom]:
            x = margin_side + ((start - 1) / chr_lengths[chrom]) * chr_width
            width = max(1, ((end - start + 1) / chr_lengths[chrom]) * chr_width)
            clamped_value = max(min(value, max_clamp), min_clamp)
            color_index = min(int((clamped_value - min_clamp) / (max_clamp - min_clamp) * len(colors)), len(colors) - 1)
            heatmap_group.add(dwg.rect(insert=(x, y), size=(width, chr_height),
                                       fill=colors[color_index], stroke='none'))
        dwg.add(heatmap_group)

        points = [(margin_side + ((start - 1) / chr_lengths[chrom]) * chr_width,
                   y - 30 - (value / max_snp_value) * 100)
                  for start, end, value in snps_data[chrom]]
        snp_line = dwg.polyline(points=points, fill='none', stroke='red', stroke_width=1.5)
        dwg.add(snp_line)

        y_axis_x = margin_side - 10
        y_top = min(pt[1] for pt in points)
        y_bottom = y - 30
        dwg.add(dwg.line(start=(y_axis_x, y_top), end=(y_axis_x, y_bottom), stroke='black', stroke_width=1))

        dwg.add(dwg.line(start=(margin_side, y_bottom), end=(margin_side + chr_width, y_bottom), stroke='black', stroke_width=1))

        dwg.add(dwg.text(chrom, insert=(margin_side - 20, y + chr_height/2 + 5),
                         font_size=20, fill='black', text_anchor='end'))

    scale_y = svg_height - margin_bottom / 2
    for i in range(6):
        x = margin_side + i * (svg_width - 2 * margin_side) / 5
        dwg.add(dwg.line(start=(x, scale_y), end=(x, scale_y + 5), stroke='black', stroke_width=1))
        dwg.add(dwg.text(f'{i * max_chr_len / 5 / 1e6:.1f} Mb', insert=(x, scale_y + 20),
                         font_size=20, text_anchor='middle'))

    heatmap_scale_width = 100
    heatmap_scale_height = 20
    for i, color in enumerate(colors):
        x = svg_width - margin_side - heatmap_scale_width + (i * heatmap_scale_width / len(colors))
        dwg.add(dwg.rect(insert=(x, margin_top / 2), size=(heatmap_scale_width / len(colors), heatmap_scale_height),
                         fill=color, stroke='none'))
    dwg.add(dwg.text('Min', insert=(svg_width - margin_side - heatmap_scale_width, margin_top / 2 + heatmap_scale_height + 15),
                     font_size=20, text_anchor='start'))
    dwg.add(dwg.text('Max', insert=(svg_width - margin_side, margin_top / 2 + heatmap_scale_height + 15),
                     font_size=20, text_anchor='end'))

    dwg.save()

create_visualization('AllChrGenes2.cir', 'CoverageCa2.cir', 'GenesCoverage_visualization.svg')

