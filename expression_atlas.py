#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Libraries
import argparse
import sys
import os
import json

# Delimiter
def detect_delimiter(line):
    if "\t" in line and line.count("\t") >= line.count(","):
        return "\t"
    elif ";" in line:
        return ";"
    elif "," in line:
        return ","
    else:
        return "\t"

# Function
def main():
    parser = argparse.ArgumentParser(description="Generate an HTML Gene Expression Atlas report with interactive charts.")
    parser.add_argument("--matrix", required=True,help="Normalized expression matrix file (TXT/TSV/CSV with header).")
    parser.add_argument("--genes", required=True,help="File containing a list of genes of interest (one ID per line, TXT with no header).")
    parser.add_argument("--annot", required=False,help="Optional annotations file (TXT/TSV/CSV with header).")
    parser.add_argument("--output", required=True,help="Name of the output HTML report file.")
    args = parser.parse_args()

    matrix_file = args.matrix
    genes_file = args.genes
    annot_file = args.annot
    output_file = args.output

    # Verify files
    if not os.path.isfile(matrix_file):
        print(f"Error: File not found: '{matrix_file}'", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.isfile(genes_file):
        print(f"Error: File not found: '{genes_file}'", file=sys.stderr)
        sys.exit(1)
    
    if annot_file and not os.path.isfile(annot_file):
        print(f"Warning: Annotation file '{annot_file}' not found; proceeding without it", file=sys.stderr)
        annot_file = None

    # Read genes
    with open(genes_file, encoding='utf-8') as f:
        gene_list = [l.strip() for l in f if l.strip()]
    gene_set = set(gene_list)

    # Read matrix
    group_to_indices = {}
    group_order = []
    expression_data = {}
    with open(matrix_file, encoding='utf-8') as mf:
        header = mf.readline().strip()
        delim = detect_delimiter(header)
        cols = header.split(delim)
        samples = cols[1:]
        for i, s in enumerate(samples):
            if s not in group_to_indices:
                group_to_indices[s] = []
                group_order.append(s)
            group_to_indices[s].append(i)
        for line in mf:
            parts = line.strip().split(delim)
            gene = parts[0]
            if gene in gene_set:
                vals = []
                for x in parts[1:]:
                    try:
                        vals.append(float(x))
                    except:
                        vals.append(None)
                expression_data[gene] = vals
                if len(expression_data) == len(gene_set):
                    break

    # Read annotations
    annot_headers = []
    annot_data = {}
    gene_col = 0
    if annot_file:
        with open(annot_file, encoding='utf-8') as af:
            first = af.readline().strip()
            sep = detect_delimiter(first)
            hdr = first.split(sep)
            annot_headers = hdr
            for idx, h in enumerate(annot_headers):
                if 'gene' in h.lower() and 'id' in h.lower():
                    gene_col = idx
                    break
            for line in af:
                parts = line.strip().split(detect_delimiter(line))
                if len(parts) > gene_col:
                    annot_data[parts[gene_col]] = parts

    # Build series
    line_series = []
    heatmap_series = []
    for g in gene_list:
        if g not in expression_data:
            continue
        vals = expression_data[g]
        means = []
        pts = []
        for s in group_order:
            idxs = group_to_indices[s]
            v = [vals[i] for i in idxs if i < len(vals) and vals[i] is not None]
            m = round(sum(v) / len(v), 2) if v else None
            means.append(m)
            pts.append({'x': s, 'y': m})
        line_series.append({'name': g, 'data': means})
        heatmap_series.append({'name': g, 'data': pts})
    heatmap_series.reverse()

    # Build replicates series
    max_reps = max(len(v) for v in group_to_indices.values())
    dot_data = {}
    for g in gene_list:
        if g not in expression_data:
            continue
        vals = expression_data[g]
        ser = []
        for r in range(max_reps):
            dpts = []
            for s in group_order:
                idxs = group_to_indices[s]
                if r < len(idxs) and idxs[r] < len(vals) and vals[idxs[r]] is not None:
                    dpts.append({'x': s, 'y': round(vals[idxs[r]], 2)})
            if dpts:
                ser.append({'name': f'Rep {r+1}', 'data': dpts})
        dot_data[g] = ser

    # Heatmap color ranges
    color_ranges = [
        {"from": 0,     "to": 0.99,    "color": "#eceff1", "name": "<1"},
        {"from": 1,     "to": 1.99,    "color": "#b3e5fc", "name": ">=1"},
        {"from": 2,     "to": 4.99,    "color": "#80cbc4", "name": ">=2"},
        {"from": 5,     "to": 9.99,    "color": "#ffee58", "name": ">=5"},
        {"from": 10,    "to": 49.99,   "color": "#ffb74d", "name": ">=10"},
        {"from": 50,    "to": 99.99,   "color": "#ff8f00", "name": ">=50"},
        {"from": 100,   "to": 199.99,  "color": "#ff4f00", "name": ">=100"},
        {"from": 200,   "to": 499.99,  "color": "#cc0000", "name": ">=200"},
        {"from": 500,   "to": 999.99,  "color": "#D72C79", "name": ">=500"},
        {"from": 1000,  "to": 4999.99, "color": "#801C5A", "name": ">=1000"},
        {"from": 5000,  "to": 50000,   "color": "#6D3917", "name": ">=5000"}
    ]

    # Write HTML
    with open(output_file, 'w', encoding='utf-8') as out:

        # Header
        out.write('<!DOCTYPE html>\n<html>\n<head>\n<meta charset="utf-8">\n')

        out.write('<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">\n')
        out.write('<link rel="stylesheet" href="https://cdn.datatables.net/1.10.24/css/dataTables.bootstrap4.min.css">\n')
        out.write('<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>\n')
        out.write('<script src="https://cdn.datatables.net/1.10.24/js/jquery.dataTables.min.js"></script>\n')
        out.write('<script src="https://cdn.datatables.net/1.10.24/js/dataTables.bootstrap4.min.js"></script>\n')
        out.write('<script src="https://cdnjs.cloudflare.com/ajax/libs/apexcharts/4.4.0/apexcharts.min.js"></script>\n')

        out.write('</head>\n<body>\n  <div class="container">\n')

        # Title
        out.write('<h1 class="mt-4 mb-4">Gene Expression Atlas</h1>')

        # Line Chart
        out.write('<div class="card mb-4">\n  <div class="card-body">\n')
        out.write('<h5 class="card-title">Line Chart</h5>\n')
        out.write('<div id="lineChart"></div>\n  </div>\n</div>\n')

        # Bar Chart
        out.write('<div class="card mb-4">\n  <div class="card-body">\n')
        out.write('<h5 class="card-title">Bar Chart</h5>\n')
        out.write('<div id="barChart"></div>\n  </div>\n</div>\n')

        # Heatmap
        out.write('<div class="card mb-4">\n  <div class="card-body">\n')
        out.write('<h5 class="card-title">Heatmap</h5>\n')
        out.write('<div id="heatmapChart"></div>\n  </div>\n</div>\n')

        # Replicates Chart
        out.write('<div class="card mb-4">\n  <div class="card-body">\n')
        out.write('<h5 class="card-title">Replicates Chart</h5>\n')
        out.write('<label for="geneSelect">Select gene:</label><br/>\n')
        out.write('<select id="geneSelect" class="form-control mb-3">\n')
        for g in gene_list:
            out.write(f'  <option value="{g}">{g}</option>\n')
        out.write('</select>\n')
        out.write('<div id="dotChart"></div>\n  </div>\n</div>\n')

        # Datatable
        out.write('<div class="card mb-4">\n  <div class="card-body">\n')
        out.write('<h5 class="card-title">Expression Data Table</h5>\n')
        out.write('<div class="table-responsive">\n')
        out.write('<table id="datatable" class="display table table-striped table-bordered" border="0">\n')
        out.write('<thead><tr>\n')
        out.write('<th>Gene ID</th>\n')
        for s in group_order:
            for r in range(len(group_to_indices[s])):
                out.write(f'<th>{s} ({r+1})</th>')
        if annot_headers:
            for i, c in enumerate(annot_headers):
                if i == gene_col:
                    continue
                out.write(f'<th>{c}</th>')
        out.write('</tr></thead><tbody>\n')
        for g in gene_list:
            if g not in expression_data:
                continue
            out.write('<tr>')
            out.write(f'<td>{g}</td>')
            for s in group_order:
                for idx in group_to_indices[s]:
                    val = expression_data[g][idx] if idx < len(expression_data[g]) else None
                    cell = round(val, 2) if val is not None else ''
                    out.write(f'<td>{cell}</td>')
            if annot_headers:
                ann = annot_data.get(g, [])
                for i, c in enumerate(annot_headers):
                    if i == gene_col: continue
                    v = ann[i] if i < len(ann) else ''
                    out.write(f'<td>{v}</td>')
            out.write('</tr>\n')
        out.write('</tbody></table>\n')
        out.write('</div></div>\n')

        # Scripts
        out.write('<script>\n')
        out.write(f'var lineSeries={json.dumps(line_series)};\n')
        out.write(f'var barSeries={json.dumps(line_series)};\n')
        out.write(f'var heatmapSeries={json.dumps(heatmap_series)};\n')
        out.write(f'var dotData={json.dumps(dot_data)};\n')
        out.write(f'var categories={json.dumps(group_order)};\n')

        # Line chart
        out.write("""
            var lineOptions={
                chart:{type:'line',height:400},
                series:lineSeries,
                xaxis:{categories:categories},
                dataLabels:{enabled:true},
                legend:{position:'top'}
            };
            var lineChart=new ApexCharts(document.querySelector('#lineChart'),lineOptions);
            lineChart.render();
        """)

        # Bar chart
        out.write("""
            var barOptions={
                chart:{type:'bar',height:400},
                series:barSeries,
                xaxis:{categories:categories},
                plotOptions:{bar:{dataLabels:{position:'top'}}},
                dataLabels:{enabled:true,style:{colors:['#fff']},background:{enabled:true,foreColor:'#000'}},
                legend:{position:'top'}
            };
            var barChart=new ApexCharts(document.querySelector('#barChart'),barOptions);
            barChart.render();
        """)

        # Heatmap
        out.write(f"""
            var heatmapOptions={{
                chart:{{type:'heatmap',height:{max(400,30*len(heatmap_series))}}},
                series:heatmapSeries,
                xaxis:{{type:'category'}},
                plotOptions:{{
                    heatmap:{{
                        useFillColorAsStroke:true,
                        colorScale:{{
                            ranges:{json.dumps(color_ranges)}
                        }}
                    }}
                }},
                dataLabels:{{enabled:true,style:{{colors:['#fff']}}}}
            }};
            var heatmapChart=new ApexCharts(document.querySelector('#heatmapChart'),heatmapOptions);
            heatmapChart.render();
        """)

        # Replicates chart
        first = gene_list[0] if gene_list else ''
        out.write(f"""
            var scatterOptions={{
                chart:{{type:'scatter',height:400}},
                series:dotData['{first}'],
                xaxis:{{type:'category',categories:categories}},
                markers:{{size:6}},
                legend:{{position:'top'}}
            }};
            var dotChart=new ApexCharts(document.querySelector('#dotChart'),scatterOptions);
            dotChart.render();
            var sel=document.querySelector('#geneSelect');
            sel.addEventListener('change',function(){{if(dotData[this.value])dotChart.updateSeries(dotData[this.value]);}});
        """)

        # Datatable
        out.write('$(function(){$("#datatable").DataTable();});\n')

        # End function
        out.write('</script>\n</div>\n</body>\n</html>')


# Call the function
if __name__ == '__main__':
    main()