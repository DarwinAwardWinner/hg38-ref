rule rulegraphs:
    input: 'rulegraph-all.png', 'dag-all.png'

rule svg_to_png:
    input: '{filename}.svg'
    output: '{filename}.png'
    shell: '''rsvg-convert -d 180 -p 180 -f png {input:q} -o {output:q}'''

rule dag_svg:
    output: 'dag-{target}.svg'
    shell: '''snakemake --nolock -f --dag {wildcards.target:q} | dot -Tsvg > {output:q}'''

rule rulegraph_svg:
    output: 'rulegraph-{target}.svg'
    shell: '''snakemake --nolock -f --rulegraph {wildcards.target:q} | dot -Tsvg > {output:q}'''
