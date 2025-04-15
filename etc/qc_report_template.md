## QC Information for {{ conf['project']['project_id'] }}

### Configuration and Summary Statistics

Config:

min_reads = {{ conf['fastq']['min_reads'] }}  
max_repeats = {{ conf['fastq']['max_repeats'] }}  
max_mismatch = {{ conf['collapse']['max_mismatch'] }}  
target_min_reads = {{ conf['vbctable']['target_min_reads'] }}  
inj_min_reads = {{ conf['vbctable']['inj_min_reads'] }}  
target_min_umi = {{  conf['vbctable']['target_min_umi'] }} 

FASTQ Processing:

{% for ( k,v ) in stats['fastq'].items() %}
  {{ k }}  = {{ v }}  
{% endfor %}

Split and Filter:

{% for ( k,v ) in stats['fastq_filter'].items() %}
  {{ k }}  = {{ v }}  
{% endfor %}

Align and Collapse:

{% for ( k,v ) in stats['collapse'].items() %}
  {{ k }}  = {{ v }}  
{% endfor %}

Read Table:

{% for ( k,v ) in stats['readtable'].items() %}
  {{ k }}  = {{ v }}  
{% endfor %}

VBC Table:

{% for ( k,v ) in stats['vbctable'].items() %}
  {{ k }}  = {{ v }}  
{% endfor %}

{% for bid, bdict in stats['matrices'].items() %}
### Brain {{ bid }}

{% for (k,v) in bdict.items() %}
{{ k }}  = {{ v }}  
{% endfor %}

{% endfor %}