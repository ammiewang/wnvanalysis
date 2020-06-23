# wnvanalysis
These programs were used to analyze West Nile Virus effective population/sequence divergence in the Cohan Lab paper "Delving Below the Species Level To Characterize the Ecological Diversity within the Global Virome: An Exploration of West Nile Virus".

Requirements
- python3
- BioPython
- ete3

Recommended
- pypy3

Steps:
0. Run ES2.
1. Run draw_lines.py on the log file from ES2.
2. Run ecotype_files.py on the outputs from draw_lines.py.
3. Run pd_overall.py or pd_overall_with_seqdiv.py on the outputs from ecotype_files.py.
