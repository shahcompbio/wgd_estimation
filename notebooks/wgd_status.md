---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Globals

```python
import pandas as pd
import numpy as np
import yaml
import glob
import os
import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 13})
```

# Plot WGD status

```python
def calculate_wgd_status(df:pd.DataFrame) -> pd.Series:
    """ return WGD status based on LOH (fraction) and Ploidy
    """
    wgd_status = ((2.8 - (1.1/0.85) * df['LOH']) < df['Ploidy'])# -(1.1/0.85) * LOH + 2.8 < Ploidy
    return wgd_status
```

## parse WGD status tsv

```python
metadata_path = yaml.load(open('../config/pipeline.yaml'), Loader=yaml.Loader)['metadata']
metadata = pd.read_table(metadata_path)
samples = metadata['isabl_sample_id'].unique()
```

```python
results_path = '../output/results'
```

```python
df = pd.DataFrame()
for sample in samples:
    wgd_path = f'{results_path}/{sample}.wgd.tsv'
    if not os.path.exists(wgd_path): continue
    wgd = pd.read_table(wgd_path)
    wgd.index = [sample]
    df = pd.concat([df, wgd])
```

## parse PCAWG

```python
pcawg_path = '/juno/work/shah/users/mcphera1/repos/spectrumanalysis/external/pcawg_cnv/ploidy_loh.csv.gz'
pcawg = pd.read_csv(pcawg_path)
pcawg.rename(columns={'loh': 'LOH', 'ploidy':'Ploidy'}, inplace=True)
pcawg.set_index('sample', inplace=True)
wgd_status = calculate_wgd_status(pcawg)
```

## QC

```python
# wgd_status = ((2.9 - 2 * df['LOH']) <= df['Ploidy'])# 2.9 -2*hom <= ploidy
# wgd_status = ((2.8 - (1.1/0.85) * df['LOH']) < df['Ploidy'])# -(1.1/0.85) * LOH + 2.8 < Ploidy
_wgd_status = calculate_wgd_status(df)
```

```python
(_wgd_status != df['WGD']).sum()
```

## plot

```python
# line_x = np.linspace(0.02, 0.78, 50)
line_x = np.linspace(0.01, 0.95, 50)
```

```python
line_y = 2.8 - (1.1/0.85) * line_x
```

```python
?leg.handles
```

```python
cohort = 'DLBCL'

fig, ax = plt.subplots(1, 1)
fig.set_figheight(5)
fig.set_figwidth(5)
p1 = sns.scatterplot(data=pcawg, x='LOH', y='Ploidy', color='grey', alpha=0.3, ax=ax, s=10)
p2 = sns.scatterplot(data=df, x='LOH', y='Ploidy', hue='WGD', ax=ax, alpha=0.7, s=30)
ax.plot(line_x, line_y, "k:")#color='grey', linestyle='--')
ax.set_xlabel('LOH');
ax.set_title(f'{cohort} (n={df.shape[0]})');

leg = ax.legend(title="WGD", frameon=True)
ax.add_artist(leg)
h = [plt.plot([],[], color="gray", marker="o", ms=2, ls="")[0]]
plt.legend(handles=h, labels=[''], title="PCAWG", loc=(0.695,0.6), frameon=True)
```
