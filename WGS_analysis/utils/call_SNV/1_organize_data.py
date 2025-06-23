import pandas as pd

import os
from pathlib import Path

import gzip
import io

# Create a list of files
# Get FASTQ for the iPSC project
data_locs=['/data6/wgs_16p_novogene/data/usftp21.novogene.com/01.RawData',
			'/data6/wgs_16p_novogene/data/download_03_18_2024/usftp21.novogene.com/01.RawData',
			'/data6/wgs_16p_novogene/data/download_06_01_2025/usftp21.novogene.com/01.RawData']

print('Finding files')
files=[]
for dl in data_locs:
	for path in Path(dl).rglob('*1.f*q.gz'):
		# Check the size of the file, sometimes samples are duplicated without actual reads
		if os.path.getsize(str(path.parent)+'/'+str(path.name))>5000000000:
			files.append([str(path.name), str(path.parent)+'/'])

df=pd.DataFrame(files, columns=['read1', 'file_location'])

# Restrict to only needed samples
df['Sample']=df.read1.str.split('_', expand=True)[0]
df=df[df.Sample.str.contains('C')]

# Add in read 2
df['read2']=df.read1.str.replace('1.', '2.', regex=False)

# Organize
df=df[['Sample', 'file_location', 'read1', 'read2']]

print('Add in sequencing information')

print('Read group')
# Add read group (no sample multiplexing, so each sample is a separate read group)
df['read_group']=df.Sample+'_RG'

print('Library')
# Library information (unique for each sample)
df['library_name']='Macrogen-'+df.Sample
df.loc[df.file_location.str.contains('novogene'), 'library_name']='Novogene-'+df[df.file_location.str.contains('novogene')]['Sample']

print('Platform')
# Platform unit (taken from read header)
platforms=[]
for idx, row in df.iterrows():
	filename=row['file_location']+row['read1']
	with gzip.open(filename, 'rb') as ip:
		with io.TextIOWrapper(ip, encoding='utf-8') as decoder:
			header=decoder.readline()

			ins_name=header.split(':')[0]
			ins_name=ins_name.replace('@', '')
			platforms.append(ins_name)

df['platform_unit']=platforms

# Platform name
df['platform_name']='illumina'

print('Sequencing center')
# Sequencing center
df['sequencing_center']='Macrogen'
df.loc[df.file_location.str.contains('novogene'), 'sequencing_center']='Novogene'

print(df)

# Save all info to file
df.to_csv('files/sample_info.csv', index=False)
