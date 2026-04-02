#!/bin/bash

# Download low confidence regions from UCSC genome browser
CENTROMERES=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
SDS=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz


mkdir -p files/UCSC_files
cd files/UCSC_files

wget $CENTROMERES
wget $SDS

# Other regions were downloaded from the TableBrowser (https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=3474400103_iKAqa0pIoIOqiPqrGOprSGNCHo3p&c=chr16&g=problematic)
