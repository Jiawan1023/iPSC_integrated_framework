while read LINE; do
    sample_id=`echo $LINE | cut -d" " -f1`
    echo "${sample_id}"
    outdir="/path/data/filter_for_sv/${sample_id}"
    mkdir ${outdir}
    
    # Filter out common gnomad variants
    /data7/software/bedtools/bedtools2/bin/intersectBed -a /path/data/calls_with_gene_names_bed/${sample_id}.bed -b /path/gnomad_common.bed -f 0.5 -r -c > ${outdir}/${sample_id}_remove_common_gnomad.bed
    
    # Filter for overlap with gene, exon, intron, utr3, utr5, upstream, & downstream
    /data7/software/bedtools/bedtools2/bin/intersectBed -a /path/data/calls_with_gene_names_bed/${sample_id}.bed -b /path/data/annotations/use_these/*.bed -r -wa -wb -filenames > ${outdir}/${sample_id}_gene_and_stream_overlaps.bed
    
    # Replace file paths with the actual name (i.e. ".../exon.bed" to "exon")
    sed -i 's#/gene.bed#gene#g' ${outdir}/${sample_id}_gene_and_stream_overlaps.bed
    sed -i 's#/exon.bed#exon#g' ${outdir}/${sample_id}_gene_and_stream_overlaps.bed
    sed -i 's#/intron.bed#intron#g' ${outdir}/${sample_id}_gene_and_stream_overlaps.bed
    sed -i 's#/utr3.bed#utr3#g' ${outdir}/${sample_id}_gene_and_stream_overlaps.bed
    sed -i 's#/utr5.bed#utr5#g' ${outdir}/${sample_id}_gene_and_stream_overlaps.bed
    sed -i 's#/upstream.bed#upstream#g' ${outdir}/${sample_id}_gene_and_stream_overlaps.bed
    sed -i 's#/downstream.bed#downstream#g' ${outdir}/${sample_id}_gene_and_stream_overlaps.bed


    # TODO: add append and sort values script
done < /slurm/files/6_smap.txt
