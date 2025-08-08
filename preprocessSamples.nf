process makeBamToSampleNameMap {
    input:
    path(bam_files)

    output:
    path("bam_to_sample_name_map.txt")

    script:
    """
    for f in ${bam_files.join(' ')}; do
        samplename=\$(samtools view -H "\$f" | awk '
    {
      for (i=1; i<=NF; i++)
        if (\$i ~ /^SM:/)
        {
            split(\$i, a, ":");
            print a[2];
            exit
        }
    }')
        filename=\$(basename "\$f")
        filename=\${filename%.*}
        printf "%s\t%s\n" \$filename \$samplename
    done >> bam_to_sample_name_map.txt
    """
}
