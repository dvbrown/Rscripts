calc_seqSaturation = function(singlecellexperiment) {
    require(scPipe)
    # Sequencing Saturation = 1 - (n_deduped_reads / n_reads)
    # n_deduped_reads = Number of unique (valid cell-barcode, valid UMI, gene) combinations among confidently mapped reads. 
    # n_reads = Total number of confidently mapped, valid cell-barcode, valid UMI reads.
    dup_info = (UMI_dup_info(sce))
    n_reads = sum(dup_info$count)
    non_zero = dup_info[dup_info$count > 0,]
    n_deduped_reads = dup_info[1,2]
    saturation = 1 - (n_deduped_reads / n_reads)
    saturation_percent = saturation * 100
    return (saturation_percent)
}

calc_numReads_newUMI = function(singlecellexperiment) {
    require(scPipe)
    # The inverse of the sequencing saturation can be interpreted as roughly the number of new transcripts 
    # you expect to find with one new read. If sequencing saturation is at 50%, it means that every 2 new reads
    #will result in 1 new UMI count (unique transcript) detected. 
    # In contrast, 90% sequencing saturation means that 10 new reads are necessary to obtain one new UMI count.
    dup_info = (UMI_dup_info(sce))
    n_reads = sum(dup_info$count)
    non_zero = dup_info[dup_info$count > 0,]
    n_deduped_reads = dup_info[1,2]
    unsaturation = n_deduped_reads / n_reads
    newReads = unsaturation ^-1
    return (newReads)
}