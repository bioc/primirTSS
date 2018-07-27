# primirTSS
primirTSS is a first R package to predict **pri-miRNA** Transcription Start Site.

Identifying human miRNA transcriptional start sites (TSSs) plays a significant role in understanding the transcriptional regulation of miRNA. However, **due to the quick capping of pri-miRNA** and many miRNA genes may lie in the introns or even exons of other genes, it is difficult to detect miRNA TSSs. miRNA TSSs are cell-specific. And miRNA TSSs are cell-specific, which implies the same miRNA in different cell-lines may start transcribing at different TSSs.

High throughput sequencing, like ChIP-seq, has gradually become an essential and versatile approach for us to identify and understand genomes and their transcriptional processes. 
By integrating H3k4me4 and Pol II data, parting of false positive counts after scoring can be filter out. Besides, DNase I hypersensitive sites(DHS) also imply TSSs, where miRNAs will be accessible and functionally related to transcription activities. And additionally, the expression profile of miRNA and genes in certain cell-line will be considered as well for improve fidelity. By employing all these different kinds of data, here we have developed the primirTSS package to assist users to identify miRNA TSSs in human and to provide them with related information about where miRNA genes lie in the genome, with both **command-line** and **graphical** interfaces.

Before using primirTSS, you should install **java** first. Users can learn  the details of each function in the package vignettes and the help pages.
