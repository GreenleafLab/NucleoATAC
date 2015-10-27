###Frequently Asked Questions



**1) Can I use NucleoATAC on my single-end ATAC-seq data?**

No. The method relies on the fragment size information for each sequenced fragment. The functions included with pyatac that take in bam files are also intended exclusively for paired-end data.

**2)  My fragment size distribution shows no evidence of nucleosome laddering-- there is no peak at ~150-200 bp.  Can I still use NucleoATAC?**

No, the method relies on nucleosome-spanning fragments being present in the ATAC-seq libraries.

**3) How can I interpret the \*.nucleoatac_signal.bedgraph.gz track?  Why are there regions with negative values?  What does height of signal track mean?**

The height of this signal track corresponds to how well the ATAC-seq fragmentation pattern at that location matches the expected V-plot pattern for a nucleosome relative to a background model based on sequence bias.  Negative values represent when the signal is less like V-plot than might be expected from the background model.  The background model is somewhat conservative as the fragment size distribution is based on the distribution of all fragments in the peaks supplied.  **It is important not to interpret low or even negative values in this track as indicating the lack of a nucleosome!**  The height of peaks is very much dependent on the number of fragments mapping to a location, which will be dependent on the accessibility of the linker regions.


**4) How can I interpret the \*.occ.bedgraph.gz track?**  

This track is not based on the V-plot pattern, but just the fragment size distribution at a locus.  It thus does not provide as high-resolution positioning information as the \*.nucleoatac_signal.bedgraph.gz track.  However, as the nucleosome occupancy estimate is based on the proportion of fragments estimated to be coming from the nucleosomal distribution, the height of peaks will not be confounded with accessibility like the nucleoatac_signal.  This track is thus much better for comparing nucleosome occupancy between loci or between samples.  

When the signal-to-noise ratio in a sample is low, the nucleosome occupancy estimate will be systematically biased downwards.  This will especially be true outside of accessible regions of the genome, as most fragments in such region will be noise.  Thus to confidently determine that a region is nucleosome free, the occupancy track should be considered along with some measure of coverage, such as a track of the insertion density.  


**5) There are several outputs with positions of nucleosomes-- what are the differences?**

The \*.occpeaks.bed.gz file contains positions of nucleosomes estimated from the \*.occ.bedgraph.gz track.  These will be *low resolution*.

The \*.nucpos.bed.gz file contains positions of nucleosomes estimated from the \*.nucleoatac_signal.bedgraph.gz track.  These will be higher resolution, and should enable observation of underlying sequence preferences of nucleosomes.

The \*.nucmap_combined.bed.gz file is a merge of the other two, with the positions from the nucleoatac_signal favored when there is an overlap.  Generally many more positions will be called from the occupancy track, so this combined track is a way to get a more comprehensive set of positions while still having the higher resolution when possible.

**6)  What is the \*.nfrpos.bed.gz output?**

These are regions between two nucleosomes that have low intervening nucleosome occupancy.











 




