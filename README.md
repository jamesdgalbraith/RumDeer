# RumDeer
Guide to calling variants from whole genome sequencing using the Rum Deer dataset as an example. The three steps outlined are what I used for the Soay sheep, reindeer, and sparrows.

## 01 Pseudo-autosomal region checks
In chromosomal assemblies which contain both sex chromosomes it's necessary to check if the pseudoautosomal region (PAR) is present on both sex chromosomes, and where the boundaries of the PAR are. If it is present on both sex chromosomes it should be masked on heteromorphic chromosome before reads are aligned to prevent multimapping. The red deer assembly doesn't contain the Y chromosome, so this isn't necessary here. However, independant of whether the PAR is on one or both sex chromosomes (or only the homomorphic sex chromosome is included), the variant caller needs to be told where the PAR's boundary are to ensure the remainder of the sex chromosomes can be properly called as diploid, haploid or absent.

As reference genome papers might not include details about the PAR, we'd first check if it's on both sex chromosomes by aligning them to each other using [D-GENIES](https://dgenies.toulouse.inra.fr/) (run locally). The reindeer and sheep variant calling has more detailed descriptions of this. In the sheep both sex chromosomes contained the PAR, while in the reindeer only the X did.

If the PAR's not present on both sex chromosomes, it's normally present on the Z/X. Examining read coverage/depth across the Z/X in a male/female (respecitvely) enables us to determine where the PAR begins/ends. In the case of the red deer, the reference genome doesn't contain the Y chromosome, so this step will also help us determine if the PAR is on the X (I expect it will be).

We'll align the read data for a male (seeing as deer are XY) to the reference genome, calculate read depth across the smallest autosome and X chromosome with `samtools depth`, and examine where coverage drops to half on the X chromosome.

## 02 Read alignment
